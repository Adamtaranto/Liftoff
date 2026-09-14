"""FASTA helpers: extracting gene sequences and splitting target genomes."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from pyfaidx import Fasta, FastaIndexingError, UnsupportedCompressionFormat

from liftoff.errors import InputError
from liftoff.models import Feature, LiftoverType


def open_fasta(path: str, *, first_word_keys: bool = False) -> Any:
    """Open an indexed FASTA file, creating the ``.fai`` index if needed.

    Parameters
    ----------
    path : str
        FASTA file (plain or BGZF-compressed).
    first_word_keys : bool, default False
        Key sequences by the first whitespace-delimited word of the header.

    Returns
    -------
    pyfaidx.Fasta
        Random-access FASTA reader.

    Raises
    ------
    InputError
        If the file does not exist, is compressed with plain gzip instead of
        BGZF, or cannot be indexed.
    """
    if not Path(path).is_file():
        raise InputError(f'FASTA file not found: {path}')
    try:
        if first_word_keys:
            return Fasta(path, key_function=lambda name: name.split()[0])
        return Fasta(path)
    except UnsupportedCompressionFormat as exc:
        raise InputError(
            f'{path} is gzip-compressed; compressed FASTA must use BGZF '
            '(e.g. `gunzip file.fa.gz && bgzip file.fa`) or be decompressed'
        ) from exc
    except FastaIndexingError as exc:
        raise InputError(f'could not index FASTA file {path}: {exc}') from exc


def split_target_sequence(
    target_chroms: Sequence[str], target_fasta_name: str, inter_files: str | Path
) -> Any:
    """Write each requested target chromosome to its own FASTA file.

    Parameters
    ----------
    target_chroms : sequence of str
        Chromosome names to extract. An entry equal to ``target_fasta_name``
        stands for the whole genome and is not written.
    target_fasta_name : str
        Target genome FASTA.
    inter_files : str or pathlib.Path
        Directory receiving ``<chrom>.fa`` files.

    Returns
    -------
    pyfaidx.Fasta
        Reader for the whole target genome.
    """
    # Opening the FASTA builds (or refreshes) its .fai index.
    target_fasta = open_fasta(target_fasta_name, first_word_keys=True)
    for chrom in target_chroms:
        if chrom != target_fasta_name:
            with (Path(inter_files) / f'{chrom}.fa').open('w') as out:
                out.write('>' + chrom + '\n' + str(target_fasta[chrom]))
    return target_fasta


def get_genome_size(fasta: Any) -> int:
    """Return the total length of all sequences in a FASTA file.

    Parameters
    ----------
    fasta : pyfaidx.Fasta
        Indexed FASTA reader.

    Returns
    -------
    int
        Sum of sequence lengths.
    """
    return sum(len(record) for record in fasta.values())


def features_file_name(
    chrom_name: str, reference_fasta_name: str, liftover_type: LiftoverType
) -> str:
    """Return the base name used for extracted gene sequences of a stage.

    Parameters
    ----------
    chrom_name : str
        Reference chromosome being processed, or the reference FASTA path
        when processing the whole genome.
    reference_fasta_name : str
        Reference genome FASTA path.
    liftover_type : LiftoverType
        Pipeline stage.

    Returns
    -------
    str
        Name such as ``"reference_all"`` or ``"<chrom>"``.
    """
    if chrom_name == reference_fasta_name and liftover_type == LiftoverType.CHROM_BY_CHROM:
        return 'reference_all'
    if chrom_name == reference_fasta_name and liftover_type == LiftoverType.COPIES:
        # Distinct from the first stage so each stage's SAM file has a unique
        # name (required to look up pre-computed alignments).
        return 'reference_all_copies'
    if liftover_type == LiftoverType.UNMAPPED:
        return 'unmapped_to_expected_chrom'
    if liftover_type == LiftoverType.UNPLACED:
        return 'unplaced'
    return chrom_name


def write_gene_sequences(
    parent_dict: Mapping[str, Feature],
    ref_chroms: Sequence[str],
    reference: str,
    inter_files: str | Path,
    liftover_type: LiftoverType,
    flank: float,
) -> None:
    """Extract the sequence of every top-level feature to lift.

    The coordinates of each feature are set **in place** to the extracted
    region: the annotated coordinates extended by ``flank`` times the feature
    length on each side (clipped to the chromosome), so that later coordinate
    conversions refer to the extracted sequence. The annotated coordinates are
    kept in :attr:`~liftoff.models.Feature.original_bounds`, and every call
    extends from them, so repeated stages do not widen the flank again.

    Parameters
    ----------
    parent_dict : mapping of str to Feature
        Features whose sequences are extracted.
    ref_chroms : sequence of str
        Reference chromosomes to process; the reference FASTA path stands for
        the whole genome.
    reference : str
        Reference genome FASTA.
    inter_files : str or pathlib.Path
        Directory receiving ``<name>_genes.fa`` files.
    liftover_type : LiftoverType
        Pipeline stage (controls file naming and append mode).
    flank : float
        Fraction of feature length to add on each side.

    Raises
    ------
    InputError
        If there are no features to extract.
    """
    inter_files = Path(inter_files)
    fasta = open_fasta(reference)
    if liftover_type == LiftoverType.UNPLACED:
        # Unplaced sequences are appended to one shared file; truncate it first.
        (inter_files / 'unplaced_genes.fa').open('w').close()
    sorted_parents = sorted(parent_dict.values(), key=lambda feature: feature.seqid)
    if not sorted_parents:
        raise InputError(
            'no top-level features were selected for lift-over: the annotation contains no '
            'features of the selected types (genes by default). Use --feature-types or '
            '--all-feature-types to select other types, or check --exclude-feature-types.'
        )
    mode = 'a' if liftover_type == LiftoverType.UNPLACED else 'w'
    for chrom in ref_chroms:
        name = features_file_name(chrom, reference, liftover_type)
        with (inter_files / f'{name}_genes.fa').open(mode) as fasta_out:
            _write_chrom_gene_sequences(chrom, reference, fasta, sorted_parents, fasta_out, flank)


def _write_chrom_gene_sequences(
    chrom_name: str,
    reference_fasta_name: str,
    fasta: Any,
    parents: Sequence[Feature],
    fasta_out: Any,
    flank: float,
) -> None:
    """Write sequences of ``parents`` located on ``chrom_name`` (or all)."""
    whole_genome = chrom_name == reference_fasta_name
    current_chrom = parents[0].seqid if whole_genome else chrom_name
    chrom_seq: str = fasta[current_chrom][:].seq
    for parent in parents:
        if whole_genome and parent.seqid != current_chrom:
            # Parents are sorted by sequence, so load each chromosome once.
            current_chrom = parent.seqid
            chrom_seq = fasta[current_chrom][:].seq
        if parent.seqid == chrom_name or whole_genome:
            if parent.original_bounds is None:
                parent.original_bounds = (parent.start, parent.end)
            start, end = parent.original_bounds
            gene_length = end - start + 1
            parent.start = round(max(1, start - flank * gene_length))
            parent.end = round(min(end + flank * gene_length, len(chrom_seq)))
            parent_seq = chrom_seq[parent.start - 1 : parent.end]
            fasta_out.write('>' + parent.id + '\n' + parent_seq + '\n')
