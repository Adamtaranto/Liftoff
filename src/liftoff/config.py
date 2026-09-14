"""Typed run configuration.

:class:`LiftoffConfig` replaces the ``argparse.Namespace`` that earlier versions
threaded through (and mutated inside) every function. It is populated by the
command line interface but can equally be constructed directly when using
Liftoff as a library.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
import shlex

from liftoff.errors import ConfigError

#: minimap2 options Liftoff relies on. Each entry is ``(flag, value)``; a flag
#: is appended with its default value when the user did not set it.
REQUIRED_MM2_OPTIONS: tuple[tuple[str, str | None], ...] = (
    ('-a', None),  # SAM output
    ('--end-bonus', '5'),  # favour end-to-end gene alignments
    ('--eqx', None),  # =/X CIGAR operations distinguish matches and mismatches
    ('-N', '50'),  # retain many secondary alignments
    ('-p', '0.5'),  # ...even if they score well below the primary
)

#: Default value of ``mm2_options``.
DEFAULT_MM2_OPTIONS = '-a --end-bonus 5 --eqx -N 50 -p 0.5'

#: Value written to indicate that output goes to standard output.
STDOUT = 'stdout'


def merge_minimap2_options(user_options: str) -> list[str]:
    """Combine user-supplied minimap2 options with those Liftoff requires.

    Parameters
    ----------
    user_options : str
        Space-delimited minimap2 arguments. Shell-style quoting is honoured.

    Returns
    -------
    list of str
        Tokenised arguments: the user's options followed by any required
        options that were not already specified.

    Examples
    --------
    >>> merge_minimap2_options('-r 2k -z 5000')
    ['-r', '2k', '-z', '5000', '-a', '--end-bonus', '5', '--eqx', '-N', '50', '-p', '0.5']
    """
    tokens = shlex.split(user_options)
    # Compare whole tokens (not substrings) so that e.g. "--end-bonus" does
    # not satisfy "-N", and "-ax" style combined flags are not mis-detected.
    present = set(tokens)
    for flag, value in REQUIRED_MM2_OPTIONS:
        if flag not in present:
            tokens.append(flag)
            if value is not None:
                tokens.append(value)
    return tokens


@dataclass(slots=True)
class LiftoffConfig:
    """All options controlling a Liftoff run.

    Attributes
    ----------
    target : str
        Target genome FASTA to lift features to.
    reference : str
        Reference genome FASTA to lift features from.
    gff : str or None, default None
        Reference annotation (GFF3 or GTF). Exactly one of ``gff`` and ``db``
        must be given.
    db : str or None, default None
        Pre-built gffutils feature database.
    output : str, default "stdout"
        Output annotation path, or ``"stdout"``.
    unmapped : str, default "unmapped_features.txt"
        File listing IDs of features that could not be mapped.
    exclude_partial : bool, default False
        Report partial/low identity mappings as unmapped instead of writing
        them to the output annotation.
    intermediate_dir : str, default "intermediate_files"
        Directory for intermediate FASTA and SAM files.
    mm2_options : str, default DEFAULT_MM2_OPTIONS
        Additional minimap2 arguments.
    min_coverage : float, default 0.5
        Minimum alignment coverage of child features for a feature to be
        considered mapped.
    min_identity : float, default 0.5
        Minimum sequence identity of child features for a feature to be
        considered mapped.
    distance_factor : float, default 2.0
        Alignment blocks further apart in the target than this factor times
        their distance in the reference are not chained.
    flank : float, default 0.0
        Fraction of gene length added as flanking sequence on each side before
        alignment.
    threads : int, default 1
        Number of parallel alignment processes.
    minimap2 : str or None, default None
        Path to the minimap2 executable; searched on ``PATH`` when omitted.
    feature_types : str or None, default None
        File listing additional top-level feature types to lift.
    all_feature_types : bool, default False
        Lift every top-level feature type in the annotation instead of only
        genes (and ``feature_types``).
    exclude_feature_types : tuple of str, default ()
        Top-level feature types never to lift (e.g. ``'region'``). Applied
        after ``feature_types`` / ``all_feature_types``.
    infer_genes : bool, default False
        Let gffutils infer gene features from transcripts.
    infer_transcripts : bool, default False
        Let gffutils infer transcript features from exons/CDS.
    chroms : str or None, default None
        Comma separated file of corresponding reference,target chromosomes.
    unplaced : str or None, default None
        File of unplaced reference sequence names (requires ``chroms``).
    copies : bool, default False
        Search for additional gene copies in the target.
    copy_identity : float, default 1.0
        Minimum child feature identity for an extra copy.
    max_overlap : float, default 0.1
        Maximum fraction two lifted features may overlap.
    mismatch : int, default 2
        Mismatch penalty used when scoring alignment chains.
    gap_open : int, default 2
        Gap open penalty used when scoring alignment chains.
    gap_extend : int, default 1
        Gap extension penalty used when scoring alignment chains.
    polish : bool, default False
        Re-align exons to repair broken coding sequences.
    cds : bool, default True
        Annotate the ORF status of every lifted CDS.
    alignments : str or None, default None
        Directory of pre-computed SAM files; when given minimap2 is not run.
    command_line : str, default ""
        Command line recorded in the output header.
    """

    target: str
    reference: str
    gff: str | None = None
    db: str | None = None
    output: str = STDOUT
    unmapped: str = 'unmapped_features.txt'
    exclude_partial: bool = False
    intermediate_dir: str = 'intermediate_files'
    mm2_options: str = DEFAULT_MM2_OPTIONS
    min_coverage: float = 0.5
    min_identity: float = 0.5
    distance_factor: float = 2.0
    flank: float = 0.0
    threads: int = 1
    minimap2: str | None = None
    feature_types: str | None = None
    all_feature_types: bool = False
    exclude_feature_types: tuple[str, ...] = ()
    infer_genes: bool = False
    infer_transcripts: bool = False
    chroms: str | None = None
    unplaced: str | None = None
    copies: bool = False
    copy_identity: float = 1.0
    max_overlap: float = 0.1
    mismatch: int = 2
    gap_open: int = 2
    gap_extend: int = 1
    polish: bool = False
    cds: bool = True
    alignments: str | None = None
    command_line: str = field(default='')

    def __post_init__(self) -> None:
        self.validate()

    def validate(self) -> None:
        """Check option combinations and value ranges.

        Raises
        ------
        ConfigError
            If any option is invalid.
        """
        if (self.gff is None) == (self.db is None):
            raise ConfigError('exactly one of a GFF/GTF file or a feature database is required')
        for name in ('min_coverage', 'min_identity', 'copy_identity', 'max_overlap', 'flank'):
            value = getattr(self, name)
            if not 0.0 <= value <= 1.0:
                raise ConfigError(f'{name} must be between 0 and 1 (got {value})')
        if self.min_identity > self.copy_identity:
            raise ConfigError('copy identity must be greater than or equal to minimum identity')
        if self.all_feature_types and self.feature_types is not None:
            raise ConfigError('all feature types cannot be combined with a feature types file')
        if isinstance(self.exclude_feature_types, str):
            raise ConfigError('exclude_feature_types must be a sequence of feature type names')
        if any(not name.strip() or name != name.strip() for name in self.exclude_feature_types):
            raise ConfigError('excluded feature types must be non-empty names without spaces')
        # Normalise to a de-duplicated tuple, preserving the given order.
        self.exclude_feature_types = tuple(dict.fromkeys(self.exclude_feature_types))
        if self.chroms is None and self.unplaced is not None:
            raise ConfigError('unplaced sequences can only be used together with chroms')
        if self.threads < 1:
            raise ConfigError(f'threads must be at least 1 (got {self.threads})')
        if self.distance_factor <= 0:
            raise ConfigError(f'distance factor must be positive (got {self.distance_factor})')
        for name in ('mismatch', 'gap_open', 'gap_extend'):
            if getattr(self, name) < 0:
                raise ConfigError(f'{name} penalty must not be negative')

    @property
    def minimap2_arguments(self) -> list[str]:
        """Tokenised minimap2 options including Liftoff's required flags.

        Returns
        -------
        list of str
            Arguments passed to minimap2.
        """
        return merge_minimap2_options(self.mm2_options)

    @property
    def intermediate_path(self) -> Path:
        """Directory for intermediate files.

        Returns
        -------
        pathlib.Path
            The intermediate directory.
        """
        return Path(self.intermediate_dir)
