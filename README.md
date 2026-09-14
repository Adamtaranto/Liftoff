# Liftoff
![PyPI - Downloads](https://img.shields.io/pypi/dm/liftoff?style=plastic)
[![Conda](https://img.shields.io/conda/dn/bioconda/liftoff?label=bioconda-install&style=plastic)](https://anaconda.org/bioconda/liftoff)
[![CI](https://github.com/agshumate/Liftoff/actions/workflows/CI.yml/badge.svg)](https://github.com/agshumate/Liftoff/actions/workflows/CI.yml)
![Stars](https://img.shields.io/github/stars/agshumate/Liftoff?style=plastic)

<img width="200" src="docs/images/liftoff_logo_current_dark_mode4.svg" alt="Liftoff logo">
Liftoff is a tool that accurately maps annotations in GFF or GTF between assemblies of the same, or closely-related species. Unlike current coordinate lift-over tools which require a pre-generated “chain” file as input, Liftoff is a standalone tool that takes two genome assemblies and a reference annotation as input and outputs an annotation of the target genome. Liftoff uses Minimap2 [(Li, 2018)](https://academic.oup.com/bioinformatics/article/34/18/3094/4994778) to align the gene sequences from a reference genome to the target genome. Rather than aligning whole genomes, aligning only the gene sequences allows genes to be lifted over even if there are many structural differences between the two genomes. For each gene, Liftoff finds the alignments of the exons that maximize sequence identity while preserving the transcript and gene structure.  If two genes incorrectly map to overlapping loci, Liftoff determines which gene is most-likely mis-mapped, and attempts to re-map it. Liftoff can also find additional gene copies present in the target assembly that are not annotated in the reference.



### Getting Started

#### Installation

Liftoff requires Python 3.12 or newer and [minimap2](https://github.com/lh3/minimap2)
(>= 2.28) on your `PATH` (or pass its location with `--minimap2`).

> **Note:** the releases currently published on Bioconda and PyPI are Liftoff 1.6.x, which
> use the older command line described in [Migrating from Liftoff 1.x](#migrating-from-liftoff-1x).
> Until a new release is published, install this version from source.

Install from source into a conda environment that also provides minimap2:

```
git clone https://github.com/agshumate/Liftoff liftoff
cd liftoff
mamba env create -f environment.yml
mamba activate liftoff
```

Or install with pip into an existing environment (minimap2 must be installed separately):

```
pip install "git+https://github.com/agshumate/Liftoff"
```

The conda environment also contains the development tools; run `pre-commit install` to enable
the git hooks and `pytest` to run the test suite.

#### Quick start

Lift the annotation of `reference.fa` onto `target.fa`:

```
liftoff -g reference.gff3 -o target.gff3 -u unmapped.txt target.fa reference.fa
```

This writes the lifted annotation to `target.gff3` and the IDs of features that could not be
lifted to `unmapped.txt`. The feature database built from the annotation is saved as
`reference.gff3_db` and can be reused in later runs, and intermediate FASTA/SAM files are kept
in `intermediate_files/`. Without `-o` the annotation is written to standard output.

Reuse the feature database, lift chromosome by chromosome, place genes from unplaced
reference sequences, search for extra gene copies and run four alignments in parallel:

```
liftoff --db reference.gff3_db -o target.gff3 -u unmapped.txt \
    --chroms chroms.txt --unplaced unplaced.txt \
    --copies --copy-identity 0.95 -p 4 \
    target.fa reference.fa
```

Repair CDSs broken by the lift-over (writes `target.gff3` and `target.gff3_polished`):

```
liftoff -g reference.gff3 -o target.gff3 --polish target.fa reference.fa
```

Lift every top-level feature type in the annotation (not only genes), except chromosome-length
`region` records and centromeres:

```
liftoff -g reference.gff3 -o target.gff3 --all-feature-types \
    --exclude-feature-types region,centromere target.fa reference.fa
```

Lift additional top-level feature types listed in a file and pass extra options to minimap2:

```
liftoff -g reference.gff3 -o target.gff3 -f feature_types.txt \
    --mm2-options="-r 2k -z 5000" target.fa reference.fa
```

### Usage
```
usage: liftoff [-h] [-V] (-g GFF | --db DB) [-o FILE] [-u FILE] [--exclude-partial]
               [--intermediate-dir DIR] [--mm2-options STR] [--minimap2 PATH] [--alignments DIR]
               [-a A] [-s S] [-d D] [--flank F] [-p P] [-f FILE | --all-feature-types]
               [--exclude-feature-types TYPES] [--infer-genes] [--infer-transcripts]
               [--chroms FILE] [--unplaced FILE] [--copies] [--copy-identity SC] [--max-overlap O]
               [--mismatch M] [--gap-open GO] [--gap-extend GE] [--polish] [--cds | --no-cds]
               [-v | -q]
               target reference

Lift features from one genome assembly to another.

options:
  -h, --help            show this help message and exit
  -V, --version         show program's version number and exit

Required input (sequences):
  target                target FASTA genome to lift features to
  reference             reference FASTA genome to lift features from

Required input (annotation):
  -g, --gff GFF         annotation file to lift over in GFF3 or GTF format
  --db DB               feature database built by a previous run (written next to the GFF as
                        <GFF>_db)

Output:
  -o, --output FILE     write lifted annotation to FILE in the input format; 'stdout' writes to
                        the terminal (default: stdout)
  -u, --unmapped FILE   write IDs of unmapped features to FILE (default: unmapped_features.txt)
  --exclude-partial     report mappings below --min-coverage/--min-identity as unmapped instead of
                        writing them with partial_mapping=True / low_identity=True
  --intermediate-dir DIR
                        directory for intermediate FASTA and SAM files (default:
                        intermediate_files)

Alignment:
  --mm2-options STR     space-delimited minimap2 options; use --mm2-options='-r 2k -z 5000' so
                        values starting with '-' are not parsed as Liftoff options (default: -a
                        --end-bonus 5 --eqx -N 50 -p 0.5)
  --minimap2 PATH       minimap2 executable (default: search PATH)
  --alignments DIR      use pre-computed SAM files from DIR instead of running minimap2
  -a, --min-coverage A  designate a feature mapped only if its child features align with coverage
                        >= A (default: 0.5)
  -s, --min-identity S  designate a feature mapped only if its child features align with identity
                        >= S (default: 0.5)
  -d, --distance-factor D
                        alignment blocks further apart in the target than D times their distance
                        in the reference are not chained (default: 2.0)
  --flank F             flanking sequence to align, as a fraction [0-1] of gene length (default:
                        0.0)

Lift-over settings:
  -p, --threads P       parallel alignment processes (default: 1)
  -f, --feature-types FILE
                        file listing additional top-level feature types to lift (one per line)
  --all-feature-types   lift every top-level feature type in the annotation, not only genes
  --exclude-feature-types TYPES
                        comma-separated top-level feature types not to lift, e.g.
                        region,centromere (may be repeated)
  --infer-genes         annotation only contains transcripts and exon/CDS features
  --infer-transcripts   annotation only contains genes and exon/CDS features
  --chroms FILE         comma separated file of corresponding reference,target chromosomes
  --unplaced FILE       file of unplaced reference sequence names to lift after --chroms
  --copies              look for extra gene copies in the target genome
  --copy-identity SC    with --copies, minimum exon/CDS identity for an extra copy (>= --min-
                        identity) (default: 1.0)
  --max-overlap O       maximum fraction [0-1] of overlap allowed between two features (default:
                        0.1)
  --mismatch M          mismatch penalty in exons (default: 2)
  --gap-open GO         gap open penalty in exons (default: 2)
  --gap-extend GE       gap extend penalty in exons (default: 1)
  --polish              re-align exons to repair broken CDS; writes <output> and <output>_polished
  --cds, --no-cds       annotate CDS status (partial, missing start/stop, in-frame stop codon)
                        (default: True)

Logging:
  -v, --verbose         show debug messages and minimap2 output
  -q, --quiet           only show warnings and errors
```

#### Migrating from Liftoff 1.x

This version uses conventional `--long-option` names. The old single-dash long options
have been renamed as follows (single-letter options are unchanged):

| Liftoff 1.x | Current |
|---|---|
| `-db DB` | `--db DB` |
| `-exclude_partial` | `--exclude-partial` |
| `-dir DIR` | `--intermediate-dir DIR` |
| `-mm2_options="..."` | `--mm2-options="..."` |
| `-flank F` | `--flank F` |
| `-m PATH` | `--minimap2 PATH` |
| `-infer_genes` / `-infer_transcripts` | `--infer-genes` / `--infer-transcripts` |
| `-chroms TXT` / `-unplaced TXT` | `--chroms FILE` / `--unplaced FILE` |
| `-copies` / `-sc SC` | `--copies` / `--copy-identity SC` |
| `-overlap O` | `--max-overlap O` |
| `-mismatch` / `-gap_open` / `-gap_extend` | `--mismatch` / `--gap-open` / `--gap-extend` |
| `-polish` | `--polish` |
| `-cds` (could not be disabled) | `--cds` / `--no-cds` |
| — | `--all-feature-types` (lift every top-level feature type) |
| — | `--exclude-feature-types TYPES` (top-level feature types not to lift) |
| — | `--alignments DIR` (use pre-computed SAM files) |
| — | `-v/--verbose`, `-q/--quiet` |

Every short option (`-g`, `-o`, `-u`, `-a`, `-s`, `-d`, `-p`, `-f`, `-V`) also has a long
form, e.g. `-g/--gff`, `-a/--min-coverage`, `-p/--threads`. Progress messages are now
written to standard error through Python's `logging` module.

### Input
The only required inputs are the reference genome sequence(fasta format), the target genome sequence(fasta format) and the reference annotation or feature database. If an annotation file is provided with the `-g/--gff` argument, a feature database will be built automatically and can be used for future lift overs by providing the `--db` argument. Annotation files may be gzip-compressed. FASTA files may be compressed with BGZF (`bgzip`), but not plain `gzip`. All input files must be in writable directories, because indices (`.fai`, the minimap2 `.mmi` index) and the feature database (`<annotation>_db`) are created next to them.

### Feature Types
By default, 'gene' features and all child features of genes (i.e. transcripts, mRNA, exons, CDS, UTRs) will be lifted over. The `-f/--feature-types` option can be used to specify a file containing a list of additional parent feature types you wish to lift-over. Example of a feature types file would be the following:

```
biological_region
miRNA
repeat_element
```

Alternatively, `--all-feature-types` lifts every top-level feature in the annotation — any feature without a parent, whatever its type — together with its child features (it cannot be combined with `-f`). This also includes features such as chromosome-length `region` records, and features too short to align with minimap2 (for example yeast centromere elements of about 25 bp) are reported as unmapped. Lifting additional feature types does not change how genes are lifted.

Use `--exclude-feature-types` with a comma-separated list (the option may be repeated) to skip top-level feature types, whichever way the types were selected. For example, `--all-feature-types --exclude-feature-types region,centromere` lifts everything except chromosome regions and centromeres, which avoids aligning whole chromosomes as single features. Excluded features and their children are neither lifted nor listed as unmapped. Only top-level types can be excluded; Liftoff warns about any listed type that does not occur as a top-level feature (for example a misspelling, or a child type such as `exon`).

### Feature IDs
Feature IDs may not contain spaces and must be unique: two different features may not share an `ID`. When building the feature database from a GFF3 file, Liftoff checks this and stops with an error listing the conflicting IDs and their line numbers. A single feature spanning several lines (for example a CDS split across exons) may repeat its `ID` on each line, provided those lines have the same sequence, feature type, strand and `Parent`. GTF files carry no `ID` attribute and are not checked, and neither are pre-built databases passed with `--db`.

### Sequence Identity and Alignment Coverage
A gene will be considered mapped successfully if the alignment coverage and sequence identity in the child features (usually exons/CDS) is >= 50%. This can be changed with the `-a/--min-coverage` and `-s/--min-identity` options. By default, genes that map below these thresholds will be included in the gff file with partial_mapping=True and low_identity=True in the last column. To exclude these partial/low identity mappings from the final GFF use `--exclude-partial`, and these genes will instead be written to the unmapped_features.txt file. The sequence identity and alignment coverage is reported in the final column of the output GFF for each gene.

### Minimap2 parameters
By default liftoff uses the following parameters for the minimap2 alignments `-a --eqx --end-bonus 5  -N 50 -p 0.5`
`-a` and `--eqx` specify that the output should be in SAM format with the cigar string including "=" for matches and "X" for mismatches (opposed to the default SAM format using 'M' for both). The `-N` and `-p` parameters specified allow for more secondary alignments to be considered which is helpful in the resolution of multi-gene families. The `--end-bonus` parameter favors end-to-end alignments of the gene over soft clipping a mismatched base at the start or end of the alignment. For example if the stop codon of the reference gene is TAA and the stop codon of the target gene is TAG, without the end-bonus parameter, this alignment and subsequent annotation would be truncated by 1 base.

The user may wish to change the minimap2 parameters for their specific data. This can be done with the `--mm2-options` option and a quoted string of options to add or change; the required options above are added automatically if they are missing. Join the value to the option name with "=" (as below) so that a value starting with "-" is never mistaken for a Liftoff option. For more divergent species in particular, increasing the `-r` and `-z` parameters may improve results (see Minimap2 documentation for more details). An example of changing these with `--mm2-options` would be

```
--mm2-options="-r 2k -z 5000"
```
### Polishing Exon/CDS Annotations
With the `--polish` option Liftoff will re-align the exons in attempt to restore proper coding sequences in cases where the lift-over resulted in start/stop codon loss or introduced an in-frame stop codon. This will increase the run time but offers improvements in preserving proper CDS annotations. With the polish option, 2 output files are created: `{output}` contains the annotations prior to the polishing step and `{output}_polished` contains the annotations after being polished (when writing to standard output, both are written there in turn). A polished gene replaces the original lift-over only if it has more valid ORFs, or the same number with higher sequence identity (or equal identity and higher coverage).

### CDS Trimming
When a lifted CDS no longer forms a valid ORF (for example because an indel in the target introduces a frameshift or premature stop codon) but still contains an intact ORF of at least 180 bases, the CDS features of that transcript are trimmed to the longest such ORF (searched in all three reading frames) and their phases are recomputed. CDS features falling entirely outside the ORF are removed. The `valid_ORF` and codon attributes then describe the trimmed CDS, while `matches_ref_protein` compares the untrimmed translation with the reference protein. Transcripts without such an ORF keep their full lifted CDS and are candidates for `--polish`.

### Gene Structure in Cross-Species Lift-over
Liftoff works best when the gene structure (i.e intron size) is similar in the reference and target genomes. When genes differ significantly in size, the alignments are more fragmented and often small exons at the beginning or end of the gene are not aligned. Adding and aligning some percentage of flanking sequence to the gene with the `--flank` option can improve this in some cases. Additionally increasing the `-d/--distance-factor` option will allow mappings where the genes are much larger in the target genome than in the reference.

### Chromosome by Chromosome Lift-over
By default, all genes will be aligned to the entire target assembly. However, for chromosome-scale assemblies of the same species, the `--chroms` option can be used to perform the lift-over chromosome by chromosome which improves accuracy. After the chromosome by chromosome lift over is complete, any genes that did not map will be aligned to the whole genome. This is strongly recommended for repetitive/polyploid genomes where there are many similar genes on different chromosomes. This option can be enabled by providing a comma separated file `chroms.txt` with corresponding chromosome names with the `--chroms` argument. Each line of the file should follow {ref_chrom_name},{target_chrom_name} for each pair of corresponding chromosomes. For example, a lift over from a Genbank human assembly to a Refseq human assembly would have the following `chroms.txt` file.
 ```
chr1,NC_000001.10
chr2,NC_000002.11
chr3,NC_000003.11
chr4,NC_000004.11
chr5,NC_000005.9
chr6,NC_000006.11
chr7,NC_000007.13
chr8,NC_000008.10
chr9,NC_000009.11
chr10,NC_000010.10
chr11,NC_000011.9
chr12,NC_000012.11
chr13,NC_000013.10
chr14,NC_000014.8
chr15,NC_000015.9
chr16,NC_000016.9
chr17,NC_000017.10
chr18,NC_000018.9
chr19,NC_000019.9
chr20,NC_000020.10
chr21,NC_000021.8
chr22,NC_000022.10
chrX,NC_000023.10
chrY,NC_000024.9
```

#### Unplaced Genes
A file of unplaced reference sequence names (one per line) can be provided with the `--unplaced` option, which requires `--chroms`. With this option, genes from these unplaced contigs in the reference will be mapped to the whole target assembly after the genes on the main chromosomes in `chroms.txt` have been mapped. For example:
```
chrUn_KI270302v1
chrUn_KI270304v1
```


### Extra Gene Copies
With the `--copies` option, Liftoff will look for extra copies of genes that are not annotated in the reference after the initial lift over. A gene copy will only be annotated at a locus if it does not overlap another annotated feature. By default, exons/CDSs must have 100% sequence identity; this threshold can be lowered with `--copy-identity` (which must be at least `--min-identity`). Extra gene copies keep the ID of the reference gene with a `_{copy_number}` suffix (also applied to their child features and `Parent`/`*_id` attributes) and are tagged with `extra_copy_number={copy_number}` in the last column of the GFF file.

### Output
The output is a file in the same format as the reference annotation (GFF3 or GTF) for the target genome and a file with the IDs of unmapped genes (`-u`, default `unmapped_features.txt`). The 9th column of the target annotation will contain the same information as the original reference plus the following

#### Genes:
```
sequence_ID: The sequence identity of the gene compared to the reference in exon regions
coverage: The alignment coverage of the gene in exon regions
valid_ORFs: The number of valid ORFs annotated within the gene
copy_num_ID: The gene ID with its copy number appended, e.g. gene1_0
partial_mapping: True if the coverage is below --min-coverage
low_identity: True if the sequence identity is below --min-identity
```

#### Transcripts:
```
valid_ORF: Indicates the CDS annotation properly starts with a start codon, ends with a stop codon,
           and does not have any in-frame stop codons (after any CDS trimming, see above).
matches_ref_protein: Indicates the translated CDS matches the reference CDS exactly
partial_ORF: Indicates the CDS encodes fewer than 3 amino acids
missing_start_codon: Indicates the CDS does not begin with a start codon
missing_stop_codon: Indicates the CDS does not end with a stop codon
inframe_stop_codon: Indicates the CDS has an inframe stop codon.
```

#### All features:
```
extra_copy_number: The copy number increase of this feature compared to the reference.
extra_copy_number=0 means this is the original reference gene.
```

### Alignment backends and WebAssembly (Pyodide)

All of Liftoff's Python dependencies are pure Python or available in
[Pyodide](https://pyodide.org) (version 314.0 or newer, which provides numpy >= 2.4,
biopython >= 1.87 and networkx >= 3.6), so Liftoff can run in the browser or under Node.js.
minimap2 cannot be launched as a subprocess there, so the alignment step is pluggable:

* By default, `Minimap2Aligner` runs the `minimap2` executable.
* `--alignments DIR` (`PrecomputedSamAligner`) reads SAM files produced elsewhere. Each
  pipeline stage expects a file with a specific name, such as `reference_all_to_target_all.sam`.
  If one is missing, Liftoff stops with an error showing the minimap2 command that produces it
  from the gene sequences it has just written to the intermediate directory; run that command
  and start Liftoff again (runs with several stages, e.g. `--copies`, may need this repeated).
* From Python, `CallableAligner` delegates each alignment to your own function — for example
  one that calls a WebAssembly build of minimap2. This example runs the minimap2 executable:

```python
import subprocess

from liftoff.align.base import AlignmentJob
from liftoff.align.external import CallableAligner
from liftoff.config import LiftoffConfig
from liftoff.pipeline import run_liftoff


def align(job: AlignmentJob) -> None:
    """Align job.features_fasta to job.target_fasta, writing SAM to job.output_sam."""
    with job.output_sam.open('w') as sam:
        subprocess.run(
            ['minimap2', *job.minimap2_options, str(job.target_fasta), str(job.features_fasta)],
            stdout=sam,
            check=True,
        )


config = LiftoffConfig(
    target='target.fa',
    reference='reference.fa',
    gff='reference.gff3',
    output='target.gff3',
    unmapped='unmapped.txt',
)
result = run_liftoff(config, aligner=CallableAligner(align))
print(len(result.lifted_features), 'lifted;', len(result.unmapped_features), 'unmapped')
```

Polishing (`--polish`) uses a built-in NumPy implementation of semi-global alignment that
reproduces the parasail library's results exactly, so no compiled extensions are needed.

#### Feature database (SQLite) under WebAssembly

The reference annotation is loaded into an SQLite database by
[gffutils](https://github.com/daler/gffutils) through Python's standard `sqlite3` module,
which Pyodide includes. The test suite, including database creation and complete lift-overs
from pre-computed alignments, runs under Pyodide on Node.js. When planning a WebAssembly
deployment, keep in mind:

* **Writable storage.** The database is written next to the annotation as `<annotation>_db`,
  and FASTA indices are written next to the FASTA files. In the browser, place the input files
  in Pyodide's in-memory file system (for example with `pyodide.FS.writeFile`) rather than on a
  read-only mount.
* **Memory.** The in-memory file system lives inside the WebAssembly heap, which is limited to
  4 GB and often less in browsers. The database is roughly 2.3 times the size of the uncompressed
  GFF3 (the 6.7 MB yeast annotation produces a 15.7 MB database), and genome FASTA files,
  indices and intermediate files add to that. Small and medium genomes are practical in the
  browser; large (e.g. mammalian) annotations are better processed under Node.js, where files
  can be stored on disk.
* **Persistence.** Files in the in-memory file system disappear when the page is reloaded. To
  reuse a database with `--db`, save it to persistent storage (such as Emscripten's IDBFS or
  a directory mounted with `pyodide.mountNativeFS`).
* **Speed.** Building the yeast database takes about 2.5 s under Pyodide compared with 0.7 s
  natively. SQLite is used from a single thread, so the lack of threads in WebAssembly does not
  affect it.

## Future Development

**Performance of the pure-Python parasail implementation.** Liftoff no longer depends on
[parasail](https://github.com/jeffdaily/parasail): `liftoff.align.semiglobal` re-implements
`sg_dx_trace_scan_sat` in NumPy and reproduces parasail's scores and tracebacks exactly
(verified against recorded parasail output). It is, however, roughly 20 times slower than
parasail's SIMD C code for the same alignments (25 exon alignments from the chromosome I test
data take about 1.05 s versus 0.055 s), and its traceback table uses one byte per cell, so long
exons need a lot of memory (a 15 kb exon aligned to a 15 kb region needs about 225 MB). The
overall effect on typical runs is small, but it grows with the number of genes that need
polishing. Future work on this implementation could include:

* reducing per-column overhead in the NumPy dynamic-programming loop (e.g. anti-diagonal or
  striped vectorisation and fewer temporary arrays) and the size of the traceback table;
* optionally using parasail (or another compiled aligner) when it is installed, while keeping the
  pure-Python implementation as the fallback for Pyodide;
* publishing the implementation as a standalone, parasail-compatible package so it can be tested
  against more of parasail's alignment modes and reused by other pure-Python tools.

## Known Issues
Extracting gene sequences from bgzipped files is possible but much slower. It is recommended to decompress bgzipped FASTA files first.

## Citation
If you use Liftoff in your work please cite<br/>

Shumate, Alaina, and Steven L. Salzberg. 2020. “Liftoff: Accurate Mapping of Gene Annotations.” Bioinformatics , December. https://doi.org/10.1093/bioinformatics/btaa1016.
