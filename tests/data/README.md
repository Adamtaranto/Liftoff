# Test data

| Directory | Contents |
|---|---|
| `yeast/` | *S. cerevisiae* R64 genome (BGZF) and annotation, chromosome correspondence (`chroms.txt`) and unplaced sequence list used by the whole-genome tests. |
| `chr1/` | Chromosome I (BK006935.2) subset: `ref.fa`, `annotation.gff3`, and `target_mutated.fa`, a copy with the CDS mutations listed in `mutations.tsv` (so that polishing has broken ORFs to repair). |
| `expected/` | Expected outputs. GFF header lines (`#...`) and the SAM `@PG` line were removed. `yeast_*` and `chr1_polish_polished_without_orf_trimming.gff3` were recorded with the original Liftoff 1.6.3 code (commit 94d6073, using parasail, pysam and minimap2 2.28). The other `chr1_*.gff3` files were regenerated after fixing CDS trimming to the longest ORF, which 1.6.3 never applied; each trimmed CDS was checked to translate to the longest intact ORF of the original lifted CDS, with correct phases. The yeast outputs are unaffected by that fix. |
| `parasail/` | `sg_dx_*.json.gz`: inputs and exact `parasail.sg_dx_trace_scan_sat` results (score, end position, traceback strings) recorded during polishing and for synthetic edge cases. `polish_trace_chr1.json.gz`: the per-gene `polish.sam` files and polished candidate features from the 1.6.3 chr1 polishing run (compared with ORF trimming disabled, which reproduces 1.6.3's choice of genes to polish). |

The scripts in `scripts/test_data/` document how these files were generated.

The pre-refactor expected files (`*_expected_basic.gff`, `*_expected_advanced.gff`) had become
stale: they predated upstream changes that write the CDS phase (column 8) and preserve
percent-encoded attribute values. The original test compared only the unmapped feature lists,
because its GFF comparison loop never executed. The files in `expected/` were regenerated
from the unmodified 1.6.3 code; apart from those two known differences they agree with the
old files feature for feature.
