"""Build the chromosome I subset in ``tests/data/chr1``.

Extracts chromosome I (BK006935.2) and its annotation from the yeast test data
and writes a copy of the sequence with deterministic CDS mutations (indels,
premature stops, substitutions and one ``N``) so that polishing is exercised::

    python make_chr1_data.py tests/data/yeast OUT_DIR

Output files were renamed to ``ref.fa``, ``target_mutated.fa``,
``annotation.gff3`` and ``mutations.tsv`` when committed.
"""

import gzip
from pathlib import Path
import random
import sys

CHROM = 'BK006935.2'
src, out = Path(sys.argv[1]), Path(sys.argv[2])
out.mkdir(parents=True, exist_ok=True)

# --- reference sequence -------------------------------------------------------
seq_lines: list[str] = []
keep = False
with gzip.open(src / 'GCA_000146045.2_R64_genomic.fna.gz', 'rt') as fh:
    for line in fh:
        if line.startswith('>'):
            keep = line[1:].split()[0] == CHROM
            continue
        if keep:
            seq_lines.append(line.strip())
ref_seq = ''.join(seq_lines)


def write_fasta(path: Path, name: str, seq: str) -> None:
    with open(path, 'w') as fh:
        fh.write(f'>{name}\n')
        for i in range(0, len(seq), 80):
            fh.write(seq[i : i + 80] + '\n')


write_fasta(out / 'chr1_ref.fa', CHROM, ref_seq)

# --- annotation subset ----------------------------------------------------------
gff_lines = []
cds: list[tuple[int, int, str]] = []
with gzip.open(src / 'GCA_000146045.2_R64_genomic.gff.gz', 'rt') as fh:
    for line in fh:
        if line.startswith('##gff-version'):
            gff_lines.append(line)
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if f[0] != CHROM:
            continue
        gff_lines.append(line)
        if f[2] == 'CDS':
            cds.append((int(f[3]), int(f[4]), f[8].split('Parent=')[1].split(';')[0]))
with open(out / 'chr1.gff3', 'w') as fh:
    fh.writelines(gff_lines)

# --- mutate CDSs ------------------------------------------------------------------
rng = random.Random(20260914)
by_parent: dict[str, list[tuple[int, int]]] = {}
for s, e, p in cds:
    by_parent.setdefault(p, []).append((s, e))
parents = sorted(by_parent)
chosen = [p for i, p in enumerate(parents) if i % 3 == 0]
mutations: list[tuple[int, str, str]] = []  # (1-based pos, kind, payload)
kinds = ['del1', 'ins2', 'stop', 'del4', 'ins1', 'sub']
for n, p in enumerate(chosen):
    s, e = max(by_parent[p], key=lambda x: x[1] - x[0])
    if e - s < 60:
        continue
    pos = rng.randint(s + 20, e - 20)
    mutations.append((pos, kinds[n % len(kinds)], ''.join(rng.choice('ACGT') for _ in range(2))))

seq = list(ref_seq)
for pos, kind, payload in sorted(mutations, reverse=True):
    i = pos - 1
    if kind == 'del1':
        del seq[i]
    elif kind == 'del4':
        del seq[i : i + 4]
    elif kind == 'ins2':
        seq[i:i] = list(payload)
    elif kind == 'ins1':
        seq[i:i] = [payload[0]]
    elif kind == 'stop':
        seq[i : i + 3] = list('TAA')
    elif kind == 'sub':
        seq[i] = {'A': 'C', 'C': 'G', 'G': 'T', 'T': 'A'}.get(seq[i], 'A')
# Add an ambiguous base inside one CDS to exercise the unknown-character mapping.
if mutations:
    seq[mutations[0][0] + 5] = 'N'
write_fasta(out / 'chr1_target_mutated.fa', CHROM, ''.join(seq))
with open(out / 'chr1_mutations.tsv', 'w') as fh:
    fh.write('position\tkind\tpayload\n')
    for m in sorted(mutations):
        fh.write('\t'.join(map(str, m)) + '\n')
print(len(mutations), 'mutations over', len(chosen), 'transcripts')
