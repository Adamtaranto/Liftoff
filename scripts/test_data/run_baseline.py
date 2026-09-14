"""Record reference outputs with the ORIGINAL Liftoff 1.6.3 (parasail + pysam).

This script documents how ``tests/data/expected`` and
``tests/data/parasail/*trace*`` were generated. It must be run from a checkout
of Liftoff 1.6.3 (commit 94d6073, where tests live in ``liftoff/tests``) in an
environment with ``parasail-python``, ``pysam`` and ``minimap2`` installed::

    python run_baseline.py OUT_DIR basic advanced chr1_basic chr1_copies chr1_polish

``chr1_*`` runs expect the output of ``make_chr1_data.py`` in
``OUT_DIR/../chr1data``. Header lines (``#``) of the GFF outputs and the
``@PG`` line of SAM files were stripped before committing.
"""

import gzip
import json
from pathlib import Path
import shutil
import sys

import parasail

import liftoff.polish
import liftoff.run_liftoff

T = Path('liftoff/tests')
OUT = Path(sys.argv[1])
OUT.mkdir(parents=True, exist_ok=True)
fa = str(T / 'GCA_000146045.2_R64_genomic.fna.gz')
gff = str(T / 'GCA_000146045.2_R64_genomic.gff.gz')

cases = []
_orig = parasail.sg_dx_trace_scan_sat


def wrapped(s1, s2, o, e, matrix):
    r = _orig(s1, s2, o, e, matrix)
    cases.append(
        {
            's1': s1,
            's2': s2,
            'open': o,
            'extend': e,
            'score': r.score,
            'query': r.traceback.query,
            'ref': r.traceback.ref,
            'end_query': r.end_query,
            'end_ref': r.end_ref,
        }
    )
    return r


polish_sams = []
check_cds_calls = []
_orig_pa = liftoff.polish.polish_annotations
_orig_cc = liftoff.run_liftoff.check_cds


def wrapped_pa(feature_list, ref_faidx, target_faidx, args, fh, target_feature):
    r = _orig_pa(feature_list, ref_faidx, target_faidx, args, fh, target_feature)
    if r:
        sam_text = Path(args.dir, 'polish.sam').read_text()
        polish_sams.append({'feature': target_feature, 'sam': sam_text})
    return r


def wrapped_cc(feature_list, feature_hierarchy, args):
    _orig_cc(feature_list, feature_hierarchy, args)
    rows = [
        [key, f.id, f.featuretype, f.seqid, f.strand, f.start, f.end, f.frame, dict(f.attributes)]
        for key in sorted(feature_list)
        for f in feature_list[key]
    ]
    check_cds_calls.append(rows)


def main():
    liftoff.polish.parasail.sg_dx_trace_scan_sat = wrapped
    liftoff.run_liftoff.polish.polish_annotations = wrapped_pa
    liftoff.run_liftoff.check_cds = wrapped_cc
    runs = {
        'basic': [],
        'advanced': [
            '-chroms',
            str(T / 'chroms.txt'),
            '-unplaced',
            str(T / 'unplaced.txt'),
            '-copies',
        ],
        'polish': ['-polish'],
    }
    C = Path(OUT.parent / 'chr1data')
    chr1 = {
        'chr1_basic': [],
        'chr1_polish': ['-polish'],
        'chr1_copies': ['-copies', '-sc', '0.95'],
    }
    for name in sys.argv[2:]:
        d = OUT / name
        d.mkdir(exist_ok=True)
        if name in chr1:
            args = [
                '-g',
                str(d / 'chr1.gff3'),
                '-o',
                str(d / 'out.gff3'),
                '-u',
                str(d / 'unmapped.txt'),
                '-dir',
                str(d / 'inter'),
                *chr1[name],
                str(d / 'target.fa'),
                str(d / 'ref.fa'),
            ]
            shutil.copy(C / 'chr1.gff3', d / 'chr1.gff3')
            shutil.copy(C / 'chr1_ref.fa', d / 'ref.fa')
            shutil.copy(C / 'chr1_target_mutated.fa', d / 'target.fa')
        else:
            args = [
                '-g',
                gff,
                '-o',
                str(d / 'out.gff3'),
                '-u',
                str(d / 'unmapped.txt'),
                '-dir',
                str(d / 'inter'),
                *runs[name],
                fa,
                fa,
            ]
        liftoff.run_liftoff.main(args)

    if cases:
        with gzip.open(OUT / 'sg_dx_cases.json.gz', 'wt') as fh:
            json.dump(cases, fh)
        print('recorded', len(cases), 'parasail calls')
    if polish_sams:
        with gzip.open(OUT / 'polish_trace.json.gz', 'wt') as fh:
            json.dump({'polish_sams': polish_sams, 'check_cds_calls': check_cds_calls[-1:]}, fh)
        print('recorded', len(polish_sams), 'polish sams;', len(check_cds_calls), 'check_cds calls')


if __name__ == '__main__':
    main()
