"""Build the synthetic parasail corpus ``tests/data/parasail/sg_dx_synthetic.json.gz``.

Run with the ORIGINAL Liftoff 1.6.3 importable and ``parasail-python`` installed::

    python make_synthetic.py sg_dx_synthetic.json.gz
"""

import gzip
import json
import random
import sys

import parasail

from liftoff.polish import make_scoring_matrix

rng = random.Random(7)
matrix = make_scoring_matrix(3)
cases = []


def rec(s1, s2, o=10, e=1, tag=''):
    r = parasail.sg_dx_trace_scan_sat(s1, s2, o, e, matrix)
    cases.append(
        {
            'tag': tag,
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


def mutate(s, rate):
    out = []
    for c in s:
        x = rng.random()
        if x < rate / 3:
            continue
        if x < 2 * rate / 3:
            out.append(c)
            out.append(rng.choice('acgt'))
            continue
        if x < rate:
            out.append(rng.choice('acgt'))
            continue
        out.append(c)
    return ''.join(out)


def capitalise(s):
    # random upper-case stretches, mimicking CDS/splice-site capitalisation
    s = list(s)
    for _ in range(rng.randint(0, 3)):
        a = rng.randrange(len(s))
        b = min(len(s), a + rng.randint(1, 60))
        for i in range(a, b):
            s[i] = s[i].upper()
    return ''.join(s)


# random realistic cases
for _ in range(300):
    L = rng.choice([1, 2, 3, 5, 10, 30, 80, 150, 400])
    core = ''.join(rng.choice('acgt') for _ in range(L))
    s1 = capitalise(core)
    flank = rng.randint(0, 40)
    left = ''.join(rng.choice('acgt') for _ in range(flank))
    right = ''.join(rng.choice('acgt') for _ in range(rng.randint(0, 40)))
    s2 = left + mutate(core, rng.choice([0, 0.02, 0.1, 0.3])) + right
    if not s2:
        s2 = 'a'
    rec(s1, s2, tag='random')
# tie-heavy repeats
for unit in ['a', 'at', 'acg', 'aacc']:
    for k in [1, 3, 10, 25]:
        for j in [0, 1, 2, 5]:
            rec((unit * k).upper(), unit * (k + j), tag='repeat')
            rec(unit * (k + j), unit * k, tag='repeat_rev')
# unknown characters: N, n, IUPAC, '*'
for _ in range(40):
    core = ''.join(rng.choice('acgtn') for _ in range(rng.randint(5, 60)))
    s2 = ''.join(rng.choice('acgtnNRY*') if rng.random() < 0.1 else c for c in core)
    rec(capitalise(core), s2, tag='ambiguous')
# completely dissimilar and very different lengths
for _ in range(20):
    rec('A' * rng.randint(1, 50), 'c' * rng.randint(1, 50), tag='dissimilar')
    rec(
        ''.join(rng.choice('ACGT') for _ in range(rng.randint(100, 200))),
        ''.join(rng.choice('acgt') for _ in range(rng.randint(1, 20))),
        tag='long_short',
    )
# other gap penalties (general port correctness)
for o, e in [(1, 1), (3, 1), (5, 2), (0, 0), (12, 3)]:
    for _ in range(10):
        core = ''.join(rng.choice('acgt') for _ in range(rng.randint(10, 120)))
        rec(capitalise(core), mutate(core, 0.15), o, e, tag=f'gaps_{o}_{e}')
# long case that forces 16-bit (8-bit saturates)
core = ''.join(rng.choice('acgt') for _ in range(3000))
rec(
    capitalise(core),
    ''.join(rng.choice('acgt') for _ in range(200)) + mutate(core, 0.05),
    tag='long16',
)

with gzip.open(sys.argv[1], 'wt') as fh:
    json.dump(cases, fh)
print(len(cases))
