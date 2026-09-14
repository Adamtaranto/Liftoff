"""Lift-over of chromosome I onto a CDS-mutated copy.

The minimap2 tests exercise the full native pipeline. The pre-computed tests
feed the SAM file recorded from minimap2 through the alternative backends, so
they also run where minimap2 is unavailable (e.g. under Pyodide).
"""

from __future__ import annotations

from collections.abc import Callable
import gzip
import json
from pathlib import Path
import shutil

import pytest

from liftoff import cds
from liftoff import polish as polish_module
from liftoff.align.base import AlignmentJob
from liftoff.align.external import CallableAligner
from liftoff.cli import main
from liftoff.config import LiftoffConfig
from liftoff.errors import AlignmentError
from liftoff.pipeline import run_liftoff

SAM_NAME = 'reference_all_to_target_all.sam'


def _cli_args(inputs: dict[str, Path], workdir: Path, *extra: str) -> list[str]:
    return [
        '-q',
        '--gff',
        str(inputs['gff']),
        '--output',
        str(workdir / 'out.gff3'),
        '--unmapped',
        str(workdir / 'unmapped.txt'),
        '--intermediate-dir',
        str(workdir / 'intermediate'),
        *extra,
        str(inputs['target']),
        str(inputs['reference']),
    ]


def _config(inputs: dict[str, Path], workdir: Path, **overrides: object) -> LiftoffConfig:
    return LiftoffConfig(
        target=str(inputs['target']),
        reference=str(inputs['reference']),
        gff=str(inputs['gff']),
        output=str(workdir / 'out.gff3'),
        unmapped=str(workdir / 'unmapped.txt'),
        intermediate_dir=str(workdir / 'intermediate'),
        **overrides,  # type: ignore[arg-type]
    )


@pytest.fixture
def precomputed_dir(tmp_path: Path, data_dir: Path) -> Path:
    directory = tmp_path / 'alignments'
    directory.mkdir()
    shutil.copy(data_dir / 'expected' / f'chr1_basic_{SAM_NAME}', directory / SAM_NAME)
    return directory


@pytest.mark.minimap2
@pytest.mark.parametrize(
    ('run', 'extra'),
    [('basic', []), ('copies', ['--copies', '--copy-identity', '0.95'])],
)
def test_chr1_with_minimap2(
    run: str,
    extra: list[str],
    tmp_path: Path,
    data_dir: Path,
    chr1_inputs: dict[str, Path],
    same_annotation: Callable[[Path, Path], None],
) -> None:
    assert main(_cli_args(chr1_inputs, tmp_path, *extra)) == 0
    same_annotation(tmp_path / 'out.gff3', data_dir / 'expected' / f'chr1_{run}.gff3')
    expected_unmapped = (data_dir / 'expected' / f'chr1_{run}_unmapped.txt').read_text()
    assert (tmp_path / 'unmapped.txt').read_text() == expected_unmapped


@pytest.mark.minimap2
def test_chr1_polish_with_minimap2(
    tmp_path: Path,
    data_dir: Path,
    chr1_inputs: dict[str, Path],
    same_annotation: Callable[[Path, Path], None],
) -> None:
    assert main(_cli_args(chr1_inputs, tmp_path, '--polish', '--threads', '2')) == 0
    same_annotation(tmp_path / 'out.gff3', data_dir / 'expected' / 'chr1_polish.gff3')
    same_annotation(
        tmp_path / 'out.gff3_polished', data_dir / 'expected' / 'chr1_polish_polished.gff3'
    )


def test_chr1_precomputed_alignments(
    tmp_path: Path,
    data_dir: Path,
    chr1_inputs: dict[str, Path],
    precomputed_dir: Path,
    same_annotation: Callable[[Path, Path], None],
) -> None:
    config = _config(chr1_inputs, tmp_path, alignments=str(precomputed_dir))
    result = run_liftoff(config)
    assert result.output_paths == [config.output]
    same_annotation(tmp_path / 'out.gff3', data_dir / 'expected' / 'chr1_basic.gff3')
    assert result.unmapped_features == []


def test_chr1_callable_aligner(
    tmp_path: Path,
    data_dir: Path,
    chr1_inputs: dict[str, Path],
    precomputed_dir: Path,
    same_annotation: Callable[[Path, Path], None],
) -> None:
    jobs: list[AlignmentJob] = []

    def fake_minimap2(job: AlignmentJob) -> None:
        jobs.append(job)
        assert job.features_fasta.is_file()
        assert '--eqx' in job.minimap2_options
        shutil.copy(precomputed_dir / job.output_sam.name, job.output_sam)

    run_liftoff(_config(chr1_inputs, tmp_path, threads=4), CallableAligner(fake_minimap2))
    assert len(jobs) == 1
    same_annotation(tmp_path / 'out.gff3', data_dir / 'expected' / 'chr1_basic.gff3')


def test_missing_precomputed_alignment_is_reported(
    tmp_path: Path, chr1_inputs: dict[str, Path]
) -> None:
    empty = tmp_path / 'empty'
    empty.mkdir()
    with pytest.raises(AlignmentError, match=SAM_NAME):
        run_liftoff(_config(chr1_inputs, tmp_path, alignments=str(empty)))


def test_chr1_polish_matches_parasail_trace(
    tmp_path: Path,
    data_dir: Path,
    chr1_inputs: dict[str, Path],
    precomputed_dir: Path,
    monkeypatch: pytest.MonkeyPatch,
    same_annotation: Callable[[Path, Path], None],
) -> None:
    """Polishing reproduces the SAM records and candidate genes of parasail-based Liftoff.

    The trace was recorded with Liftoff 1.6.3, which never trimmed broken CDSs
    to their longest ORF. Trimming is disabled here so the same 23 genes are
    polished, verifying the NumPy aligner (and the non-mutating splice-site
    handling) against parasail on every exon alignment.
    """
    with gzip.open(data_dir / 'parasail' / 'polish_trace_chr1.json.gz', 'rt') as handle:
        trace = json.load(handle)

    polish_sams: list[dict[str, str]] = []
    original_polish = polish_module.polish_annotations

    def recording_polish(*args: object, **kwargs: object) -> bool:
        polished = original_polish(*args, **kwargs)  # type: ignore[arg-type]
        if polished:
            intermediate = args[3]
            assert isinstance(intermediate, Path)
            sam_text = (intermediate / polish_module.POLISH_SAM_NAME).read_text()
            polish_sams.append({'feature': str(args[5]), 'sam': sam_text})
        return polished

    check_calls: list[list[list[object]]] = []
    original_check = cds.check_cds

    def recording_check(feature_list, *args, **kwargs):  # type: ignore[no-untyped-def]
        original_check(feature_list, *args, **kwargs)
        check_calls.append(
            [
                [key, f.id, f.featuretype, f.seqid, f.strand, f.start, f.end, f.frame, f.attributes]
                for key in sorted(feature_list)
                for f in feature_list[key]
            ]
        )

    monkeypatch.setattr(polish_module, 'polish_annotations', recording_polish)
    monkeypatch.setattr(cds, 'check_cds', recording_check)
    monkeypatch.setattr(cds, 'find_longest_orf', lambda *_args, **_kwargs: None)

    config = _config(chr1_inputs, tmp_path, alignments=str(precomputed_dir), polish=True)
    result = run_liftoff(config)

    assert result.output_paths == [config.output, config.output + '_polished']
    assert polish_sams == trace['polish_sams']
    # The last check_cds call annotates the polished candidate genes.
    assert json.loads(json.dumps(check_calls[-1])) == trace['check_cds_calls'][-1]
    same_annotation(
        tmp_path / 'out.gff3_polished',
        data_dir / 'expected' / 'chr1_polish_polished_without_orf_trimming.gff3',
    )


def test_chr1_polish_precomputed(
    tmp_path: Path,
    data_dir: Path,
    chr1_inputs: dict[str, Path],
    precomputed_dir: Path,
    same_annotation: Callable[[Path, Path], None],
) -> None:
    config = _config(chr1_inputs, tmp_path, alignments=str(precomputed_dir), polish=True)
    run_liftoff(config)
    same_annotation(tmp_path / 'out.gff3', data_dir / 'expected' / 'chr1_polish.gff3')
    same_annotation(
        tmp_path / 'out.gff3_polished', data_dir / 'expected' / 'chr1_polish_polished.gff3'
    )


@pytest.mark.minimap2
def test_chr1_flank_preserves_coordinates(tmp_path: Path, chr1_inputs: dict[str, Path]) -> None:
    """With --flank, lifted features (including single-level ones) keep annotated coordinates."""
    feature_types = tmp_path / 'types.txt'
    feature_types.write_text('origin_of_replication\ntelomere\nlong_terminal_repeat\n')
    reference = str(chr1_inputs['reference'])
    args = _cli_args(chr1_inputs, tmp_path, '--flank', '0.2', '--copies', '-f', str(feature_types))
    args[-2:] = [reference, reference]  # lift chromosome I onto itself
    assert main(args) == 0

    top_level = ('gene', 'origin_of_replication', 'telomere', 'long_terminal_repeat')

    def coordinates(path: Path) -> dict[str, tuple[int, int]]:
        found = {}
        for line in path.read_text().splitlines():
            fields = line.split('\t')
            if line.startswith('#') or len(fields) < 9 or fields[2] not in top_level:
                continue
            attributes = dict(item.split('=', 1) for item in fields[8].split(';'))
            if attributes.get('extra_copy_number', '0') == '0':
                found[attributes['ID']] = (int(fields[3]), int(fields[4]))
        return found

    lifted = coordinates(tmp_path / 'out.gff3')
    annotated = coordinates(chr1_inputs['gff'])
    assert len(lifted) > 100
    assert {k: v for k, v in lifted.items() if annotated[k] != v} == {}
