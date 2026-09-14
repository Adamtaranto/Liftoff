"""Lift chromosome I onto a CDS-mutated copy with every alignment backend.

Each ``(backend, scenario)`` combination is run once per module; the tests
below inspect the outputs of that run.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path
import shutil

import pytest

from liftoff.align.base import AlignmentJob
from liftoff.align.external import CallableAligner
from liftoff.config import LiftoffConfig
from liftoff.errors import AlignmentError
from liftoff.pipeline import LiftoffResult, run_liftoff

from .helpers import CHR1_SAM_NAMES, Chr1Inputs, RunPaths

SCENARIOS = {
    'basic': [],
    'copies': ['--copies', '--copy-identity', '0.95'],
}


@dataclass(frozen=True)
class Chr1Run:
    scenario: str
    paths: RunPaths


@pytest.fixture(
    scope='module',
    params=[
        pytest.param('minimap2', marks=pytest.mark.minimap2),
        pytest.param('precomputed'),
    ],
)
def cli_backend(request: pytest.FixtureRequest) -> str:
    """Alignment backend selected on the command line."""
    return str(request.param)


@pytest.fixture(scope='module', params=sorted(SCENARIOS))
def chr1_cli_run(
    request: pytest.FixtureRequest,
    cli_backend: str,
    tmp_path_factory: pytest.TempPathFactory,
    copy_chr1_inputs: Callable[[Path], Chr1Inputs],
    chr1_alignments: Path,
    run_cli: Callable[..., RunPaths],
) -> Chr1Run:
    """Run the CLI once for a backend and scenario."""
    scenario = str(request.param)
    workdir = tmp_path_factory.mktemp(f'chr1_{cli_backend}_{scenario}')
    inputs = copy_chr1_inputs(workdir / 'inputs')
    extra = list(SCENARIOS[scenario])
    if cli_backend == 'precomputed':
        extra += ['--alignments', str(chr1_alignments)]
    paths = run_cli(workdir, inputs.target, inputs.reference, inputs.gff, extra)
    return Chr1Run(scenario, paths)


def test_annotation_matches_expected(
    chr1_cli_run: Chr1Run,
    expected_dir: Path,
    assert_same_annotation: Callable[[Path, Path], None],
) -> None:
    expected = expected_dir / f'chr1_{chr1_cli_run.scenario}.gff3'
    assert_same_annotation(chr1_cli_run.paths.output, expected)


def test_unmapped_features_match_expected(chr1_cli_run: Chr1Run, expected_dir: Path) -> None:
    expected = expected_dir / f'chr1_{chr1_cli_run.scenario}_unmapped.txt'
    assert chr1_cli_run.paths.unmapped.read_text() == expected.read_text()


# ---------------------------------------------------------------------------
# Library API with a callable aligner
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class ApiRun:
    result: LiftoffResult
    config: LiftoffConfig
    jobs: list[AlignmentJob]


@pytest.fixture(scope='module')
def callable_run(
    tmp_path_factory: pytest.TempPathFactory,
    copy_chr1_inputs: Callable[[Path], Chr1Inputs],
    chr1_alignments: Path,
) -> ApiRun:
    """Run the pipeline through the API, serving alignments from a callable."""
    workdir = tmp_path_factory.mktemp('chr1_callable')
    inputs = copy_chr1_inputs(workdir / 'inputs')
    paths = RunPaths(workdir)
    jobs: list[AlignmentJob] = []

    def serve_alignment(job: AlignmentJob) -> None:
        jobs.append(job)
        shutil.copy(chr1_alignments / job.output_sam.name, job.output_sam)

    config = LiftoffConfig(
        target=str(inputs.target),
        reference=str(inputs.reference),
        gff=str(inputs.gff),
        output=str(paths.output),
        unmapped=str(paths.unmapped),
        intermediate_dir=str(paths.intermediate),
        threads=4,
    )
    result = run_liftoff(config, CallableAligner(serve_alignment))
    return ApiRun(result, config, jobs)


def test_callable_aligner_output_matches_expected(
    callable_run: ApiRun,
    expected_dir: Path,
    assert_same_annotation: Callable[[Path, Path], None],
) -> None:
    assert_same_annotation(Path(callable_run.config.output), expected_dir / 'chr1_basic.gff3')


def test_callable_aligner_receives_one_whole_genome_job(callable_run: ApiRun) -> None:
    (job,) = callable_run.jobs
    assert job.output_sam.name == CHR1_SAM_NAMES[0]
    assert job.target_fasta == Path(callable_run.config.target)
    assert job.threads == 4


def test_callable_aligner_job_queries_exist(callable_run: ApiRun) -> None:
    (job,) = callable_run.jobs
    assert job.features_fasta.read_text().startswith('>')


def test_callable_aligner_job_includes_required_minimap2_options(callable_run: ApiRun) -> None:
    (job,) = callable_run.jobs
    assert {'-a', '--eqx'} <= set(job.minimap2_options)


def test_result_reports_outputs_and_no_unmapped_features(callable_run: ApiRun) -> None:
    assert callable_run.result.output_paths == [callable_run.config.output]
    assert callable_run.result.unmapped_features == []


def test_missing_precomputed_alignment_names_expected_file(
    tmp_path: Path, copy_chr1_inputs: Callable[[Path], Chr1Inputs]
) -> None:
    inputs = copy_chr1_inputs(tmp_path / 'inputs')
    empty = tmp_path / 'empty'
    empty.mkdir()
    config = LiftoffConfig(
        target=str(inputs.target),
        reference=str(inputs.reference),
        gff=str(inputs.gff),
        output=str(tmp_path / 'out.gff3'),
        unmapped=str(tmp_path / 'unmapped.txt'),
        intermediate_dir=str(tmp_path / 'intermediate'),
        alignments=str(empty),
    )

    with pytest.raises(AlignmentError, match=CHR1_SAM_NAMES[0]):
        run_liftoff(config)
