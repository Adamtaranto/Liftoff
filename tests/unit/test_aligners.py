"""Tests for alignment backends and alignment job planning/dispatch."""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path
import stat
import sys

import pytest

from liftoff.align import runner
from liftoff.align.base import AlignmentJob
from liftoff.align.external import CallableAligner, PrecomputedSamAligner
from liftoff.align.minimap2 import Minimap2Aligner
from liftoff.errors import AlignmentError
from liftoff.models import LiftoverType

JobFactory = Callable[..., AlignmentJob]
FakeMinimap2 = Callable[..., tuple[Path, Path]]

requires_subprocess = pytest.mark.skipif(
    sys.platform in ('emscripten', 'wasi', 'win32'),
    reason='needs POSIX executables and subprocess support',
)

FAKE_MINIMAP2 = """\
#!/bin/sh
# Record the arguments, then behave like minimap2 for index/map commands.
echo "$@" >> "{log}"
[ -n "{fail}" ] && {{ echo "simulated failure" >&2; exit 3; }}
while [ $# -gt 0 ]; do
  case "$1" in
    -d) touch "$2"; shift ;;
    -o) echo "@HD	VN:1.6" > "$2"; shift ;;
  esac
  shift
done
"""


@pytest.fixture
def make_job(tmp_path: Path) -> JobFactory:
    def factory(name: str = 'reference_all_to_target_all.sam', **overrides: object) -> AlignmentJob:
        target = tmp_path / 'target.fa'
        target.touch()
        fields: dict[str, object] = {
            'features_fasta': tmp_path / 'reference_all_genes.fa',
            'target_fasta': target,
            'output_sam': tmp_path / name,
            'minimap2_options': ('-a', '--eqx'),
            'threads': 2,
            'liftover_type': LiftoverType.CHROM_BY_CHROM,
        }
        fields.update(overrides)
        return AlignmentJob(**fields)  # type: ignore[arg-type]

    return factory


@pytest.fixture
def fake_minimap2(tmp_path: Path) -> FakeMinimap2:
    """Return a factory writing a fake minimap2 script and its argument log."""

    def factory(fail: bool = False) -> tuple[Path, Path]:
        log = tmp_path / 'minimap2_calls.log'
        script = tmp_path / 'minimap2'
        script.write_text(FAKE_MINIMAP2.format(log=log, fail='yes' if fail else ''))
        script.chmod(script.stat().st_mode | stat.S_IEXEC)
        return script, log

    return factory


@requires_subprocess
class TestMinimap2Aligner:
    def test_builds_index_then_maps(
        self, fake_minimap2: FakeMinimap2, make_job: JobFactory
    ) -> None:
        script, log = fake_minimap2()
        job = make_job()

        Minimap2Aligner(str(script)).align(job)

        index, mapping = log.read_text().splitlines()
        assert index == f'-d {job.target_fasta}.mmi {job.target_fasta} -a --eqx -t 2'
        assert (
            mapping
            == f'-o {job.output_sam} {job.target_fasta}.mmi {job.features_fasta} -a --eqx -t 2'
        )

    def test_reuses_existing_index(self, fake_minimap2: FakeMinimap2, make_job: JobFactory) -> None:
        script, log = fake_minimap2()
        job = make_job()
        Path(f'{job.target_fasta}.mmi').touch()

        Minimap2Aligner(str(script)).align(job)

        assert len(log.read_text().splitlines()) == 1

    def test_large_targets_use_split_prefix(
        self, fake_minimap2: FakeMinimap2, make_job: JobFactory, tmp_path: Path
    ) -> None:
        script, log = fake_minimap2()
        job = make_job(split_prefix=tmp_path / 'split')

        Minimap2Aligner(str(script)).align(job)

        (mapping,) = log.read_text().splitlines()
        assert mapping.endswith(f'--split-prefix {tmp_path / "split"} -t 2')

    def test_returns_output_path(self, fake_minimap2: FakeMinimap2, make_job: JobFactory) -> None:
        script, _ = fake_minimap2()
        job = make_job()
        assert Minimap2Aligner(str(script)).align(job) == job.output_sam

    def test_failure_raises_alignment_error_with_stderr(
        self, fake_minimap2: FakeMinimap2, make_job: JobFactory
    ) -> None:
        script, _ = fake_minimap2(fail=True)
        with pytest.raises(AlignmentError, match='simulated failure'):
            Minimap2Aligner(str(script)).align(make_job())


def test_missing_minimap2_executable_is_reported(make_job: JobFactory) -> None:
    with pytest.raises(AlignmentError, match='not found'):
        Minimap2Aligner('no-such-minimap2-executable').align(make_job())


class TestPrecomputedSamAligner:
    def test_returns_file_with_matching_name(self, tmp_path: Path, make_job: JobFactory) -> None:
        directory = tmp_path / 'sams'
        directory.mkdir()
        (directory / 'reference_all_to_target_all.sam').touch()
        aligner = PrecomputedSamAligner(directory)
        assert aligner.align(make_job()) == directory / 'reference_all_to_target_all.sam'

    def test_missing_file_error_suggests_minimap2_command(
        self, tmp_path: Path, make_job: JobFactory
    ) -> None:
        with pytest.raises(AlignmentError, match=r'minimap2 -a --eqx .*target\.fa'):
            PrecomputedSamAligner(tmp_path).align(make_job())


class TestCallableAligner:
    def test_uses_job_output_when_callable_returns_none(self, make_job: JobFactory) -> None:
        job = make_job()

        def write_sam(job: AlignmentJob) -> None:
            job.output_sam.write_text('@HD\n')

        aligner = CallableAligner(write_sam)
        assert aligner.align(job) == job.output_sam

    def test_uses_path_returned_by_callable(self, tmp_path: Path, make_job: JobFactory) -> None:
        elsewhere = tmp_path / 'elsewhere.sam'
        elsewhere.touch()
        assert CallableAligner(lambda _job: elsewhere).align(make_job()) == elsewhere

    def test_raises_when_no_sam_is_produced(self, make_job: JobFactory) -> None:
        with pytest.raises(AlignmentError, match='did not produce'):
            CallableAligner(lambda _job: None).align(make_job())

    def test_is_not_dispatched_to_worker_processes(self) -> None:
        assert CallableAligner(lambda _job: None).parallel_safe is False


class TestPlanAlignmentJobs:
    @pytest.fixture
    def genome(self, tmp_path: Path) -> Path:
        path = tmp_path / 'target.fa'
        path.write_text('>chrA\nACGTACGT\n>chrB\nGGGGCCCC\n')
        return path

    def test_whole_genome_stage_plans_one_job(self, genome: Path, tmp_path: Path) -> None:
        jobs = runner.plan_alignment_jobs(
            ['ref.fa'],
            [str(genome)],
            'ref.fa',
            str(genome),
            tmp_path,
            ['-a'],
            4,
            LiftoverType.CHROM_BY_CHROM,
        )
        assert [(job.output_sam.name, job.target_fasta, job.threads) for job in jobs] == [
            ('reference_all_to_target_all.sam', genome, 4)
        ]

    def test_chromosome_stage_plans_one_job_per_chromosome(
        self, genome: Path, tmp_path: Path
    ) -> None:
        jobs = runner.plan_alignment_jobs(
            ['refA', 'refB'],
            ['chrA', 'chrB'],
            'ref.fa',
            str(genome),
            tmp_path,
            ['-a'],
            4,
            LiftoverType.CHROM_BY_CHROM,
        )
        assert [(job.output_sam.name, job.target_fasta.name, job.threads) for job in jobs] == [
            ('refA_to_chrA.sam', 'chrA.fa', 2),
            ('refB_to_chrB.sam', 'chrB.fa', 2),
        ]

    def test_chromosome_stage_writes_target_chromosomes(self, genome: Path, tmp_path: Path) -> None:
        runner.plan_alignment_jobs(
            ['refA'],
            ['chrB'],
            'ref.fa',
            str(genome),
            tmp_path,
            ['-a'],
            1,
            LiftoverType.CHROM_BY_CHROM,
        )
        assert (tmp_path / 'chrB.fa').read_text() == '>chrB\nGGGGCCCC'

    def test_large_genomes_use_split_prefix(
        self, genome: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        monkeypatch.setattr(runner, 'MAX_SINGLE_INDEX_SIZE', 10)
        (job,) = runner.plan_alignment_jobs(
            ['ref.fa'],
            [str(genome)],
            'ref.fa',
            str(genome),
            tmp_path,
            ['-a'],
            1,
            LiftoverType.COPIES,
        )
        assert job.split_prefix == tmp_path / 'reference_all_copies_to_target_all_split'


class TestRunAlignmentJobs:
    @pytest.fixture
    def sam_directory(self, tmp_path: Path) -> Path:
        directory = tmp_path / 'sams'
        directory.mkdir()
        for name in ('a.sam', 'b.sam', 'c.sam'):
            (directory / name).touch()
        return directory

    @pytest.mark.parametrize(
        'processes', [pytest.param(1, id='serial'), pytest.param(3, id='parallel')]
    )
    def test_results_follow_job_order(
        self, sam_directory: Path, make_job: JobFactory, processes: int
    ) -> None:
        jobs = [make_job(name) for name in ('c.sam', 'a.sam', 'b.sam')]
        results = runner.run_alignment_jobs(PrecomputedSamAligner(sam_directory), jobs, processes)
        assert [path.name for path in results] == ['c.sam', 'a.sam', 'b.sam']

    def test_unsafe_aligners_run_in_process(self, make_job: JobFactory) -> None:
        seen: list[str] = []

        def record(job: AlignmentJob) -> None:
            seen.append(job.output_sam.name)
            job.output_sam.touch()

        jobs = [make_job('x.sam'), make_job('y.sam')]
        runner.run_alignment_jobs(CallableAligner(record), jobs, processes=4)
        assert seen == ['x.sam', 'y.sam']

    def test_webassembly_platforms_disable_multiprocessing(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        monkeypatch.setattr(sys, 'platform', 'emscripten')
        assert runner.can_use_multiprocessing() is False
