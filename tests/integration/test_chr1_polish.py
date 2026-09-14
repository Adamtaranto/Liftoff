"""Polishing of lifted chromosome I genes with broken coding sequences."""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, field
import gzip
import json
from pathlib import Path
from typing import Any

import pytest

from liftoff import cds, polish
from liftoff.config import LiftoffConfig
from liftoff.pipeline import LiftoffResult, run_liftoff

from .helpers import Chr1Inputs, RunPaths


@pytest.fixture(
    scope='module',
    params=[
        pytest.param('minimap2', marks=pytest.mark.minimap2),
        pytest.param('precomputed'),
    ],
)
def polish_run(
    request: pytest.FixtureRequest,
    tmp_path_factory: pytest.TempPathFactory,
    copy_chr1_inputs: Callable[[Path], Chr1Inputs],
    chr1_alignments: Path,
    run_cli: Callable[..., RunPaths],
) -> RunPaths:
    """Run ``liftoff --polish`` once per alignment backend."""
    workdir = tmp_path_factory.mktemp(f'chr1_polish_{request.param}')
    inputs = copy_chr1_inputs(workdir / 'inputs')
    extra = ['--polish', '--threads', '2']
    if request.param == 'precomputed':
        extra += ['--alignments', str(chr1_alignments)]
    return run_cli(workdir, inputs.target, inputs.reference, inputs.gff, extra)


def test_unpolished_annotation_matches_expected(
    polish_run: RunPaths,
    expected_dir: Path,
    assert_same_annotation: Callable[[Path, Path], None],
) -> None:
    assert_same_annotation(polish_run.output, expected_dir / 'chr1_polish.gff3')


def test_polished_annotation_matches_expected(
    polish_run: RunPaths,
    expected_dir: Path,
    assert_same_annotation: Callable[[Path, Path], None],
) -> None:
    assert_same_annotation(polish_run.polished, expected_dir / 'chr1_polish_polished.gff3')


def test_polish_writes_a_sam_file_for_polished_genes(polish_run: RunPaths) -> None:
    assert (polish_run.intermediate / polish.POLISH_SAM_NAME).is_file()


# ---------------------------------------------------------------------------
# Parity with the parasail-based Liftoff 1.6.3
# ---------------------------------------------------------------------------
@dataclass
class PolishRecording:
    """Observations made while polishing."""

    sam_files: list[dict[str, str]] = field(default_factory=list)
    check_cds_calls: list[list[list[Any]]] = field(default_factory=list)


@dataclass(frozen=True)
class TracedRun:
    result: LiftoffResult
    config: LiftoffConfig
    recording: PolishRecording


@pytest.fixture(scope='module')
def parasail_trace(parasail_dir: Path) -> dict[str, Any]:
    """Polishing trace recorded with Liftoff 1.6.3 on the chromosome I data."""
    with gzip.open(parasail_dir / 'polish_trace_chr1.json.gz', 'rt') as handle:
        trace: dict[str, Any] = json.load(handle)
    return trace


@pytest.fixture(scope='module')
def traced_polish_run(
    tmp_path_factory: pytest.TempPathFactory,
    copy_chr1_inputs: Callable[[Path], Chr1Inputs],
    chr1_alignments: Path,
) -> TracedRun:
    """Polish with ORF trimming disabled, recording polish SAMs and CDS checks.

    Liftoff 1.6.3 never trimmed broken CDSs to their longest ORF, so disabling
    trimming reproduces its selection of genes to polish (23 genes). Every exon
    alignment is then comparable with the parasail-era trace.
    """
    workdir = tmp_path_factory.mktemp('chr1_polish_trace')
    inputs = copy_chr1_inputs(workdir / 'inputs')
    paths = RunPaths(workdir)
    recording = PolishRecording()
    original_polish = polish.polish_annotations
    original_check_cds = cds.check_cds

    def record_polish(*args: Any, **kwargs: Any) -> bool:
        polished = original_polish(*args, **kwargs)
        if polished:
            intermediate_dir, target_feature = args[3], args[5]
            sam_text = (intermediate_dir / polish.POLISH_SAM_NAME).read_text()
            recording.sam_files.append({'feature': target_feature, 'sam': sam_text})
        return polished

    def record_check_cds(feature_list: Any, *args: Any, **kwargs: Any) -> None:
        original_check_cds(feature_list, *args, **kwargs)
        recording.check_cds_calls.append(
            [
                [key, f.id, f.featuretype, f.seqid, f.strand, f.start, f.end, f.frame, f.attributes]
                for key in sorted(feature_list)
                for f in feature_list[key]
            ]
        )

    config = LiftoffConfig(
        target=str(inputs.target),
        reference=str(inputs.reference),
        gff=str(inputs.gff),
        output=str(paths.output),
        unmapped=str(paths.unmapped),
        intermediate_dir=str(paths.intermediate),
        alignments=str(chr1_alignments),
        polish=True,
    )
    with pytest.MonkeyPatch.context() as patch:
        patch.setattr(polish, 'polish_annotations', record_polish)
        patch.setattr(cds, 'check_cds', record_check_cds)
        patch.setattr(cds, 'find_longest_orf', lambda *_args, **_kwargs: None)
        result = run_liftoff(config)
    return TracedRun(result, config, recording)


def test_traced_run_writes_unpolished_and_polished_outputs(traced_polish_run: TracedRun) -> None:
    output = traced_polish_run.config.output
    assert traced_polish_run.result.output_paths == [output, output + '_polished']


def test_polish_sam_records_match_parasail(
    traced_polish_run: TracedRun, parasail_trace: dict[str, Any]
) -> None:
    assert traced_polish_run.recording.sam_files == parasail_trace['polish_sams']


def test_polished_candidate_genes_match_parasail(
    traced_polish_run: TracedRun, parasail_trace: dict[str, Any]
) -> None:
    # The final check_cds call annotates the polished candidate genes.
    candidates = json.loads(json.dumps(traced_polish_run.recording.check_cds_calls[-1]))
    assert candidates == parasail_trace['check_cds_calls'][-1]


def test_polished_annotation_matches_liftoff_1_6_3(
    traced_polish_run: TracedRun,
    expected_dir: Path,
    assert_same_annotation: Callable[[Path, Path], None],
) -> None:
    assert_same_annotation(
        Path(traced_polish_run.config.output + '_polished'),
        expected_dir / 'chr1_polish_polished_without_orf_trimming.gff3',
    )
