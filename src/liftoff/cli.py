"""Command line interface for Liftoff."""

from __future__ import annotations

import argparse
from collections.abc import Sequence
import logging
import os
import shlex
import sys

from liftoff import __version__
from liftoff.config import DEFAULT_MM2_OPTIONS, STDOUT, LiftoffConfig
from liftoff.errors import ConfigError, LiftoffError
from liftoff.log import configure_logging, get_logger

logger = get_logger(__name__)


class _HelpFormatter(argparse.ArgumentDefaultsHelpFormatter):
    """Show defaults only for options that have a meaningful default value."""

    def _get_help_string(self, action: argparse.Action) -> str | None:
        # Identity checks, so numeric defaults such as 0.0 are still shown.
        if action.default is None or action.default is False or action.default is argparse.SUPPRESS:
            return action.help
        return super()._get_help_string(action)


def build_parser() -> argparse.ArgumentParser:
    """Create the argument parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the ``liftoff`` command.
    """
    parser = argparse.ArgumentParser(
        prog='liftoff',
        description='Lift features from one genome assembly to another.',
        formatter_class=_HelpFormatter,
    )
    parser.add_argument('-V', '--version', action='version', version=f'%(prog)s {__version__}')

    sequences = parser.add_argument_group('Required input (sequences)')
    sequences.add_argument('target', help='target FASTA genome to lift features to')
    sequences.add_argument('reference', help='reference FASTA genome to lift features from')

    annotation = parser.add_argument_group('Required input (annotation)')
    source = annotation.add_mutually_exclusive_group(required=True)
    source.add_argument(
        '-g', '--gff', metavar='GFF', help='annotation file to lift over in GFF3 or GTF format'
    )
    source.add_argument(
        '--db',
        metavar='DB',
        help='feature database built by a previous run (written next to the GFF as <GFF>_db)',
    )

    output = parser.add_argument_group('Output')
    output.add_argument(
        '-o',
        '--output',
        default=STDOUT,
        metavar='FILE',
        help="write lifted annotation to FILE in the input format; 'stdout' writes to the terminal",
    )
    output.add_argument(
        '-u',
        '--unmapped',
        default='unmapped_features.txt',
        metavar='FILE',
        help='write IDs of unmapped features to FILE',
    )
    output.add_argument(
        '--exclude-partial',
        action='store_true',
        help='report mappings below --min-coverage/--min-identity as unmapped instead of '
        'writing them with partial_mapping=True / low_identity=True',
    )
    output.add_argument(
        '--intermediate-dir',
        default='intermediate_files',
        metavar='DIR',
        help='directory for intermediate FASTA and SAM files',
    )

    alignment = parser.add_argument_group('Alignment')
    alignment.add_argument(
        '--mm2-options',
        default=DEFAULT_MM2_OPTIONS,
        metavar='STR',
        help="space-delimited minimap2 options; use --mm2-options='-r 2k -z 5000' so values "
        "starting with '-' are not parsed as Liftoff options",
    )
    alignment.add_argument(
        '--minimap2', metavar='PATH', help='minimap2 executable (default: search PATH)'
    )
    alignment.add_argument(
        '--alignments',
        metavar='DIR',
        help='use pre-computed SAM files from DIR instead of running minimap2',
    )
    alignment.add_argument(
        '-a',
        '--min-coverage',
        default=0.5,
        type=float,
        metavar='A',
        help='designate a feature mapped only if its child features align with coverage >= A',
    )
    alignment.add_argument(
        '-s',
        '--min-identity',
        default=0.5,
        type=float,
        metavar='S',
        help='designate a feature mapped only if its child features align with identity >= S',
    )
    alignment.add_argument(
        '-d',
        '--distance-factor',
        default=2.0,
        type=float,
        metavar='D',
        help='alignment blocks further apart in the target than D times their distance in the '
        'reference are not chained',
    )
    alignment.add_argument(
        '--flank',
        default=0.0,
        type=float,
        metavar='F',
        help='flanking sequence to align, as a fraction [0-1] of gene length',
    )

    misc = parser.add_argument_group('Lift-over settings')
    misc.add_argument(
        '-p', '--threads', default=1, type=int, metavar='P', help='parallel alignment processes'
    )
    feature_selection = misc.add_mutually_exclusive_group()
    feature_selection.add_argument(
        '-f',
        '--feature-types',
        metavar='FILE',
        help='file listing additional top-level feature types to lift (one per line)',
    )
    feature_selection.add_argument(
        '--all-feature-types',
        action='store_true',
        help='lift every top-level feature type in the annotation, not only genes',
    )
    misc.add_argument(
        '--infer-genes',
        action='store_true',
        help='annotation only contains transcripts and exon/CDS features',
    )
    misc.add_argument(
        '--infer-transcripts',
        action='store_true',
        help='annotation only contains genes and exon/CDS features',
    )
    misc.add_argument(
        '--chroms',
        metavar='FILE',
        help='comma separated file of corresponding reference,target chromosomes',
    )
    misc.add_argument(
        '--unplaced',
        metavar='FILE',
        help='file of unplaced reference sequence names to lift after --chroms',
    )
    misc.add_argument(
        '--copies', action='store_true', help='look for extra gene copies in the target genome'
    )
    misc.add_argument(
        '--copy-identity',
        default=1.0,
        type=float,
        metavar='SC',
        help='with --copies, minimum exon/CDS identity for an extra copy (>= --min-identity)',
    )
    misc.add_argument(
        '--max-overlap',
        default=0.1,
        type=float,
        metavar='O',
        help='maximum fraction [0-1] of overlap allowed between two features',
    )
    misc.add_argument(
        '--mismatch', default=2, type=int, metavar='M', help='mismatch penalty in exons'
    )
    misc.add_argument(
        '--gap-open', default=2, type=int, metavar='GO', help='gap open penalty in exons'
    )
    misc.add_argument(
        '--gap-extend', default=1, type=int, metavar='GE', help='gap extend penalty in exons'
    )
    misc.add_argument(
        '--polish',
        action='store_true',
        help='re-align exons to repair broken CDS; writes <output> and <output>_polished',
    )
    misc.add_argument(
        '--cds',
        action=argparse.BooleanOptionalAction,
        default=True,
        help='annotate CDS status (partial, missing start/stop, in-frame stop codon)',
    )

    logging_group = parser.add_argument_group('Logging')
    verbosity = logging_group.add_mutually_exclusive_group()
    verbosity.add_argument(
        '-v', '--verbose', action='store_true', help='show debug messages and minimap2 output'
    )
    verbosity.add_argument(
        '-q', '--quiet', action='store_true', help='only show warnings and errors'
    )
    return parser


def config_from_args(args: argparse.Namespace, argv: Sequence[str]) -> LiftoffConfig:
    """Build a :class:`LiftoffConfig` from parsed arguments.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command line.
    argv : sequence of str
        Raw arguments, recorded in the output header.

    Returns
    -------
    LiftoffConfig
        Validated configuration.
    """
    return LiftoffConfig(
        target=args.target,
        reference=args.reference,
        gff=args.gff,
        db=args.db,
        output=args.output,
        unmapped=args.unmapped,
        exclude_partial=args.exclude_partial,
        intermediate_dir=args.intermediate_dir,
        mm2_options=args.mm2_options,
        min_coverage=args.min_coverage,
        min_identity=args.min_identity,
        distance_factor=args.distance_factor,
        flank=args.flank,
        threads=args.threads,
        minimap2=args.minimap2,
        feature_types=args.feature_types,
        all_feature_types=args.all_feature_types,
        infer_genes=args.infer_genes,
        infer_transcripts=args.infer_transcripts,
        chroms=args.chroms,
        unplaced=args.unplaced,
        copies=args.copies,
        copy_identity=args.copy_identity,
        max_overlap=args.max_overlap,
        mismatch=args.mismatch,
        gap_open=args.gap_open,
        gap_extend=args.gap_extend,
        polish=args.polish,
        cds=args.cds,
        alignments=args.alignments,
        command_line=shlex.join(['liftoff', *argv]),
    )


def main(argv: Sequence[str] | None = None) -> int:
    """Run Liftoff from the command line.

    Parameters
    ----------
    argv : sequence of str or None, default None
        Arguments excluding the program name; defaults to ``sys.argv[1:]``.

    Returns
    -------
    int
        Process exit status.
    """
    argv = list(sys.argv[1:] if argv is None else argv)
    parser = build_parser()
    args = parser.parse_args(argv)
    configure_logging(
        logging.DEBUG if args.verbose else logging.WARNING if args.quiet else logging.INFO
    )
    try:
        config = config_from_args(args, argv)
    except ConfigError as exc:
        parser.error(str(exc))

    # Imported lazily so that `liftoff --help` stays fast.
    from liftoff.pipeline import run_liftoff

    try:
        run_liftoff(config)
    except LiftoffError as exc:
        logger.error('%s', exc)
        return 1
    except BrokenPipeError:
        # Output piped into a command that exited early (e.g. `liftoff ... | head`).
        _silence_stdout()
        return 1
    return 0


def _silence_stdout() -> None:
    """Point stdout at the null device so the interpreter's final flush cannot fail."""
    try:
        devnull = os.open(os.devnull, os.O_WRONLY)
        os.dup2(devnull, sys.stdout.fileno())
    except (OSError, ValueError, AttributeError):  # pragma: no cover - stdout without a fd
        pass
