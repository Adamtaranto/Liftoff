"""Load reference annotations and split them into a feature hierarchy.

Annotations are parsed with :mod:`gffutils` into an SQLite feature database.
Liftoff then classifies every feature as a top-level *parent* (e.g. gene), an
*intermediate* (e.g. transcript) or a lowest-level *child* (e.g. exon, CDS).
"""

from __future__ import annotations

from collections.abc import Iterator
from contextlib import contextmanager
import gzip
import json
import logging
from pathlib import Path
import sqlite3
from typing import IO, Any

import gffutils

from liftoff.errors import DuplicateFeatureIdError, GffSyntaxError, InputError
from liftoff.log import get_logger
from liftoff.models import Feature, FeatureHierarchy
from liftoff.utils import ParentOrder, find_parent_order

logger = get_logger(__name__)

#: Column list used for every query that materialises a :class:`Feature`.
_FEATURE_COLUMNS = 'id, seqid, source, featuretype, start, end, strand, frame, attributes'


def build_database(
    gff_file: str | None,
    db_file: str | None,
    infer_genes: bool = False,
    infer_transcripts: bool = False,
) -> gffutils.FeatureDB:
    """Create (or open) the gffutils feature database.

    Parameters
    ----------
    gff_file : str or None
        Annotation to parse. The database is written next to it as
        ``<gff_file>_db``, replacing any existing file.
    db_file : str or None
        Existing database to open instead of parsing ``gff_file``.
    infer_genes : bool, default False
        Ask gffutils to create gene features from transcripts.
    infer_transcripts : bool, default False
        Ask gffutils to create transcript features from exons/CDS.

    Returns
    -------
    gffutils.FeatureDB
        Connection to the feature database.

    Raises
    ------
    InputError
        If neither input is provided or a file does not exist.
    GffSyntaxError
        If the annotation cannot be parsed.
    DuplicateFeatureIdError
        If distinct features in the annotation share an ``ID``.
    """
    # Preserve percent-encoded characters (e.g. "%3B") exactly as written.
    gffutils.constants.ignore_url_escape_characters = True
    if db_file is not None:
        if not Path(db_file).is_file():
            raise InputError(f'feature database not found: {db_file}')
        return gffutils.FeatureDB(db_file)
    if gff_file is None:
        raise InputError('an annotation file or feature database is required')
    if not Path(gff_file).is_file():
        raise InputError(f'annotation file not found: {gff_file}')
    validate_unique_ids(gff_file)
    logger.info('building feature database from %s', gff_file)
    try:
        return gffutils.create_db(
            gff_file,
            gff_file + '_db',
            merge_strategy='create_unique',
            force=True,
            disable_infer_transcripts=not infer_transcripts,
            disable_infer_genes=not infer_genes,
            verbose=logger.isEnabledFor(logging.DEBUG),
        )
    except Exception as exc:
        line_number = find_problem_line(gff_file)
        if line_number is None:
            raise GffSyntaxError(None, f'could not parse annotation {gff_file}: {exc}') from exc
        raise GffSyntaxError(
            line_number, f'incorrect GFF/GTF syntax on line {line_number} of {gff_file}'
        ) from exc


#: Number of conflicting IDs listed in a :class:`DuplicateFeatureIdError` message.
MAX_REPORTED_DUPLICATES = 10


@contextmanager
def _open_annotation(path: str) -> Iterator[IO[str]]:
    """Open a plain or gzip-compressed annotation file as text."""
    with Path(path).open('rb') as probe:
        compressed = probe.read(2) == b'\x1f\x8b'
    with gzip.open(path, 'rt') if compressed else Path(path).open() as handle:
        yield handle


def _gff3_attributes(column: str) -> dict[str, str]:
    """Parse GFF3 ``key=value`` attributes; GTF-style attributes yield nothing."""
    attributes = {}
    for item in column.split(';'):
        key, separator, value = item.partition('=')
        if separator:
            attributes[key.strip()] = value.strip()
    return attributes


def validate_unique_ids(gff_file: str) -> None:
    """Check that distinct features in a GFF3 annotation have distinct ``ID`` values.

    GFF3 allows a single feature to span several lines that share an ``ID``
    (for example a CDS split across exons). Lines sharing an ``ID`` are
    therefore accepted when they agree on sequence, type, strand and
    ``Parent``; any other repeated ``ID`` identifies two different features.
    Such annotations are rejected because gffutils would silently rename one
    of the features, detaching it from the feature it belongs to. GTF files
    carry no ``ID`` attribute and are not checked.

    Parameters
    ----------
    gff_file : str
        Plain or gzip-compressed GFF3/GTF annotation.

    Raises
    ------
    DuplicateFeatureIdError
        If any ``ID`` is used by more than one distinct feature.
    """
    first_seen: dict[str, tuple[int, tuple[str, str, str, str]]] = {}
    duplicates: list[tuple[str, int, int]] = []
    with _open_annotation(gff_file) as handle:
        for line_number, line in enumerate(handle, start=1):
            if line.startswith('##FASTA'):
                break
            if line.startswith('#') or not line.strip():
                continue
            fields = line.rstrip('\r\n').split('\t')
            if len(fields) < 9:
                continue  # malformed lines are reported by the GFF parser
            attributes = _gff3_attributes(fields[8])
            feature_id = attributes.get('ID')
            if feature_id is None:
                continue
            feature = (fields[0], fields[2], fields[6], attributes.get('Parent', ''))
            previous = first_seen.setdefault(feature_id, (line_number, feature))
            if previous[1] != feature:
                duplicates.append((feature_id, previous[0], line_number))
    if duplicates:
        examples = '; '.join(
            f"ID '{feature_id}' on lines {first} and {second}"
            for feature_id, first, second in duplicates[:MAX_REPORTED_DUPLICATES]
        )
        more = len(duplicates) - MAX_REPORTED_DUPLICATES
        suffix = f'; and {more} more' if more > 0 else ''
        raise DuplicateFeatureIdError(
            duplicates,
            f'{gff_file}: {len(duplicates)} feature(s) reuse the ID of a different feature. '
            'IDs must be unique; lines may only share an ID when they are parts of the same '
            f'feature (same sequence, type, strand and Parent). {examples}{suffix}',
        )


def find_problem_line(gff_file: str) -> int | None:
    """Locate the first annotation line that gffutils cannot parse.

    Parameters
    ----------
    gff_file : str
        Plain-text annotation file.

    Returns
    -------
    int or None
        One-based line number, or ``None`` if every line parses on its own.
    """
    try:
        with Path(gff_file).open() as handle:
            lines = handle.readlines()
    except (OSError, UnicodeDecodeError):
        return None
    for index, line in enumerate(lines, start=1):
        if line.startswith('#'):
            continue
        try:
            gffutils.create_db(line, ':memory:', from_string=True, force=True).conn.close()
        except Exception:
            return index
    return None


def _feature_from_row(row: sqlite3.Row | tuple[Any, ...]) -> Feature:
    """Build a :class:`Feature` from a row selected with ``_FEATURE_COLUMNS``."""
    feature_id, seqid, source, featuretype, start, end, strand, frame, attributes = row
    return Feature(
        id=feature_id,
        featuretype=featuretype,
        seqid=seqid,
        source=source,
        strand=strand,
        start=start,
        end=end,
        frame=frame,
        attributes=json.loads(attributes),
    )


def separate_parents_and_children(
    feature_db: gffutils.FeatureDB, parent_types_to_lift: list[str] | None
) -> tuple[FeatureHierarchy, ParentOrder]:
    """Classify reference features into parents, intermediates and children.

    Rows are read in an explicit order — database insertion order (``rowid``)
    for features, and ``(child rowid, relation rowid)`` for children — which
    determines the order of features in the output annotation.

    Parameters
    ----------
    feature_db : gffutils.FeatureDB
        Reference feature database.
    parent_types_to_lift : list of str or None
        Feature types eligible to be lifted as top-level features; ``None``
        lifts top-level features of every type.

    Returns
    -------
    hierarchy : FeatureHierarchy
        Classified reference features.
    parent_order : ParentOrder
        Top-level features sorted by position.
    """
    cursor = feature_db.conn.cursor()
    relations = [
        (parent, child)
        for parent, child in cursor.execute(
            """SELECT r.parent, r.child FROM relations AS r
               JOIN features AS a ON a.id = r.parent
               JOIN features AS b ON b.id = r.child"""
        )
        if parent != child
    ]
    all_ids = {row[0] for row in cursor.execute('SELECT id FROM features')}
    child_ids = {child for _, child in relations}
    parent_ids = {parent for parent, _ in relations}

    # Features that are never a parent are leaves; those never a child are roots.
    lowest_children = all_ids - parent_ids
    highest_parents = all_ids - child_ids
    intermediates = child_ids & parent_ids

    parent_dict: dict[str, Feature] = {}
    child_dict: dict[str, list[Feature]] = {}
    intermediate_dict: dict[str, Feature] = {}
    _add_parents(cursor, parent_dict, child_dict, highest_parents, parent_types_to_lift)
    _add_children(cursor, feature_db, parent_dict, child_dict, lowest_children)
    _add_intermediates(cursor, feature_db, intermediate_dict, intermediates)
    parent_order = find_parent_order(list(parent_dict.values()))
    return FeatureHierarchy(parent_dict, intermediate_dict, child_dict), parent_order


def _add_parents(
    cursor: sqlite3.Cursor,
    parent_dict: dict[str, Feature],
    child_dict: dict[str, list[Feature]],
    highest_parents: set[str],
    parent_types_to_lift: list[str] | None,
) -> None:
    """Register root features whose type should be lifted."""
    for row in cursor.execute(f'SELECT {_FEATURE_COLUMNS} FROM features ORDER BY rowid').fetchall():
        if row[0] not in highest_parents:
            continue
        parent = _feature_from_row(row)
        if parent_types_to_lift is None or parent.featuretype in parent_types_to_lift:
            parent_dict[parent.id] = parent
            child_dict[parent.id] = []


def _add_children(
    cursor: sqlite3.Cursor,
    feature_db: gffutils.FeatureDB,
    parent_dict: dict[str, Feature],
    child_dict: dict[str, list[Feature]],
    lowest_children: set[str],
) -> None:
    """Attach leaf features to their top-level parents."""
    columns = ', '.join(f'f.{column.strip()}' for column in _FEATURE_COLUMNS.split(','))
    rows = cursor.execute(
        f"""SELECT r.parent, {columns} FROM relations AS r
            JOIN features AS f ON f.id = r.child
            ORDER BY f.rowid, r.rowid"""
    ).fetchall()
    added_children_ids: set[str] = set()
    for parent_id, *feature_row in rows:
        # Only relations to a lifted top-level feature are kept; relations to
        # intermediates (e.g. exon -> mRNA) are represented via "Parent".
        if feature_row[0] not in lowest_children or parent_id not in parent_dict:
            continue
        child = _feature_from_row(tuple(feature_row))
        if child.featuretype == 'intron':
            continue
        if 'Parent' not in child.attributes:
            _add_parent_tag(child, feature_db)
        child_dict[parent_id].append(child)
        added_children_ids.add(child.id)
    # Single-level features (e.g. a gene without transcripts) are their own child.
    for feature_id in sorted(lowest_children - added_children_ids):
        if feature_id in parent_dict:
            child_dict[feature_id] = [parent_dict[feature_id]]


def _add_parent_tag(feature: Feature, feature_db: gffutils.FeatureDB) -> None:
    """Set a ``Parent`` attribute from the database relations."""
    parent_id = ''
    parents = [p for p in feature_db.parents(feature.id, level=1) if feature.id != p.id]
    if not parents:
        parents = [p for p in feature_db.parents(feature.id) if feature.id != p.id]
    if parents:
        parent_id = parents[0].id
    feature.attributes['Parent'] = [parent_id]


def _add_intermediates(
    cursor: sqlite3.Cursor,
    feature_db: gffutils.FeatureDB,
    intermediate_dict: dict[str, Feature],
    intermediate_ids: set[str],
) -> None:
    """Register features that are both parents and children."""
    for row in cursor.execute(f'SELECT {_FEATURE_COLUMNS} FROM features ORDER BY rowid').fetchall():
        if row[0] not in intermediate_ids:
            continue
        feature = _feature_from_row(row)
        intermediate_dict[feature.id] = feature
        if 'Parent' not in feature.attributes:
            _add_parent_tag(feature, feature_db)


def get_feature_order(feature_db: gffutils.FeatureDB) -> dict[str, int]:
    """Rank feature types for sorting lifted child features.

    Exons sort first, then CDS, then all other types in database order.

    Parameters
    ----------
    feature_db : gffutils.FeatureDB
        Reference feature database.

    Returns
    -------
    dict of str to int
        Sort rank for each feature type.
    """
    feature_types = list(feature_db.featuretypes())
    feature_order: dict[str, int] = {}
    for preferred in ('exon', 'CDS'):
        if preferred in feature_types:
            feature_order[preferred] = len(feature_order)
    for feature_type in feature_types:
        if feature_type not in feature_order:
            feature_order[feature_type] = len(feature_order)
    return feature_order


def annotation_format(feature_db: gffutils.FeatureDB) -> str:
    """Return the dialect of the parsed annotation.

    Parameters
    ----------
    feature_db : gffutils.FeatureDB
        Reference feature database.

    Returns
    -------
    str
        ``"gff3"`` or ``"gtf"``.
    """
    return str(feature_db.dialect['fmt'])
