"""Write lifted features as GFF3 or GTF."""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from contextlib import ExitStack
from pathlib import Path
import sys
from typing import TextIO

from liftoff import __version__
from liftoff.config import STDOUT
from liftoff.models import Attributes, Feature
from liftoff.utils import LiftedFeatures, get_parent_list

#: Attributes that are re-appended (so that they appear last) on output.
_KEYS_TO_READD_AT_END = ('copy_num_ID', 'partial_mapping', 'low_identity', 'extra_copy_number')


def write_header(handle: TextIO, out_type: str, command_line: str) -> None:
    """Write the annotation header lines.

    Parameters
    ----------
    handle : TextIO
        Destination stream.
    out_type : str
        ``"gff3"`` or ``"gtf"``.
    command_line : str
        Command recorded in the header.
    """
    if out_type == 'gff3':
        handle.write('##gff-version 3\n')
    handle.write(f'# Liftoff v{__version__}\n')
    handle.write(f'# {command_line}\n')


def write_new_gff(
    lifted_features: LiftedFeatures,
    output: str,
    out_type: str,
    min_coverage: float,
    min_identity: float,
    command_line: str = '',
) -> None:
    """Write all lifted features, grouped by top-level feature.

    Top-level features are assigned copy numbers in ID order, then written
    sorted by position, each followed by its descendants depth-first.

    Parameters
    ----------
    lifted_features : dict of str to list of Feature
        Lifted feature groups keyed by copy ID.
    output : str
        Output path, or ``"stdout"``.
    out_type : str
        ``"gff3"`` or ``"gtf"``.
    min_coverage : float
        Features below this coverage are tagged ``partial_mapping=True``.
    min_identity : float
        Features below this identity are tagged ``low_identity=True``.
    command_line : str, default ""
        Command recorded in the header (defaults to ``sys.argv``).
    """
    with ExitStack() as stack:
        if output == STDOUT:
            handle: TextIO = sys.stdout
        else:
            handle = stack.enter_context(Path(output).open('w'))
        write_header(handle, out_type, command_line or ' '.join(sys.argv))
        parents = get_parent_list(lifted_features)
        parents.sort(key=lambda x: x.id)
        final_parent_list = finalize_parent_features(parents, min_coverage, min_identity)
        final_parent_list.sort(key=lambda x: (x.seqid, x.start))
        for final_parent in final_parent_list:
            child_features = lifted_features[final_parent.attributes['copy_id'][0]]
            parent_child_dict = build_parent_dict(child_features, final_parent)
            write_feature([final_parent], handle, parent_child_dict, out_type)


def finalize_parent_features(
    parents: list[Feature], min_coverage: float, min_identity: float
) -> list[Feature]:
    """Assign copy numbers and add summary attributes to top-level features.

    Parameters
    ----------
    parents : list of Feature
        Top-level features sorted by ID; copies of the same feature are
        numbered in this order starting from zero.
    min_coverage : float
        Coverage threshold for ``partial_mapping``.
    min_identity : float
        Identity threshold for ``low_identity``.

    Returns
    -------
    list of Feature
        The same features, updated in place.
    """
    copy_num_dict: dict[str, int] = {}
    for parent in parents:
        copy_num_dict[parent.id] = copy_num_dict[parent.id] + 1 if parent.id in copy_num_dict else 0
        _add_attributes(parent, copy_num_dict[parent.id], min_coverage, min_identity)
    return list(parents)


def _add_attributes(
    parent: Feature, copy_num: int, min_coverage: float, min_identity: float
) -> None:
    """Attach copy number and mapping-quality attributes to a parent."""
    attributes = parent.attributes
    if 'copy_id' not in attributes:
        # Remember the internal copy ID used as key into the lifted features.
        attributes['copy_id'] = attributes['copy_num_ID']
    for key in _KEYS_TO_READD_AT_END:
        attributes.pop(key, None)
    attributes['extra_copy_number'] = [str(copy_num)]
    attributes['copy_num_ID'] = [f'{parent.id}_{copy_num}']
    if float(attributes['coverage'][0]) < min_coverage:
        attributes['partial_mapping'] = ['True']
    if float(attributes['sequence_ID'][0]) < min_identity:
        attributes['low_identity'] = ['True']


def build_parent_dict(
    child_features: Iterable[Feature], final_parent: Feature
) -> dict[str, list[Feature]]:
    """Group features by their ``Parent`` attribute.

    Also propagates the parent's ``extra_copy_number`` to every child.

    Parameters
    ----------
    child_features : iterable of Feature
        All features of one lifted gene.
    final_parent : Feature
        The gene's top-level feature.

    Returns
    -------
    dict of str to list of Feature
        Children keyed by parent ID, in input order.
    """
    parent_child_dict: dict[str, list[Feature]] = {}
    for child in child_features:
        if 'Parent' in child.attributes:
            child.attributes['extra_copy_number'] = final_parent.attributes['extra_copy_number']
            parent_child_dict.setdefault(child.attributes['Parent'][0], []).append(child)
    return parent_child_dict


def write_feature(
    features: Iterable[Feature],
    handle: TextIO,
    parent_dict: Mapping[str, list[Feature]],
    output_type: str,
) -> None:
    """Write features followed (recursively) by their children.

    Parameters
    ----------
    features : iterable of Feature
        Features to write at this level.
    handle : TextIO
        Destination stream.
    parent_dict : mapping of str to list of Feature
        Children keyed by parent ID.
    output_type : str
        ``"gff3"`` or ``"gtf"``.
    """
    for feature in features:
        write_line(feature, handle, output_type)
        if feature.id in parent_dict:
            write_feature(parent_dict[feature.id], handle, parent_dict, output_type)


def write_line(feature: Feature, handle: TextIO, output_type: str) -> None:
    """Write a single feature line.

    Parameters
    ----------
    feature : Feature
        Feature to write.
    handle : TextIO
        Destination stream.
    output_type : str
        ``"gff3"`` or ``"gtf"``.
    """
    if feature.attributes['extra_copy_number'][0] != '0':
        attr_dict = edit_copy_ids(feature)
    else:
        attr_dict = feature.attributes
    line = (
        make_gff_line(attr_dict, feature)
        if output_type == 'gff3'
        else make_gtf_line(attr_dict, feature)
    )
    handle.write(line + '\n')


def make_gff_line(attr_dict: Attributes, feature: Feature) -> str:
    """Format a GFF3 line with ``ID`` as the first attribute.

    Parameters
    ----------
    attr_dict : dict of str to list of str
        Attributes to write.
    feature : Feature
        Feature providing the other columns.

    Returns
    -------
    str
        The line without a trailing newline.
    """
    attributes = [f'ID={attr_dict["ID"][0]}']
    attributes.extend(
        f'{key}={",".join(values)}'
        for key, values in attr_dict.items()
        if key not in ('copy_id', 'ID')
    )
    columns = (
        feature.seqid,
        feature.source,
        feature.featuretype,
        str(feature.start),
        str(feature.end),
        '.',
        feature.strand,
        feature.frame,
        ';'.join(attributes),
    )
    return '\t'.join(columns)


def edit_copy_ids(feature: Feature) -> Attributes:
    """Append the copy number to identifier attributes of an extra copy.

    Parameters
    ----------
    feature : Feature
        Feature of an extra gene copy.

    Returns
    -------
    dict of str to list of str
        A copy of the attributes with ``ID``, ``Parent`` and ``*_id`` values
        suffixed by ``_<copy number>``.
    """
    new_attr_dict = feature.attributes.copy()
    copy_num = feature.attributes['extra_copy_number'][0]
    for attr, values in feature.attributes.items():
        if attr.endswith('_id') or attr in ('ID', 'Parent'):
            new_attr_dict[attr] = [f'{values[0]}_{copy_num}']
    return new_attr_dict


def make_gtf_line(attr_dict: Attributes, feature: Feature) -> str:
    """Format a GTF line.

    Parameters
    ----------
    attr_dict : dict of str to list of str
        Attributes to write.
    feature : Feature
        Feature providing the other columns.

    Returns
    -------
    str
        The line without a trailing newline.
    """
    attributes = ''.join(
        f'{key} "{",".join(values)}"; '
        for key, values in attr_dict.items()
        if key != 'copy_id' and len(values) > 0
    )
    columns = (
        feature.seqid,
        feature.source,
        feature.featuretype,
        str(feature.start),
        str(feature.end),
        '.',
        feature.strand,
        '.',
        attributes,
    )
    return '\t'.join(columns)
