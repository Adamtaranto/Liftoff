"""Resolve lifted features that incorrectly overlap each other.

When two lifted genes overlap in the target but not in the reference, one of
them is probably mis-mapped (typically a member of a gene family mapped onto a
paralogue). The less likely mapping is removed and re-lifted while avoiding
already occupied loci, repeating until no disallowed overlaps remain.
"""

from __future__ import annotations

from collections.abc import Mapping

import gffutils

from liftoff.align.blocks import AlignedSegments
from liftoff.config import LiftoffConfig
from liftoff.intervals import Interval, IntervalIndex
from liftoff.log import get_logger
from liftoff.mapping.lift import lift_all_features
from liftoff.models import Feature, FeatureHierarchy
from liftoff.utils import (
    FeatureLocations,
    LiftedFeatures,
    ParentOrder,
    convert_id_to_original,
    find_nonoverlapping_upstream_neighbor,
    find_overlaps,
    find_parent_order,
    get_parent_list,
)

logger = get_logger(__name__)


def fix_incorrectly_overlapping_features(
    all_lifted_features: LiftedFeatures,
    features_to_check: LiftedFeatures,
    all_aligned_segs: AlignedSegments,
    unmapped_features: list[Feature],
    threshold: float,
    feature_hierarchy: FeatureHierarchy,
    feature_db: gffutils.FeatureDB,
    ref_parent_order: ParentOrder,
    seq_id_threshold: float,
    config: LiftoffConfig,
    max_overlap: float,
) -> None:
    """Detect overlapping lifted features and re-map the likely mis-mapped ones.

    Parameters
    ----------
    all_lifted_features : dict of str to list of Feature
        All lifted features (updated in place).
    features_to_check : dict of str to list of Feature
        Features to test for overlaps.
    all_aligned_segs : dict of str to list of AlignedSegment
        Alignments of every feature, used for re-mapping.
    unmapped_features : list of Feature
        Receives features that cannot be placed without overlaps.
    threshold : float
        Minimum coverage when re-mapping.
    feature_hierarchy : FeatureHierarchy
        Reference features.
    feature_db : gffutils.FeatureDB
        Reference feature database.
    ref_parent_order : ParentOrder
        Reference top-level features in genomic order.
    seq_id_threshold : float
        Minimum identity when re-mapping.
    config : LiftoffConfig
        Run configuration.
    max_overlap : float
        Maximum allowed overlap fraction between non-copy features.
    """
    features_to_remap, feature_locations = _check_homologues(
        all_lifted_features,
        features_to_check,
        feature_hierarchy.parents,
        ref_parent_order,
        max_overlap,
    )
    _resolve_overlapping_homologues(
        all_aligned_segs,
        all_lifted_features,
        features_to_remap,
        unmapped_features,
        threshold,
        feature_hierarchy,
        feature_db,
        ref_parent_order,
        seq_id_threshold,
        feature_locations,
        config,
        max_overlap,
    )


def _check_homologues(
    all_lifted_features: LiftedFeatures,
    lifted_features_to_check: LiftedFeatures,
    parent_dict: Mapping[str, Feature],
    ref_parent_order: ParentOrder,
    max_overlap: float,
) -> tuple[set[str], FeatureLocations]:
    """Return copy IDs that should be re-mapped and the current feature index."""
    all_feature_list = get_parent_list(all_lifted_features)
    features_to_check_list = get_parent_list(lifted_features_to_check)
    target_parent_order = find_parent_order(all_feature_list)
    remap_features: set[str] = set()
    feature_locations = build_interval_list(all_feature_list)
    for feature in features_to_check_list:
        copy_id = feature.attributes['copy_num_ID'][0]
        overlaps = find_overlaps(
            feature.start - 1,
            feature.end - 1,
            feature.seqid,
            feature.strand,
            copy_id,
            feature_locations,
            parent_dict,
            all_lifted_features,
            max_overlap,
        )
        for _start, _end, overlap_copy_id, overlap_feature in overlaps:
            if overlap_copy_id != copy_id:
                feature_to_remap = _find_feature_to_remap(
                    feature, overlap_feature, ref_parent_order, target_parent_order, remap_features
                )
                remap_features.add(feature_to_remap.attributes['copy_num_ID'][0])
    return remap_features, feature_locations


def build_interval_list(features: list[Feature]) -> FeatureLocations:
    """Index lifted top-level features by zero-based target coordinates.

    Parameters
    ----------
    features : list of Feature
        Lifted top-level features carrying a ``copy_num_ID`` attribute.

    Returns
    -------
    IntervalIndex of Feature
        Index keyed by copy ID.
    """
    return IntervalIndex(
        Interval(feature.start - 1, feature.end - 1, feature.attributes['copy_num_ID'][0], feature)
        for feature in features
    )


def _find_feature_to_remap(
    feature: Feature,
    overlap_feature: Feature,
    ref_parent_order: ParentOrder,
    target_parent_order: ParentOrder,
    remap_features: set[str],
) -> Feature:
    """Decide which of two overlapping features is more likely mis-mapped.

    Criteria, in order: extra copies yield to originals; a feature already
    scheduled for re-mapping is chosen again; lower sequence identity loses;
    the feature further from its reference neighbour loses; finally the
    shorter feature loses.
    """
    feature_is_copy = _is_copy(feature)
    overlap_feature_is_copy = _is_copy(overlap_feature)
    if feature_is_copy and not overlap_feature_is_copy:
        return feature
    if overlap_feature_is_copy and not feature_is_copy:
        return overlap_feature
    if feature.attributes['copy_num_ID'][0] in remap_features:
        return feature
    if overlap_feature.attributes['copy_num_ID'][0] in remap_features:
        return overlap_feature
    # Scores are 1 - identity, so a lower score means higher identity.
    if feature.score < overlap_feature.score:
        return overlap_feature
    if overlap_feature.score < feature.score:
        return feature
    out_of_order = _find_out_of_order_feature(
        feature_is_copy,
        ref_parent_order,
        target_parent_order,
        feature,
        overlap_feature_is_copy,
        overlap_feature,
    )
    if out_of_order is not None:
        return out_of_order
    if _is_shorter(feature, overlap_feature):
        return feature
    return overlap_feature


def _is_copy(feature: Feature) -> bool:
    """Return whether a lifted feature is an extra copy (copy tag other than ``_0``)."""
    return feature.attributes['copy_num_ID'][0][-2:] != '_0'


def _find_out_of_order_feature(
    feature_is_copy: bool,
    ref_parent_order: ParentOrder,
    target_parent_order: ParentOrder,
    feature: Feature,
    overlap_feature_is_copy: bool,
    overlap_feature: Feature,
) -> Feature | None:
    """Return the feature whose target placement disagrees more with the reference order."""
    distance_feature = (
        None
        if feature_is_copy
        else _check_order(ref_parent_order, target_parent_order, feature, overlap_feature)
    )
    distance_overlap = (
        None
        if overlap_feature_is_copy
        else _check_order(ref_parent_order, target_parent_order, overlap_feature, feature)
    )
    if distance_feature is not None and distance_overlap is not None:
        return overlap_feature if distance_feature < distance_overlap else feature
    return None


def _check_order(
    ref_parent_order: ParentOrder,
    target_parent_order: ParentOrder,
    feature: Feature,
    overlap_feature: Feature,
) -> int | None:
    """Distance in the target between a feature and its reference upstream neighbour."""
    ref_neighbor_upstream = find_nonoverlapping_upstream_neighbor(ref_parent_order, feature.id)
    if ref_neighbor_upstream is None:
        return None
    feature_location = target_parent_order.index(feature.id)
    neighbor_location = target_parent_order.index(ref_neighbor_upstream)
    if feature_location is None or neighbor_location is None:
        return None
    low, high = sorted((feature_location, neighbor_location))
    distance = high - low
    # The overlapping feature itself does not count as an intervening gene.
    if overlap_feature.id in target_parent_order.ids[low : high + 1]:
        distance -= 1
    return distance


def _is_shorter(feature1: Feature, feature2: Feature) -> bool:
    """Return whether ``feature1`` spans fewer bases than ``feature2``."""
    return (feature1.end - feature1.start) < (feature2.end - feature2.start)


def _resolve_overlapping_homologues(
    all_aligned_segs: AlignedSegments,
    lifted_feature_list: LiftedFeatures,
    features_to_remap: set[str],
    unmapped_features: list[Feature],
    threshold: float,
    feature_hierarchy: FeatureHierarchy,
    feature_db: gffutils.FeatureDB,
    ref_parent_order: ParentOrder,
    seq_id_threshold: float,
    feature_locations: FeatureLocations,
    config: LiftoffConfig,
    max_overlap: float,
) -> None:
    """Re-lift overlapping features repeatedly, avoiding occupied loci."""
    max_iterations = 10 * len(features_to_remap)
    iteration = 0
    while features_to_remap:
        iteration += 1
        if iteration > max_iterations:
            break
        logger.debug('re-mapping %d overlapping features', len(features_to_remap))
        # Sorted for deterministic results across runs (set order of strings
        # depends on hash randomisation).
        ordered = sorted(features_to_remap)
        aligned_segs_for_remap: AlignedSegments = {}
        for copy_id in ordered:
            del lifted_feature_list[copy_id]
            aligned_segs_for_remap[copy_id] = all_aligned_segs[copy_id]
        lift_all_features(
            aligned_segs_for_remap,
            threshold,
            feature_db,
            feature_hierarchy,
            unmapped_features,
            lifted_feature_list,
            seq_id_threshold,
            feature_locations,
            config,
            ref_parent_order,
        )
        remapped = {
            copy_id: lifted_feature_list[copy_id]
            for copy_id in ordered
            if copy_id in lifted_feature_list
        }
        features_to_remap, feature_locations = _check_homologues(
            lifted_feature_list, remapped, feature_hierarchy.parents, ref_parent_order, max_overlap
        )
    for copy_id in sorted(features_to_remap):
        unmapped_features.append(feature_hierarchy.parents[convert_id_to_original(copy_id)])
        del lifted_feature_list[copy_id]
