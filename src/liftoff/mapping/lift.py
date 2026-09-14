"""Lift every aligned gene: choose its best mapping and rebuild its features."""

from __future__ import annotations

from collections.abc import Mapping

import gffutils

from liftoff.align.blocks import AlignedSegments
from liftoff.config import LiftoffConfig
from liftoff.io.gff_db import get_feature_order
from liftoff.mapping.best_mapping import find_best_mapping
from liftoff.mapping.merge import merge_lifted_features
from liftoff.models import AlignedSegment, Feature, FeatureHierarchy
from liftoff.utils import (
    FeatureLocations,
    LiftedFeatures,
    ParentOrder,
    convert_id_to_original,
    find_nonoverlapping_upstream_neighbor,
)


def lift_all_features(
    alns: AlignedSegments,
    threshold: float,
    feature_db: gffutils.FeatureDB,
    feature_hierarchy: FeatureHierarchy,
    unmapped_features: list[Feature],
    lifted_feature_list: LiftedFeatures,
    seq_id_threshold: float,
    feature_locations: FeatureLocations | None,
    config: LiftoffConfig,
    ref_parent_order: ParentOrder,
) -> None:
    """Lift all aligned features in reference genome order.

    Processing in genome order means a gene's upstream neighbour has usually
    been lifted already, which helps choose among alternative alignments.

    Parameters
    ----------
    alns : dict of str to list of AlignedSegment
        Aligned blocks keyed by copy-aware feature name.
    threshold : float
        Minimum coverage.
    feature_db : gffutils.FeatureDB
        Reference feature database.
    feature_hierarchy : FeatureHierarchy
        Reference features.
    unmapped_features : list of Feature
        Receives features that could not be lifted.
    lifted_feature_list : dict of str to list of Feature
        Receives successfully lifted features (updated in place).
    seq_id_threshold : float
        Minimum identity.
    feature_locations : IntervalIndex or None
        Lifted features to avoid overlapping, or ``None``.
    config : LiftoffConfig
        Run configuration.
    ref_parent_order : ParentOrder
        Reference top-level features in genomic order.
    """
    features_to_lift = feature_hierarchy.parents
    feature_order = get_feature_order(feature_db)
    for alignment in _sort_alignments(features_to_lift, alns):
        previous_start, previous_seq, previous_ref_start = _find_neighbor_location(
            features_to_lift, alignment, lifted_feature_list, ref_parent_order
        )
        lifted_features, parent_name = _lift_single_feature(
            threshold,
            feature_order,
            features_to_lift,
            feature_hierarchy,
            previous_start,
            previous_ref_start,
            previous_seq,
            unmapped_features,
            alignment,
            seq_id_threshold,
            feature_locations,
            lifted_feature_list,
            config,
        )
        if lifted_features:
            lifted_feature_list[parent_name] = lifted_features


def _sort_alignments(
    parent_dict: Mapping[str, Feature], alignments: AlignedSegments
) -> list[list[AlignedSegment]]:
    """Order alignment groups by the reference position of their feature."""
    parent_list = [
        parent_dict[convert_id_to_original(blocks[0].query_name)] for blocks in alignments.values()
    ]
    parent_list.sort(key=lambda x: (x.seqid, x.start))
    # Copies share a parent; the last occurrence determines their rank.
    order_dict = {parent.id: order for order, parent in enumerate(parent_list)}
    values = list(alignments.values())
    values.sort(key=lambda x: order_dict[convert_id_to_original(x[0].query_name)])
    return values


def _find_neighbor_location(
    ref_parents: Mapping[str, Feature],
    alignment: list[AlignedSegment],
    lifted_feature_list: LiftedFeatures,
    ref_parent_order: ParentOrder,
) -> tuple[int, str, int]:
    """Return target start, target sequence and reference start of the upstream neighbour."""
    ref_feature = ref_parents[convert_id_to_original(alignment[0].query_name)]
    ref_neighbor_name = find_nonoverlapping_upstream_neighbor(ref_parent_order, ref_feature.id)
    if ref_neighbor_name is not None:
        ref_neighbor_key = ref_neighbor_name + '_0'
        if ref_neighbor_key in lifted_feature_list:
            lifted_neighbor = lifted_feature_list[ref_neighbor_key][0]
            return (
                lifted_neighbor.start,
                lifted_neighbor.seqid,
                ref_parents[ref_neighbor_name].start,
            )
    return 0, '', 0


def _lift_single_feature(
    threshold: float,
    feature_order: Mapping[str, int],
    features_to_lift: Mapping[str, Feature],
    feature_hierarchy: FeatureHierarchy,
    previous_feature_start: int,
    previous_feature_ref_start: int,
    previous_gene_seq: str,
    unmapped_features: list[Feature],
    aligned_feature: list[AlignedSegment],
    seq_id_threshold: float,
    feature_locations: FeatureLocations | None,
    lifted_features_list: LiftedFeatures,
    config: LiftoffConfig,
) -> tuple[list[Feature], str]:
    """Lift one copy of a feature from its aligned blocks."""
    new_parent_name = aligned_feature[0].query_name
    parent = features_to_lift[convert_id_to_original(new_parent_name)]
    lifted_children, alignment_coverage, seq_id = find_best_mapping(
        aligned_feature,
        parent.end - parent.start + 1,
        parent,
        feature_hierarchy,
        previous_feature_start,
        previous_feature_ref_start,
        previous_gene_seq,
        feature_locations,
        lifted_features_list,
        config,
    )
    lifted_features = merge_lifted_features(
        lifted_children,
        parent,
        unmapped_features,
        threshold,
        new_parent_name,
        feature_order,
        feature_hierarchy,
        alignment_coverage,
        seq_id,
        seq_id_threshold,
    )
    return lifted_features, new_parent_name
