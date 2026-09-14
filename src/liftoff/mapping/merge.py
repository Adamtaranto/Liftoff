"""Rebuild the full feature tree of a lifted gene from its lifted children."""

from __future__ import annotations

from collections.abc import Mapping

from liftoff.models import Feature, FeatureHierarchy


def merge_lifted_features(
    mapped_children: Mapping[str, Feature],
    parent: Feature,
    unmapped_features: list[Feature],
    aln_cov_threshold: float,
    copy_id: str,
    feature_order: Mapping[str, int],
    feature_hierarchy: FeatureHierarchy,
    aln_cov: float,
    seq_id: float,
    seq_id_threshold: float,
) -> list[Feature]:
    """Create intermediate and top-level features spanning lifted children.

    Parameters
    ----------
    mapped_children : mapping of str to Feature
        Lifted lowest-level features.
    parent : Feature
        Reference top-level feature.
    unmapped_features : list of Feature
        Receives ``parent`` if the mapping fails the thresholds.
    aln_cov_threshold : float
        Minimum coverage.
    copy_id : str
        Copy-aware ID of this lifted copy.
    feature_order : mapping of str to int
        Sort rank of each feature type.
    feature_hierarchy : FeatureHierarchy
        Reference features.
    aln_cov : float
        Coverage of this mapping.
    seq_id : float
        Sequence identity of this mapping.
    seq_id_threshold : float
        Minimum identity.

    Returns
    -------
    list of Feature
        Lifted features with the top-level feature first, or an empty list if
        the mapping was rejected.
    """
    if len(mapped_children) == 0:
        unmapped_features.append(parent)
        return []
    feature_list: dict[str, Feature] = {}
    non_parents: list[tuple[Feature, str]] = []
    top_target_feature: Feature | None = None
    for child_feature in mapped_children.values():
        feature_list[child_feature.id] = child_feature
        if child_feature.id == parent.id:
            top_target_feature = child_feature
        else:
            non_parents.append((child_feature, child_feature.attributes['Parent'][0]))
    # Walk up the hierarchy one level at a time until the top-level feature
    # has been created.
    while non_parents:
        non_parents, top_target_feature = _create_parents(
            non_parents, parent, feature_hierarchy, feature_list
        )
    return _process_final_features_list(
        feature_list,
        top_target_feature,
        seq_id,
        seq_id_threshold,
        aln_cov,
        aln_cov_threshold,
        unmapped_features,
        parent,
        feature_order,
        copy_id,
    )


def _create_parents(
    non_parents: list[tuple[Feature, str]],
    top_ref_parent: Feature,
    feature_hierarchy: FeatureHierarchy,
    feature_list: dict[str, Feature],
) -> tuple[list[tuple[Feature, str]], Feature | None]:
    """Create the next level of parents for ``non_parents``."""
    added_parent_ids: set[str] = set()
    new_non_parents: list[tuple[Feature, str]] = []
    top_target_feature: Feature | None = None
    for _child, parent_id in non_parents:
        if parent_id in added_parent_ids:
            continue
        target_parent = _make_new_parent(feature_list, parent_id, feature_hierarchy)
        if target_parent.id == top_ref_parent.id:
            top_target_feature = target_parent
        else:
            new_non_parents.append((target_parent, target_parent.attributes['Parent'][0]))
        added_parent_ids.add(parent_id)
    return new_non_parents, top_target_feature


def _make_new_parent(
    feature_list: dict[str, Feature], parent_id: str, feature_hierarchy: FeatureHierarchy
) -> Feature:
    """Create a lifted parent spanning all lifted features that name it as parent."""
    children = [
        feature
        for feature in feature_list.values()
        if 'Parent' in feature.attributes and feature.attributes['Parent'][0] == parent_id
    ]
    if parent_id in feature_hierarchy.parents:
        ref_parent = feature_hierarchy.parents[parent_id]
    else:
        ref_parent = feature_hierarchy.intermediates[parent_id]
    target_parent = Feature(
        id=ref_parent.id,
        featuretype=ref_parent.featuretype,
        seqid=children[0].seqid,
        source='Liftoff',
        strand=children[0].strand,
        start=min(child.start for child in children),
        end=max(child.end for child in children),
        frame=ref_parent.frame,
        attributes=dict(ref_parent.attributes),
    )
    feature_list[target_parent.id] = target_parent
    return target_parent


def _process_final_features_list(
    feature_list: dict[str, Feature],
    top_target_feature: Feature | None,
    seq_id: float,
    seq_id_threshold: float,
    aln_cov: float,
    aln_cov_threshold: float,
    unmapped_features: list[Feature],
    parent: Feature,
    feature_order: Mapping[str, int],
    copy_id: str,
) -> list[Feature]:
    """Order features, apply thresholds and annotate the top-level feature."""
    final_features = [
        feature for feature in feature_list.values() if feature is not top_target_feature
    ]
    # Two stable sorts: by position, then by feature type rank.
    final_features.sort(key=lambda x: (x.seqid, x.start))
    final_features.sort(key=lambda x: feature_order[x.featuretype])
    if aln_cov < aln_cov_threshold or seq_id < seq_id_threshold:
        unmapped_features.append(parent)
        return []
    if top_target_feature is None:
        raise ValueError(f'lifted features of {parent.id} do not include the top-level feature')
    final_features.insert(0, top_target_feature)
    top_target_feature.score = 1 - seq_id
    top_target_feature.attributes['copy_num_ID'] = [copy_id]
    # Values are truncated (not rounded) to five characters.
    top_target_feature.attributes['coverage'] = [str(aln_cov)[0:5]]
    top_target_feature.attributes['sequence_ID'] = [str(seq_id)[0:5]]
    return final_features
