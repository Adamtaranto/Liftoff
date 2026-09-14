"""Choose the best chain of aligned blocks for a gene.

Aligned blocks become nodes of a directed acyclic graph. Edges connect blocks
that can follow one another in a consistent gene model (same target sequence
and strand, increasing coordinates, bounded distance, no overlap with other
lifted genes). Node weights penalise mismatches inside exons, and edge weights
penalise exon bases left unaligned between two blocks. The minimum-weight path
from a start sentinel to an end sentinel is the lifted gene model, whose child
feature coordinates are then converted to the target.

The insertion order of nodes and edges is significant: it breaks ties between
equal-weight paths, so the construction below follows earlier versions exactly.
"""

from __future__ import annotations

from bisect import bisect_left, bisect_right
from collections.abc import Mapping
from itertools import pairwise

import networkx as nx

from liftoff.config import LiftoffConfig
from liftoff.models import AlignedSegment, Feature, FeatureHierarchy
from liftoff.utils import (
    FeatureLocations,
    LiftedFeatures,
    count_overlap,
    find_overlaps,
    get_relative_child_coord,
    get_strand,
    merge_children_intervals,
)

_START_ID = -1
_END_ID = -2


def find_best_mapping(
    alignments: list[AlignedSegment],
    query_length: int,
    parent: Feature,
    feature_hierarchy: FeatureHierarchy,
    previous_feature_start: int,
    previous_feature_ref_start: int,
    previous_gene_seq: str,
    feature_locations: FeatureLocations | None,
    lifted_features_list: LiftedFeatures,
    config: LiftoffConfig,
) -> tuple[dict[str, Feature], float, float]:
    """Lift one gene by finding its best alignment chain.

    Parameters
    ----------
    alignments : list of AlignedSegment
        All blocks for one copy of the gene. The list is not modified, but the
        blocks on the chosen path have their start trimmed where they overlap.
    query_length : int
        Length of the extracted gene sequence.
    parent : Feature
        Reference top-level feature.
    feature_hierarchy : FeatureHierarchy
        Reference features.
    previous_feature_start : int
        Target start of the lifted upstream neighbour (0 if none).
    previous_feature_ref_start : int
        Reference start of the upstream neighbour (0 if none).
    previous_gene_seq : str
        Target sequence of the lifted upstream neighbour (``""`` if none).
    feature_locations : IntervalIndex or None
        Already lifted features to avoid; ``None`` disables overlap checks.
    lifted_features_list : dict of str to list of Feature
        Currently lifted features.
    config : LiftoffConfig
        Run configuration (penalties and distance factor).

    Returns
    -------
    mapped_children : dict of str to Feature
        Lifted child features keyed by ID (empty if nothing mapped).
    coverage : float
        Fraction of child feature bases aligned.
    identity : float
        Sequence identity over child features.
    """
    children = feature_hierarchy.children[parent.id]
    children_coords = merge_children_intervals(children)
    node_dict, aln_graph = _initialize_graph()
    head_nodes = _add_single_alignments(
        node_dict,
        aln_graph,
        alignments,
        children_coords,
        parent,
        previous_feature_start,
        previous_feature_ref_start,
        previous_gene_seq,
        feature_locations,
        feature_hierarchy.parents,
        lifted_features_list,
        config,
    )
    for head_node in head_nodes:
        _add_edges(
            head_node,
            node_dict,
            aln_graph,
            parent,
            children_coords,
            feature_locations,
            feature_hierarchy.parents,
            lifted_features_list,
            config,
        )
    _add_target_node(aln_graph, node_dict, query_length, children_coords, parent, config)
    shortest_path_nodes = _find_shortest_path(node_dict, aln_graph)
    if not shortest_path_nodes:
        return {}, 0, 0
    return _convert_all_children_coords(shortest_path_nodes, children, parent)


def _initialize_graph() -> tuple[dict[int, AlignedSegment], nx.DiGraph[int]]:
    """Create the graph with the start sentinel as node 0."""
    aln_graph: nx.DiGraph[int] = nx.DiGraph()
    node_dict = {
        0: AlignedSegment(_START_ID, 'start', 'start', -1, -1, -1, -1, False, []),
    }
    aln_graph.add_node(0)
    return node_dict, aln_graph


def _add_single_alignments(
    node_dict: dict[int, AlignedSegment],
    aln_graph: nx.DiGraph[int],
    alignments: list[AlignedSegment],
    children_coords: list[list[int]],
    parent: Feature,
    previous_feature_start: int,
    previous_feature_ref_start: int,
    previous_gene_seq: str,
    feature_locations: FeatureLocations | None,
    parent_dict: Mapping[str, Feature],
    lifted_features_list: LiftedFeatures,
    config: LiftoffConfig,
) -> list[int]:
    """Add blocks as nodes, linking consecutive blocks of the same alignment.

    Returns the node numbers that begin each alignment's chain (head nodes).
    """
    if feature_locations is not None:
        alignments = _remove_alignments_with_overlap(
            alignments, feature_locations, parent, parent_dict, lifted_features_list
        )
    alignments.sort(key=lambda x: x.query_block_start)
    alignments = _sort_alignments(
        parent, previous_feature_start, previous_feature_ref_start, previous_gene_seq, alignments
    )
    node_num = 1
    previous_aln_id: int | None = None
    head_nodes: list[int] = []
    previous_node = 0
    for aln in alignments:
        if aln.aln_id != previous_aln_id:
            # A new SAM record starts: its first block hangs off the start node.
            previous_node = 0
            previous_aln_id = aln.aln_id
        if _is_valid_alignment(
            previous_node,
            aln,
            node_dict,
            parent,
            feature_locations,
            parent_dict,
            lifted_features_list,
        ):
            if previous_node == 0:
                head_nodes.append(node_num)
            node_dict[node_num] = aln
            aln_graph.add_node(
                node_num, weight=_get_node_weight(aln, children_coords, parent, config)
            )
            edge_weight = _get_edge_weight(
                node_dict[previous_node], aln, children_coords, parent, config
            )
            aln_graph.add_edge(previous_node, node_num, cost=edge_weight)
            previous_node = node_num
            node_num += 1
    return head_nodes


def _sort_alignments(
    parent: Feature,
    previous_feature_start: int,
    previous_feature_ref_start: int,
    previous_gene_seq: str,
    alignments: list[AlignedSegment],
) -> list[AlignedSegment]:
    """Order alignments so the most plausible (by neighbour distance) come first.

    Blocks are ranked by whether they are on the neighbour's sequence and by
    how well their distance from the neighbour matches the reference; each SAM
    record then keeps the rank of its best block, and blocks within a record
    stay in query order.
    """
    expected_distance = parent.start - previous_feature_ref_start
    ranked = sorted(
        alignments,
        key=lambda x: (
            previous_gene_seq != x.reference_name,
            abs(
                expected_distance
                - (x.reference_block_start - x.query_block_start - previous_feature_start)
            ),
        ),
    )
    order_dict: dict[int, int] = {}
    for aln in ranked:
        order_dict.setdefault(aln.aln_id, len(order_dict))
    return sorted(ranked, key=lambda x: (order_dict[x.aln_id], x.query_block_start))


def _remove_alignments_with_overlap(
    alignments: list[AlignedSegment],
    feature_locations: FeatureLocations,
    parent: Feature,
    parent_dict: Mapping[str, Feature],
    lifted_features_list: LiftedFeatures,
) -> list[AlignedSegment]:
    """Drop whole alignments whose span overlaps another lifted feature."""
    # NOTE: iteration over a set of ints mirrors earlier versions; the result
    # is re-sorted by the caller, so only ties can be affected.
    aln_ids = {aln.aln_id for aln in alignments}
    alignments_to_keep: list[AlignedSegment] = []
    for aln_id in aln_ids:
        group = [aln for aln in alignments if aln.aln_id == aln_id]
        min_ref_start = min(aln.reference_block_start for aln in group)
        max_ref_end = max(aln.reference_block_end for aln in group)
        overlaps = find_overlaps(
            min_ref_start,
            max_ref_end,
            group[0].reference_name,
            get_strand(group[0], parent),
            group[0].query_name,
            feature_locations,
            parent_dict,
            lifted_features_list,
            0,
        )
        if not overlaps:
            alignments_to_keep += group
    return alignments_to_keep


def _is_valid_alignment(
    previous_node: int,
    aln: AlignedSegment,
    node_dict: dict[int, AlignedSegment],
    parent: Feature,
    feature_locations: FeatureLocations | None,
    parent_dict: Mapping[str, Feature],
    lifted_features_list: LiftedFeatures,
) -> bool:
    """Return whether a block may be added after ``previous_node``."""
    if aln.reference_block_start == -1:
        return False
    return previous_node == 0 or not _spans_overlap_region(
        node_dict[previous_node], aln, parent, feature_locations, parent_dict, lifted_features_list
    )


def _spans_overlap_region(
    from_node: AlignedSegment,
    to_node: AlignedSegment,
    parent: Feature,
    feature_locations: FeatureLocations | None,
    parent_dict: Mapping[str, Feature],
    lifted_features_list: LiftedFeatures,
) -> bool:
    """Return whether the gap between two blocks contains another lifted feature."""
    if from_node.reference_name != to_node.reference_name or feature_locations is None:
        return False
    node_overlap = _get_node_overlap(from_node, to_node)
    overlaps = find_overlaps(
        from_node.reference_block_end,
        to_node.reference_block_start + node_overlap,
        from_node.reference_name,
        get_strand(from_node, parent),
        from_node.query_name,
        feature_locations,
        parent_dict,
        lifted_features_list,
        0,
    )
    return len(overlaps) > 0


def _get_node_overlap(from_node: AlignedSegment, to_node: AlignedSegment) -> int:
    """Return how many query bases ``to_node`` re-aligns from ``from_node``."""
    if from_node.reference_name == 'start' or to_node.reference_name == 'start':
        return 0
    return max(0, from_node.query_block_end - to_node.query_block_start + 1)


def _count_in_range(sorted_positions: list[int], low: int, high: int) -> int:
    """Count sorted positions within the closed interval ``[low, high]``."""
    return bisect_right(sorted_positions, high) - bisect_left(sorted_positions, low)


def _get_node_weight(
    aln: AlignedSegment, children_coords: list[list[int]], parent: Feature, config: LiftoffConfig
) -> float:
    """Mismatch penalty for bases of ``aln`` that fall inside child features."""
    weight = 0
    for child_start_abs, child_end_abs in children_coords:
        relative_start = get_relative_child_coord(parent, child_start_abs, aln.is_reverse)
        relative_end = get_relative_child_coord(parent, child_end_abs, aln.is_reverse)
        low, high = min(relative_start, relative_end), max(relative_start, relative_end)
        weight += _count_in_range(aln.mismatches, low, high) * config.mismatch
    return weight


def _get_edge_weight(
    from_node: AlignedSegment,
    to_node: AlignedSegment,
    children_coords: list[list[int]],
    parent: Feature,
    config: LiftoffConfig,
) -> float:
    """Affine gap penalty for child-feature bases skipped between two blocks."""
    node_overlap = _get_node_overlap(from_node, to_node)
    unaligned_range = (from_node.query_block_end + 1, to_node.query_block_start + node_overlap - 1)
    is_reverse = to_node.is_reverse if from_node.reference_name == 'start' else from_node.is_reverse
    unaligned_exon_bases = 0
    for child_start_abs, child_end_abs in children_coords:
        relative_start = get_relative_child_coord(parent, child_start_abs, is_reverse)
        relative_end = get_relative_child_coord(parent, child_end_abs, is_reverse)
        child_start, child_end = (
            min(relative_start, relative_end),
            max(relative_start, relative_end),
        )
        overlap = count_overlap(child_start, child_end, min(unaligned_range), max(unaligned_range))
        if (
            overlap == 1
            and unaligned_range[0] == unaligned_range[1] + 1
            and from_node.reference_name == to_node.reference_name
        ):
            # Adjacent in the query: penalise the gap opened in the target.
            unaligned_exon_bases += (
                (to_node.reference_block_start + node_overlap) - from_node.reference_block_end - 1
            ) * config.gap_extend
        else:
            unaligned_exon_bases += max(0, overlap) * config.gap_extend
    if unaligned_exon_bases > 0:
        unaligned_exon_bases += config.gap_open - config.gap_extend
    return unaligned_exon_bases


def _add_edges(
    head_node: int,
    node_dict: dict[int, AlignedSegment],
    aln_graph: nx.DiGraph[int],
    parent: Feature,
    children_coords: list[list[int]],
    feature_locations: FeatureLocations | None,
    parent_dict: Mapping[str, Feature],
    lifted_features_list: LiftedFeatures,
    config: LiftoffConfig,
) -> None:
    """Connect every compatible block to the head of another alignment."""
    for node_name, from_node in node_dict.items():
        if _is_valid_edge(
            node_name,
            head_node,
            parent,
            feature_locations,
            parent_dict,
            lifted_features_list,
            config,
            node_dict,
        ):
            edge_weight = _get_edge_weight(
                from_node, node_dict[head_node], children_coords, parent, config
            )
            aln_graph.add_edge(node_name, head_node, cost=edge_weight)


def _is_valid_edge(
    from_node_name: int,
    to_node_name: int,
    parent: Feature,
    feature_locations: FeatureLocations | None,
    parent_dict: Mapping[str, Feature],
    lifted_features_list: LiftedFeatures,
    config: LiftoffConfig,
    node_dict: dict[int, AlignedSegment],
) -> bool:
    """Return whether ``to_node`` may directly follow ``from_node`` in a chain."""
    from_node, to_node = node_dict[from_node_name], node_dict[to_node_name]
    if from_node.aln_id == _START_ID:
        # The start sentinel is already linked to every head node.
        return False
    if from_node.aln_id == to_node.aln_id:
        return False
    if from_node.query_block_end >= to_node.query_block_end:
        return False
    if from_node.is_reverse != to_node.is_reverse:
        return False
    if from_node.reference_name != to_node.reference_name:
        return False
    expected_distance = to_node.query_block_end - from_node.query_block_start
    actual_distance = to_node.reference_block_end - from_node.reference_block_start
    if to_node.reference_block_start < from_node.reference_block_end:
        return False
    if actual_distance > config.distance_factor * expected_distance:
        return False
    return not _spans_overlap_region(
        from_node, to_node, parent, feature_locations, parent_dict, lifted_features_list
    )


def _add_target_node(
    aln_graph: nx.DiGraph[int],
    node_dict: dict[int, AlignedSegment],
    query_length: int,
    children_coords: list[list[int]],
    parent: Feature,
    config: LiftoffConfig,
) -> None:
    """Add the end sentinel and connect every node without successors to it."""
    end_node = max(node_dict) + 1
    node_dict[end_node] = AlignedSegment(
        _END_ID, 'end', 'end', query_length, query_length, query_length, query_length, False, []
    )
    aln_graph.add_node(end_node)
    for node_name in node_dict:
        if node_name != end_node and not any(True for _ in aln_graph.successors(node_name)):
            edge_weight = _get_edge_weight(
                node_dict[node_name], node_dict[end_node], children_coords, parent, config
            )
            aln_graph.add_edge(node_name, end_node, cost=edge_weight)


def _find_shortest_path(
    node_dict: dict[int, AlignedSegment], aln_graph: nx.DiGraph[int]
) -> list[AlignedSegment]:
    """Return blocks on the minimum-weight start-to-end path (sentinels excluded)."""

    def weight(u: int, v: int, data: dict[str, float]) -> float:
        # Each node's weight is split across its incoming and outgoing edges.
        nodes = aln_graph.nodes
        node_weights: float = nodes[u].get('weight', 0) / 2 + nodes[v].get('weight', 0) / 2
        return node_weights + data.get('cost', 1)

    path = nx.shortest_path(aln_graph, source=0, target=len(node_dict) - 1, weight=weight)
    shortest_path_nodes = [node_dict[node] for node in path[1:-1]]
    _trim_path_boundaries(shortest_path_nodes)
    return shortest_path_nodes


def _trim_path_boundaries(shortest_path_nodes: list[AlignedSegment]) -> None:
    """Remove query bases shared between consecutive blocks from the later block."""
    for from_node, to_node in pairwise(shortest_path_nodes):
        node_overlap = _get_node_overlap(from_node, to_node)
        to_node.query_block_start += node_overlap
        to_node.reference_block_start += node_overlap


def _convert_all_children_coords(
    shortest_path_nodes: list[AlignedSegment], children: list[Feature], parent: Feature
) -> tuple[dict[str, Feature], float, float]:
    """Project child features through the chosen blocks and score the mapping."""
    shortest_path_nodes.sort(key=lambda x: x.query_block_start)
    mapped_children: dict[str, Feature] = {}
    total_bases = mismatches = insertions = deletions = 0
    first_node = shortest_path_nodes[0]
    for child in children:
        # Annotated bounds: a single-level feature is its own child, and its
        # start/end may include the alignment flank.
        child_start, child_end = child.annotated_bounds
        total_bases += child_end - child_start + 1
        nearest_start, nearest_end, relative_start, relative_end = (
            _find_nearest_aligned_start_and_end(child_start, child_end, shortest_path_nodes, parent)
        )
        if nearest_start == -1 or nearest_end == -1:
            deletions += child_end - child_start + 1
            continue
        lifted_start, start_node = _convert_coord(nearest_start, shortest_path_nodes)
        lifted_end, end_node = _convert_coord(nearest_end, shortest_path_nodes)
        deletions += _find_deletions(start_node, end_node, shortest_path_nodes)
        deletions += (nearest_start - relative_start) + (relative_end - nearest_end)
        mismatches += _find_mismatched_bases(child_start, child_end, shortest_path_nodes, parent)
        insertions += _find_insertions(start_node, end_node, shortest_path_nodes)
        if 'ID' not in child.attributes:
            child.attributes['ID'] = [child.id]
        new_child = Feature(
            id=child.id,
            featuretype=child.featuretype,
            seqid=first_node.reference_name,
            source='Liftoff',
            strand=get_strand(first_node, parent),
            start=min(lifted_start, lifted_end) + 1,
            end=max(lifted_start, lifted_end) + 1,
            frame=child.frame,
            attributes=dict(child.attributes),
        )
        mapped_children[new_child.id] = new_child
    alignment_length = total_bases + insertions
    coverage = (total_bases - deletions) / total_bases
    identity = (alignment_length - insertions - mismatches - deletions) / alignment_length
    return mapped_children, coverage, identity


def _relative_bounds(start: int, end: int, parent: Feature, is_reverse: bool) -> tuple[int, int]:
    """Return a child's ``(start, end)`` offsets within the aligned parent sequence."""
    coord1 = get_relative_child_coord(parent, start, is_reverse)
    coord2 = get_relative_child_coord(parent, end, is_reverse)
    return min(coord1, coord2), max(coord1, coord2)


def _find_nearest_aligned_start_and_end(
    child_start: int, child_end: int, shortest_path_nodes: list[AlignedSegment], parent: Feature
) -> tuple[int, int, int, int]:
    """Find the aligned query positions closest to a child's ends (-1 if none)."""
    relative_start, relative_end = _relative_bounds(
        child_start, child_end, parent, shortest_path_nodes[0].is_reverse
    )
    nearest_start = _find_nearest_aligned_start(relative_start, relative_end, shortest_path_nodes)
    nearest_end = _find_nearest_aligned_end(shortest_path_nodes, relative_end, relative_start)
    return nearest_start, nearest_end, relative_start, relative_end


def _find_nearest_aligned_start(
    relative_start: int, relative_end: int, shortest_path_nodes: list[AlignedSegment]
) -> int:
    """First aligned query position at or after the child's start, within the child."""
    for node in shortest_path_nodes:
        if relative_start <= node.query_block_end:
            if relative_start >= node.query_block_start:
                return relative_start
            if node.query_block_start < relative_end:
                return node.query_block_start
            return -1
    return -1


def _find_nearest_aligned_end(
    shortest_path_nodes: list[AlignedSegment], relative_end: int, relative_start: int
) -> int:
    """Last aligned query position at or before the child's end, within the child."""
    nearest_end = -1
    node = shortest_path_nodes[-1]
    for i, node in enumerate(shortest_path_nodes):
        if relative_end <= node.query_block_end:
            if relative_end >= node.query_block_start:
                nearest_end = relative_end
            elif i > 0 and shortest_path_nodes[i - 1].query_block_end > relative_start:
                nearest_end = shortest_path_nodes[i - 1].query_block_end
            break
    if nearest_end == -1 and relative_start < node.query_block_end < relative_end:
        nearest_end = node.query_block_end
    return nearest_end


def _convert_coord(
    relative_coord: int, shortest_path_nodes: list[AlignedSegment]
) -> tuple[int, int]:
    """Map a query offset to a target coordinate via the last block containing it."""
    lifted_coord = 0
    node_index = 0
    for i, node in enumerate(shortest_path_nodes):
        if node.query_block_start <= relative_coord <= node.query_block_end:
            lifted_coord = node.reference_block_start + (relative_coord - node.query_block_start)
            node_index = i
    return lifted_coord, node_index


def _find_deletions(start: int, end: int, shortest_path_nodes: list[AlignedSegment]) -> int:
    """Query bases skipped between blocks ``start`` and ``end``."""
    return sum(
        shortest_path_nodes[i].query_block_start - shortest_path_nodes[i - 1].query_block_end - 1
        for i in range(start + 1, end + 1)
    )


def _find_insertions(start: int, end: int, shortest_path_nodes: list[AlignedSegment]) -> int:
    """Target bases skipped between blocks ``start`` and ``end``."""
    return sum(
        shortest_path_nodes[i].reference_block_start
        - shortest_path_nodes[i - 1].reference_block_end
        - 1
        for i in range(start + 1, end + 1)
    )


def _find_mismatched_bases(
    start: int, end: int, shortest_path_nodes: list[AlignedSegment], parent: Feature
) -> int:
    """Count mismatches of all path blocks within a child feature."""
    relative_start, relative_end = _relative_bounds(
        start, end, parent, shortest_path_nodes[0].is_reverse
    )
    return sum(
        _count_in_range(node.mismatches, relative_start, relative_end)
        for node in shortest_path_nodes
    )
