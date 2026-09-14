"""Coordinate, identifier and overlap helpers shared by pipeline stages."""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field

from liftoff.intervals import IntervalIndex
from liftoff.models import AlignedSegment, Feature

#: Maps a copy-aware feature ID (e.g. ``"gene1_0"``) to its lifted features,
#: with the top-level feature first.
LiftedFeatures = dict[str, list[Feature]]

#: Interval index of lifted top-level features keyed by copy ID.
FeatureLocations = IntervalIndex[Feature]


def count_overlap(start1: int, end1: int, start2: int, end2: int) -> int:
    """Return the number of positions shared by two closed intervals.

    Parameters
    ----------
    start1, end1 : int
        Inclusive bounds of the first interval.
    start2, end2 : int
        Inclusive bounds of the second interval.

    Returns
    -------
    int
        Overlap length; zero or negative when the intervals are disjoint.
    """
    return min(end1, end2) - max(start1, start2) + 1


def get_relative_child_coord(parent: Feature, coord: int, is_reverse: bool) -> int:
    """Convert a genomic coordinate to an offset within the parent sequence.

    Parameters
    ----------
    parent : Feature
        Top-level feature whose sequence was aligned.
    coord : int
        One-based genomic coordinate inside ``parent``.
    is_reverse : bool
        If true the offset is measured from the parent's end, matching the
        orientation of a reverse-strand alignment.

    Returns
    -------
    int
        Zero-based offset from the start (or end) of the parent.
    """
    if is_reverse:
        return parent.end - coord
    return coord - parent.start


def merge_children_intervals(children: Iterable[Feature]) -> list[list[int]]:
    """Merge overlapping child feature intervals.

    Annotated coordinates are used, so a single-level feature that is its own
    child does not include the alignment flank.

    Parameters
    ----------
    children : iterable of Feature
        Features to merge.

    Returns
    -------
    list of list of int
        Sorted, non-overlapping ``[start, end]`` intervals.
    """
    return merge_intervals(child.annotated_bounds for child in children)


def merge_intervals(intervals: Iterable[Sequence[int]]) -> list[list[int]]:
    """Merge overlapping closed intervals.

    Parameters
    ----------
    intervals : iterable of sequence of int
        ``(start, end)`` pairs.

    Returns
    -------
    list of list of int
        Sorted, non-overlapping ``[start, end]`` intervals.

    Examples
    --------
    >>> merge_intervals([(5, 9), (1, 3), (2, 4)])
    [[1, 4], [5, 9]]
    """
    ordered: list[list[int]] = sorted(
        ([start, end] for start, end in intervals), key=lambda iv: iv[0]
    )
    if not ordered:
        return []
    merged = [ordered[0]]
    for current in ordered:
        previous = merged[-1]
        if current[0] <= previous[1]:
            previous[1] = max(previous[1], current[1])
        else:
            merged.append(current)
    return merged


def get_parent_list(feature_list: Mapping[str, Sequence[Feature]]) -> list[Feature]:
    """Return the top-level feature of every non-empty lifted feature group.

    Parameters
    ----------
    feature_list : mapping of str to sequence of Feature
        Lifted feature groups; the first element of each is the top-level
        feature.

    Returns
    -------
    list of Feature
        Top-level features in mapping order.
    """
    return [features[0] for features in feature_list.values() if len(features) > 0]


def clear_scores(feature_list: LiftedFeatures, parent_dict: Mapping[str, Feature]) -> None:
    """Reset scores of lifted top-level features before another lift-over round.

    Parameters
    ----------
    feature_list : dict of str to list of Feature
        Lifted features to update in place.
    parent_dict : mapping of str to Feature
        Reference top-level features; only features whose ID is a key are
        reset.
    """
    for features in feature_list.values():
        for feature in features:
            if feature.id in parent_dict:
                feature.score = -1


@dataclass(slots=True)
class ParentOrder:
    """Top-level features sorted by genomic position.

    Replaces the object-dtype NumPy array used previously; lookups by ID are
    now constant time instead of a full scan.

    Attributes
    ----------
    features : list of Feature
        Features sorted by ``(seqid, start)``.
    """

    features: list[Feature]
    ids: list[str] = field(init=False)
    _first_index: dict[str, int] = field(init=False)

    def __post_init__(self) -> None:
        self.ids = [feature.id for feature in self.features]
        self._first_index = {}
        for index, feature_id in enumerate(self.ids):
            # Duplicate IDs (extra copies) resolve to their first occurrence.
            self._first_index.setdefault(feature_id, index)

    @classmethod
    def from_features(cls, parents: list[Feature]) -> ParentOrder:
        """Sort features by position and build an index.

        Parameters
        ----------
        parents : list of Feature
            Features to order. The list is sorted in place.

        Returns
        -------
        ParentOrder
            The ordered features.
        """
        parents.sort(key=lambda x: (x.seqid, x.start))
        return cls(parents)

    def index(self, feature_id: str) -> int | None:
        """Return the position of the first feature with ``feature_id``.

        Parameters
        ----------
        feature_id : str
            Feature ID to look up.

        Returns
        -------
        int or None
            Zero-based position, or ``None`` if the ID is absent.
        """
        return self._first_index.get(feature_id)

    def __len__(self) -> int:
        return len(self.features)


def find_parent_order(parents: list[Feature]) -> ParentOrder:
    """Sort top-level features by chromosome and start coordinate.

    Parameters
    ----------
    parents : list of Feature
        Features to order; sorted in place.

    Returns
    -------
    ParentOrder
        Position-ordered features with an ID index.
    """
    return ParentOrder.from_features(parents)


def convert_id_to_original(feature_id: str) -> str:
    """Strip the copy suffix (and any ``_frag`` suffix) from a lifted feature ID.

    Parameters
    ----------
    feature_id : str
        Copy-aware ID such as ``"gene1_0"`` or ``"gene1_2_frag1"``.

    Returns
    -------
    str
        The reference feature ID (``"gene1"``).
    """
    frag_split = feature_id.split('_frag')[0]
    copy_tag_len = len(frag_split.split('_')[-1])
    return frag_split[: -copy_tag_len - 1]


def get_copy_tag(feature_id: str) -> str:
    """Return the copy suffix of a lifted feature ID, including the underscore.

    Parameters
    ----------
    feature_id : str
        Copy-aware ID such as ``"gene1_0"``.

    Returns
    -------
    str
        The suffix, e.g. ``"_0"``.
    """
    copy_tag_len = len(feature_id.split('_')[-1])
    return feature_id[-copy_tag_len - 1 :]


def get_strand(aln: AlignedSegment, parent: Feature) -> str:
    """Return the strand of a lifted feature given its alignment orientation.

    Parameters
    ----------
    aln : AlignedSegment
        Alignment block of the feature.
    parent : Feature
        Reference top-level feature.

    Returns
    -------
    str
        ``"+"`` or ``"-"`` (or the parent's strand when unchanged).
    """
    if aln.is_reverse:
        return '+' if parent.strand == '-' else '-'
    return parent.strand


def find_overlaps(
    start: int,
    end: int,
    chrm: str,
    strand: str,
    feature_name: str,
    intervals: FeatureLocations,
    parent_dict: Mapping[str, Feature],
    lifted_features_list: Mapping[str, object],
    max_overlap_non_copies: float,
) -> list[tuple[int, int, str, Feature]]:
    """Find lifted features that overlap a candidate location too much.

    Overlaps are ignored when the two features also overlap in the reference
    annotation, when the other feature has since been removed from the lifted
    set, or when they share less than the allowed fraction of the shorter
    feature. Extra copies may not overlap anything.

    Parameters
    ----------
    start, end : int
        Zero-based inclusive candidate location.
    chrm : str
        Target sequence name.
    strand : str
        Candidate strand.
    feature_name : str
        Copy-aware ID of the candidate feature.
    intervals : IntervalIndex of Feature
        Locations of already lifted features.
    parent_dict : mapping of str to Feature
        Reference top-level features.
    lifted_features_list : mapping
        Currently lifted feature groups (only keys are used).
    max_overlap_non_copies : float
        Allowed overlap fraction when neither feature is an extra copy.

    Returns
    -------
    list of tuple
        ``(start, end, copy_id, feature)`` for each disallowed overlap.
    """
    incorrect_overlaps: list[tuple[int, int, str, Feature]] = []
    is_copy = get_copy_tag(feature_name) != '_0'
    for overlap in intervals.find(start, end):
        other = overlap.value
        if other.seqid != chrm or other.strand != strand:
            continue
        shortest_feature_length = min(end - start, overlap.end - overlap.start) + 1
        overlap_amount = count_overlap(start, end, overlap.start, overlap.end)
        ref_feature = parent_dict[convert_id_to_original(feature_name)]
        ref_overlap_feature = parent_dict[convert_id_to_original(overlap.key)]
        max_overlap = (
            0.0 if is_copy or get_copy_tag(overlap.key) != '_0' else max_overlap_non_copies
        )
        if (
            not overlaps_in_ref_annotation(ref_feature, ref_overlap_feature)
            and overlap.key in lifted_features_list
            and overlap_amount / shortest_feature_length > max_overlap
        ):
            incorrect_overlaps.append((overlap.start, overlap.end, overlap.key, other))
    return incorrect_overlaps


def overlaps_in_ref_annotation(ref_feature1: Feature, ref_feature2: Feature) -> bool:
    """Return whether two distinct reference features overlap on the same strand.

    Parameters
    ----------
    ref_feature1, ref_feature2 : Feature
        Reference features to compare.

    Returns
    -------
    bool
        ``True`` if the features differ and overlap in the reference.
    """
    if ref_feature1.seqid != ref_feature2.seqid:
        return False
    if ref_feature1.strand != ref_feature2.strand:
        return False
    if ref_feature1.id == ref_feature2.id:
        return False
    return (
        count_overlap(ref_feature1.start, ref_feature1.end, ref_feature2.start, ref_feature2.end)
        > 0
    )


def find_nonoverlapping_upstream_neighbor(
    parent_order: ParentOrder, feature_name: str
) -> str | None:
    """Return the ID of the feature immediately upstream on the same sequence.

    Parameters
    ----------
    parent_order : ParentOrder
        Reference features in genomic order.
    feature_name : str
        ID of the feature whose neighbour is wanted.

    Returns
    -------
    str or None
        Neighbour ID, or ``None`` if the feature is first on its sequence.
    """
    feature_location = parent_order.index(feature_name)
    if feature_location is None:
        raise KeyError(feature_name)
    neighbor_index = feature_location - 1
    if neighbor_index < 0:
        return None
    neighbor = parent_order.features[neighbor_index]
    feature = parent_order.features[feature_location]
    if neighbor.seqid != feature.seqid:
        return None
    return parent_order.ids[neighbor_index]
