"""Type aliases for factory fixtures defined in ``tests/conftest.py``."""

from __future__ import annotations

from collections.abc import Callable

from liftoff.models import Feature

#: ``make_feature(feature_id, featuretype, start, end, *, strand, seqid, frame, parent, ...)``
FeatureFactory = Callable[..., Feature]

#: ``make_sam_line(qname, flag, rname, pos, cigar, seq)``
SamLineFactory = Callable[..., str]
