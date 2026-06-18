"""IRfinder-mdl: junction-anchored intron retention from spliced BAM alignments.

Given a reference annotation (GTF) and a sorted+indexed BAM of spliced read
alignments, for every annotated intron the tool counts, at each boundary
(5'/donor and 3'/acceptor), reads that:

  - SPLICE the boundary (a CIGAR `N` op anchored at the boundary within `jitter`
    bp on the reference, with at least `anchor` matched bases on the exonic
    side); and
  - RETAIN the boundary (at least `anchor` matched bases on both the exonic and
    intronic side of the boundary, and no `N` op crossing the boundary window).

Reads that lie entirely inside an intron contribute nothing — only boundary-
crossing reads are evidence.  This mirrors the IRFinder model but is implemented
from scratch in Python for long-read spliced BAMs.
"""

from .version import __version__

__all__ = ["__version__"]
