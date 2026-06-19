"""Junction-anchored intron retention counting from a sorted+indexed BAM.

For each annotated intron at (s, e) on a chromosome (0-based half-open):

  * SPLICE at the 5' boundary  -- a read has a CIGAR `N` block whose **start**
    is within `jitter` bp of `s`, supported by at least `anchor` matched bp
    immediately upstream (5' exonic side).
  * SPLICE at the 3' boundary  -- analogous: an `N` block whose **end** is
    within `jitter` of `e`, supported by `anchor` matched bp immediately
    downstream (3' exonic side).
  * RETENTION at the 5' boundary -- the read aligns continuously through
    position `s`, with at least `anchor` matched bp on the exonic side
    (`[s-anchor, s)`) AND at least `anchor` matched bp on the intronic side
    (`[s, s+anchor)`), and no `N` block crosses the anchor window.
  * RETENTION at the 3' boundary -- analogous around `e`.

Reads that lie entirely inside an intron (do not cross either boundary with the
required anchor) contribute nothing.  This is the explicit ask: no fully-
intronic "coverage" reads are taken as evidence of retention.

A single read may contribute to several intron counts (long reads frequently
span several genes).  Within one (read, intron) pair the read can land in
multiple boolean buckets at once -- e.g. spliced at one end and retained at
the other -- and each bucket counts independently.
"""

from __future__ import annotations

import multiprocessing as mp
import os
import re
import sys
from collections import defaultdict
from dataclasses import dataclass
from typing import Iterable

import pysam
from intervaltree import IntervalTree

from .gtf import Intron


# ---------------------------------------------------------------------------
# CIGAR opcodes (pysam.AlignedSegment.cigartuples convention)
# ---------------------------------------------------------------------------
CIGAR_M, CIGAR_I, CIGAR_D, CIGAR_N, CIGAR_S, CIGAR_H, CIGAR_P, CIGAR_EQ, CIGAR_X = range(9)
_MATCHED = frozenset({CIGAR_M, CIGAR_EQ, CIGAR_X})


# ---------------------------------------------------------------------------
# Primary-assembly chromosome filter
# ---------------------------------------------------------------------------
# Matches the nuclear primary chromosomes in both UCSC/GENCODE (`chr1`, `chrX`)
# and Ensembl (`1`, `X`) naming, and *excludes*:
#   - the mitochondrion (`chrM`, `chrMT`, `M`, `MT`) -- it has no spliced
#     introns, so it contributes nothing but noise to IR;
#   - every unplaced/alt contig: GENCODE/Ensembl bare contigs (`GL000008.2`,
#     `KI270442.1`) and UCSC decoys (`chrUn_*`, `chr1_KI270706v1_random`,
#     `chr19_KI270890v1_alt`).
# The optional `chr` prefix plus the `$`-anchored alternation is what rejects
# the decoys (their `_`/`Un` suffixes never match).
_PRIMARY_CHROM_RE = re.compile(r"^(?:chr)?(?:[0-9]+|X|Y)$")


def is_primary_chrom(chrom: str) -> bool:
    """True for nuclear primary-assembly chromosomes (chr1-N, X, Y, with or
    without the `chr` prefix); False for the mitochondrion and every
    unplaced/alt contig."""
    return _PRIMARY_CHROM_RE.match(chrom) is not None


# ---------------------------------------------------------------------------
# Parameters and counts
# ---------------------------------------------------------------------------
@dataclass(slots=True, frozen=True)
class QuantParams:
    anchor: int = 8       # bp of matched alignment required on each anchor side
    jitter: int = 3       # bp tolerance for an N op aligning to an intron boundary
    min_mapq: int = 1     # exclude MAPQ=0 multimappers by default
    # Default = 0x900 = secondary(0x100) | supplementary(0x800).  Unmapped
    # reads are already excluded via the explicit `is_unmapped` check below;
    # callers may add 0x400 (duplicate) or 0x200 (QC fail) if they want.
    exclude_flags: int = 0x900


@dataclass(slots=True)
class IntronCounts:
    """Counts populated per annotated intron.

    `crossing_reads` is the unique-read total: every read that contributed at
    least one of the boolean signals below increments it once.  Most other
    fields are independent boolean tallies over the same read set, so a read
    may contribute to several of them simultaneously.

    The `interior_reads` / `intron_coverage_bp` fields carry the additional
    intronic-depth signal added in v0.0.1 (see README "Depth-augmented mode"):
    they let downstream code compute an IRFinder-S-style IR ratio that
    incorporates reads sitting inside the intron without requiring them to
    reach a boundary.
    """
    crossing_reads: int = 0
    splice_left: int = 0
    splice_right: int = 0
    splice_exact: int = 0
    retain_left: int = 0
    retain_right: int = 0
    retain_both: int = 0   # crosses both boundaries AND retains both
    mixed: int = 0         # splices one boundary, retains the other

    # Depth-augmented signal (no boundary required).
    interior_reads: int = 0      # reads fully inside [s, e) with no N op intersecting [s, e)
    intron_coverage_bp: int = 0  # sum of M bp inside [s, e) from reads that don't splice through this intron


# ---------------------------------------------------------------------------
# CIGAR walker
# ---------------------------------------------------------------------------
def parse_cigar(reference_start: int, cigartuples):
    """Return (matched_intervals, n_intervals) for an aligned read.

    `matched_intervals` is a list of `(ref_start, ref_end)` 0-based half-open
    intervals where the read carries aligned bases (M/=/X).  `n_intervals` is
    the same shape for `N` (skipped reference) ops.  Both lists are sorted
    ascending by `ref_start` because CIGAR ops are in alignment order.

    Insertions (I), soft-clips (S), hard-clips (H), and padding (P) consume
    read bases without advancing the reference, so they don't appear in the
    output.  Deletions (D) advance the reference but are not "matched" -- they
    create a gap in `matched_intervals`.
    """
    matched: list[tuple[int, int]] = []
    n_skips: list[tuple[int, int]] = []
    ref_pos = reference_start
    for op, length in cigartuples or ():
        if op in _MATCHED:
            matched.append((ref_pos, ref_pos + length))
            ref_pos += length
        elif op == CIGAR_N:
            n_skips.append((ref_pos, ref_pos + length))
            ref_pos += length
        elif op == CIGAR_D:
            ref_pos += length
        # I, S, H, P: no reference advance
    return matched, n_skips


def matched_bp_in(window_start: int, window_end: int, matched_intervals) -> int:
    """Number of reference positions in [window_start, window_end) covered by
    `matched_intervals`.  Linear in the number of overlapping intervals; for
    long-read CIGARs this is typically 1-2."""
    total = 0
    for ms, me in matched_intervals:
        if me <= window_start:
            continue
        if ms >= window_end:
            break  # matched_intervals is sorted by start
        total += min(me, window_end) - max(ms, window_start)
    return total


# ---------------------------------------------------------------------------
# Per-(read, intron) classifier
# ---------------------------------------------------------------------------
def classify_read_vs_intron(
    intron_start: int,
    intron_end: int,
    matched_intervals,
    n_skips,
    anchor: int,
    jitter: int,
):
    """Return a dict of signals contributed by this read against this intron,
    or None if the read contributes nothing.

    A read contributes if **any** of the following is true:
      * it splices either boundary (anchor-passing N op);
      * it retains either boundary (anchor-passing M alignment);
      * it carries at least one matched bp inside `[s, e)` AND no `N` op
        intersects `[s, e)` -- this is the depth-augmented signal that lets
        fully-interior reads (and partial-overlap retention reads) contribute
        to `intron_coverage_bp`.
    """
    if not matched_intervals:
        return None
    read_start = matched_intervals[0][0]
    read_end = matched_intervals[-1][1]
    s, e = intron_start, intron_end
    k, j = anchor, jitter

    # Fast reject: read is entirely outside the intron + anchor window.
    if read_end <= s - k or read_start >= e + k:
        return None

    # Splice at each boundary
    splice_left = splice_right = splice_exact = False
    for ns_start, ns_end in n_skips:
        sl_ok = (
            abs(ns_start - s) <= j
            and matched_bp_in(ns_start - k, ns_start, matched_intervals) >= k
        )
        sr_ok = (
            abs(ns_end - e) <= j
            and matched_bp_in(ns_end, ns_end + k, matched_intervals) >= k
        )
        if sl_ok:
            splice_left = True
        if sr_ok:
            splice_right = True
        if sl_ok and sr_ok:
            splice_exact = True  # one N op covers both ends of this intron

    # Retention at left boundary `s`
    retain_left = False
    if read_start <= s - k and read_end >= s + k:
        if (
            matched_bp_in(s - k, s, matched_intervals) >= k
            and matched_bp_in(s, s + k, matched_intervals) >= k
        ):
            if not any(ns < s + k and ne > s - k for ns, ne in n_skips):
                retain_left = True

    # Retention at right boundary `e`
    retain_right = False
    if read_start <= e - k and read_end >= e + k:
        if (
            matched_bp_in(e - k, e, matched_intervals) >= k
            and matched_bp_in(e, e + k, matched_intervals) >= k
        ):
            if not any(ns < e + k and ne > e - k for ns, ne in n_skips):
                retain_right = True

    # Depth-augmented signal: M bp inside [s, e), gated on the read not
    # using this intron as a splice donor/acceptor.  Reads that splice an
    # *unrelated* intron entirely outside [s, e) still contribute their
    # intronic M bp.
    n_intersects_intron = any(ns < e and ne > s for ns, ne in n_skips)
    intronic_bp = 0 if n_intersects_intron else matched_bp_in(s, e, matched_intervals)
    # Interior: read's alignment is fully inside [s, e), no N intersecting.
    interior = (
        not n_intersects_intron
        and read_start >= s
        and read_end <= e
        and intronic_bp > 0
    )

    crosses_left = splice_left or retain_left
    crosses_right = splice_right or retain_right
    if not (crosses_left or crosses_right) and intronic_bp == 0:
        # Read overlaps [s-k, e+k] but contributes no signal -- e.g. sits in
        # the upstream/downstream flank without crossing a boundary.
        return None

    return {
        "splice_left": splice_left,
        "splice_right": splice_right,
        "splice_exact": splice_exact,
        "retain_left": retain_left,
        "retain_right": retain_right,
        "interior": interior,
        "intronic_bp": intronic_bp,
    }


def update_counts(counts: IntronCounts, f: dict) -> None:
    counts.crossing_reads += 1
    counts.splice_left  += f["splice_left"]
    counts.splice_right += f["splice_right"]
    counts.splice_exact += f["splice_exact"]
    counts.retain_left  += f["retain_left"]
    counts.retain_right += f["retain_right"]
    if f["retain_left"] and f["retain_right"]:
        counts.retain_both += 1
    elif (f["retain_left"] and f["splice_right"]) or (f["retain_right"] and f["splice_left"]):
        counts.mixed += 1
    if f["interior"]:
        counts.interior_reads += 1
    counts.intron_coverage_bp += f["intronic_bp"]

# ---------------------------------------------------------------------------
# Per-chromosome worker
# ---------------------------------------------------------------------------
def quantify_chromosome(
    bam_path: str,
    chrom: str,
    introns: list[Intron],
    params: QuantParams,
) -> list[IntronCounts]:
    """Count over every aligned read on one chromosome.

    Reads are streamed in coordinate order via `pysam.fetch(chrom)`, so this
    function never materialises the full BAM in memory; only one read at a
    time, plus the intron index for this chromosome.
    """
    bam = pysam.AlignmentFile(bam_path, "rb")
    try:
        tree = IntervalTree()
        for idx, intron in enumerate(introns):
            if intron.end > intron.start:
                tree.addi(intron.start, intron.end, idx)
        counts: list[IntronCounts] = [IntronCounts() for _ in introns]

        for read in bam.fetch(chrom):
            if read.flag & params.exclude_flags:
                continue
            if read.is_unmapped:
                continue
            if read.mapping_quality < params.min_mapq:
                continue

            matched, n_skips = parse_cigar(read.reference_start, read.cigartuples)
            if not matched:
                continue
            r_start = matched[0][0]
            r_end = matched[-1][1]

            # Pad query by anchor so we don't miss reads that just barely
            # reach an intron boundary from outside the intron interval.
            pad = params.anchor
            # IntervalTree.overlap(begin, end): items with begin < end_q AND end > begin_q
            for iv in tree.overlap(r_start - pad, r_end + pad):
                idx = iv.data
                intron = introns[idx]
                flags = classify_read_vs_intron(
                    intron.start, intron.end,
                    matched, n_skips,
                    params.anchor, params.jitter,
                )
                if flags is None:
                    continue
                update_counts(counts[idx], flags)
        return counts
    finally:
        bam.close()


def _worker(args):
    bam_path, chrom, introns, params = args
    return chrom, quantify_chromosome(bam_path, chrom, introns, params)


def quantify_introns(
    bam_path: str,
    introns: list[Intron],
    params: QuantParams,
    threads: int = 1,
    *,
    chroms: Iterable[str] | None = None,
    log=print,
) -> list[IntronCounts]:
    """Quantify every intron across the chromosomes available in the BAM.

    Output order matches `introns`.  Introns on chromosomes not in the BAM
    header (e.g. unplaced contigs absent from the alignment) are returned with
    zero counts.
    """
    bam = pysam.AlignmentFile(bam_path, "rb")
    try:
        bam_chroms = set(bam.references)
    finally:
        bam.close()

    by_chrom: dict[str, list[Intron]] = defaultdict(list)
    idx_by_chrom: dict[str, list[int]] = defaultdict(list)
    for i, intron in enumerate(introns):
        by_chrom[intron.chrom].append(intron)
        idx_by_chrom[intron.chrom].append(i)

    requested = set(chroms) if chroms is not None else set(by_chrom)
    targets = [c for c in by_chrom if c in requested and c in bam_chroms]
    # Sort by intron count desc -> longer chromosomes start first under a pool.
    targets.sort(key=lambda c: -len(by_chrom[c]))
    skipped = sorted(set(by_chrom) - set(targets))
    if skipped:
        log(f"[quant] skipping {len(skipped)} chrom(s) absent from BAM or filtered: "
            f"{skipped[:6]}{'...' if len(skipped) > 6 else ''}", file=sys.stderr)

    counts: list[IntronCounts] = [IntronCounts() for _ in introns]
    tasks = [(bam_path, c, by_chrom[c], params) for c in targets]

    if threads <= 1 or len(tasks) <= 1:
        for task in tasks:
            chrom, c_list = _worker(task)
            for cnt, idx in zip(c_list, idx_by_chrom[chrom]):
                counts[idx] = cnt
            log(f"[quant] {chrom} done ({len(c_list):,} introns)", file=sys.stderr)
    else:
        ctx = mp.get_context("spawn")
        with ctx.Pool(min(threads, len(tasks))) as pool:
            for chrom, c_list in pool.imap_unordered(_worker, tasks):
                for cnt, idx in zip(c_list, idx_by_chrom[chrom]):
                    counts[idx] = cnt
                log(f"[quant] {chrom} done ({len(c_list):,} introns)", file=sys.stderr)

    return counts


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------
QUANT_TSV_HEADER = [
    "chrom", "start", "end", "strand", "intron_id",
    "gene_ids", "gene_names", "gene_types",
    "intron_length", "exon_overlap",
    "crossing_reads",
    "splice_left", "splice_right", "splice_exact",
    "retain_left", "retain_right", "retain_both", "mixed",
    "interior_reads", "intron_coverage_bp", "intron_mean_depth",
    "ir_ratio_left", "ir_ratio_right",
    "ir_ratio_junction",  # junction-anchored only
    "ir_ratio_depth",     # IRFinder-S-style: depth competes with max splice
]


def _safe_ratio(retention: float, splice: float) -> str:
    """Return the IR ratio `retention / (retention + splice)` as a fixed-
    precision string, or '.' when there is no splice evidence to compare
    against (`splice == 0`).

    The denominator's `splice` term is required to be non-zero by design: an
    IR ratio is a *relative* measure of retention vs splicing, and if no
    splicing was observed at this intron / boundary we cannot tell whether
    that is because the intron is fully retained or simply because the locus
    is not expressed enough to produce splice evidence.  Outputting '.' for
    these cases prevents lowly-expressed introns with stray interior reads
    from artificially inflating the IR distribution at 1.0.
    """
    if splice <= 0:
        return "."
    return f"{retention / (retention + splice):.6f}"


def _fmt_depth(depth: float) -> str:
    return f"{depth:.4f}"


def write_quant_tsv(
    introns: list[Intron],
    counts: list[IntronCounts],
    path: str,
) -> None:
    import gzip
    from .gtf import intron_id
    fh = gzip.open(path, "wt") if path.endswith(".gz") else open(path, "wt")
    try:
        fh.write("\t".join(QUANT_TSV_HEADER) + "\n")
        for intron, c in zip(introns, counts):
            intron_length = intron.end - intron.start
            mean_depth = (
                c.intron_coverage_bp / intron_length if intron_length > 0 else 0.0
            )
            max_splice = max(c.splice_left, c.splice_right)
            ir_left  = _safe_ratio(c.retain_left,  c.splice_left)
            ir_right = _safe_ratio(c.retain_right, c.splice_right)
            ir_junction = _safe_ratio(
                c.retain_left + c.retain_right,
                c.splice_left + c.splice_right,
            )
            # IRFinder-S form: mean intron depth competes with the larger of
            # the two splice-junction counts.
            ir_depth = _safe_ratio(mean_depth, max_splice)
            fh.write("\t".join([
                intron.chrom,
                str(intron.start),
                str(intron.end),
                intron.strand,
                intron_id(intron),
                ",".join(sorted(intron.gene_ids)) if intron.gene_ids else "",
                ",".join(sorted(intron.gene_names)) if intron.gene_names else "",
                ",".join(sorted(intron.gene_types)) if intron.gene_types else "",
                str(intron_length),
                "1" if intron.exon_overlap else "0",
                str(c.crossing_reads),
                str(c.splice_left), str(c.splice_right), str(c.splice_exact),
                str(c.retain_left), str(c.retain_right),
                str(c.retain_both), str(c.mixed),
                str(c.interior_reads), str(c.intron_coverage_bp),
                _fmt_depth(mean_depth),
                ir_left, ir_right,
                ir_junction, ir_depth,
            ]) + "\n")
    finally:
        fh.close()
