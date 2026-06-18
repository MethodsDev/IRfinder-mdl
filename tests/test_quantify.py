"""Unit tests for the CIGAR walker, anchor-window math, and per-(read, intron)
classifier.  These pin the semantics described in irfinder_mdl/quantify.py."""

from __future__ import annotations

import pytest

from irfinder_mdl.quantify import (
    CIGAR_M, CIGAR_I, CIGAR_D, CIGAR_N, CIGAR_S, CIGAR_H,
    _safe_ratio,
    classify_read_vs_intron,
    matched_bp_in,
    parse_cigar,
)


# ---------------------------------------------------------------------------
# parse_cigar
# ---------------------------------------------------------------------------
class TestParseCigar:
    def test_simple_match(self):
        matched, n = parse_cigar(100, [(CIGAR_M, 50)])
        assert matched == [(100, 150)]
        assert n == []

    def test_match_skip_match(self):
        matched, n = parse_cigar(100, [(CIGAR_M, 50), (CIGAR_N, 1000), (CIGAR_M, 30)])
        assert matched == [(100, 150), (1150, 1180)]
        assert n == [(150, 1150)]

    def test_soft_clip_does_not_consume_reference(self):
        matched, n = parse_cigar(100, [(CIGAR_S, 20), (CIGAR_M, 50)])
        assert matched == [(100, 150)]
        assert n == []

    def test_insertion_does_not_consume_reference(self):
        matched, n = parse_cigar(
            100, [(CIGAR_M, 20), (CIGAR_I, 5), (CIGAR_M, 30)]
        )
        assert matched == [(100, 120), (120, 150)]
        assert n == []

    def test_deletion_consumes_reference_but_not_matched(self):
        matched, n = parse_cigar(
            100, [(CIGAR_M, 20), (CIGAR_D, 5), (CIGAR_M, 30)]
        )
        # D advances ref but matched intervals remain disjoint
        assert matched == [(100, 120), (125, 155)]
        assert n == []

    def test_hard_clip_does_not_consume_reference(self):
        matched, n = parse_cigar(100, [(CIGAR_H, 20), (CIGAR_M, 50)])
        assert matched == [(100, 150)]

    def test_two_introns(self):
        matched, n = parse_cigar(
            100,
            [(CIGAR_M, 50), (CIGAR_N, 200), (CIGAR_M, 30),
             (CIGAR_N, 500), (CIGAR_M, 40)],
        )
        assert matched == [(100, 150), (350, 380), (880, 920)]
        assert n == [(150, 350), (380, 880)]


# ---------------------------------------------------------------------------
# matched_bp_in
# ---------------------------------------------------------------------------
class TestMatchedBpIn:
    def test_fully_inside(self):
        assert matched_bp_in(110, 120, [(100, 150)]) == 10

    def test_window_outside(self):
        assert matched_bp_in(200, 210, [(100, 150)]) == 0

    def test_partial_overlap_left(self):
        assert matched_bp_in(90, 110, [(100, 150)]) == 10

    def test_partial_overlap_right(self):
        assert matched_bp_in(140, 160, [(100, 150)]) == 10

    def test_multiple_intervals(self):
        intervals = [(100, 150), (200, 250)]
        # window covers part of both
        assert matched_bp_in(140, 210, intervals) == 10 + 10

    def test_gap_in_window_returns_only_covered(self):
        intervals = [(100, 150), (200, 250)]
        # window entirely inside the gap
        assert matched_bp_in(160, 190, intervals) == 0

    def test_empty(self):
        assert matched_bp_in(100, 200, []) == 0


# ---------------------------------------------------------------------------
# classify_read_vs_intron
# ---------------------------------------------------------------------------
# Convention used in these tests: intron at [1000, 2000) (1 kb), anchor=8,
# jitter=3.  Reads are described by their CIGAR string-equivalent.
INTRON_S, INTRON_E = 1000, 2000
K, J = 8, 3


def _classify(matched, n_skips, *, intron=(INTRON_S, INTRON_E), k=K, j=J):
    return classify_read_vs_intron(intron[0], intron[1], matched, n_skips, k, j)


class TestClassifierSplice:
    def test_exact_splice_both_boundaries(self):
        # Read: 100M_1000N_100M starting at 900 -> intron at [1000, 2000)
        matched = [(900, 1000), (2000, 2100)]
        n_skips = [(1000, 2000)]
        f = _classify(matched, n_skips)
        assert f["splice_left"] is True
        assert f["splice_right"] is True
        assert f["splice_exact"] is True
        assert f["retain_left"] is False
        assert f["retain_right"] is False

    def test_splice_within_jitter(self):
        # N starts at 1002 (jitter=2 <= 3) and ends at 2001 (jitter=1)
        matched = [(900, 1002), (2001, 2100)]
        n_skips = [(1002, 2001)]
        f = _classify(matched, n_skips)
        assert f["splice_left"] is True
        assert f["splice_right"] is True
        assert f["splice_exact"] is True

    def test_splice_outside_jitter_rejected(self):
        # N starts at 1005 (off by 5 > 3) -> no splice_left
        matched = [(900, 1005), (2000, 2100)]
        n_skips = [(1005, 2000)]
        f = _classify(matched, n_skips)
        assert f["splice_left"] is False
        assert f["splice_right"] is True
        assert f["splice_exact"] is False

    def test_insufficient_left_anchor(self):
        # Only 7 bp of M before the N (k=8 required)
        matched = [(993, 1000), (2000, 2100)]
        n_skips = [(1000, 2000)]
        f = _classify(matched, n_skips)
        assert f["splice_left"] is False
        assert f["splice_right"] is True

    def test_insufficient_right_anchor(self):
        # Only 7 bp of M after the N
        matched = [(900, 1000), (2000, 2007)]
        n_skips = [(1000, 2000)]
        f = _classify(matched, n_skips)
        assert f["splice_left"] is True
        assert f["splice_right"] is False
        assert f["splice_exact"] is False

    def test_truncated_read_only_left_junction(self):
        # Read ends right after the N: 100M_1000N (no downstream exon at all)
        matched = [(900, 1000)]
        n_skips = [(1000, 2000)]
        f = _classify(matched, n_skips)
        assert f["splice_left"] is True
        assert f["splice_right"] is False  # no M anchor after N


class TestClassifierRetention:
    def test_full_retention_through_intron(self):
        # One big M block from before to after
        matched = [(900, 2100)]
        n_skips = []
        f = _classify(matched, n_skips)
        assert f["retain_left"] is True
        assert f["retain_right"] is True
        assert f["splice_left"] is False
        assert f["splice_right"] is False

    def test_retain_left_only_truncated_inside_intron(self):
        # Read covers the 5' exon and the first half of the intron
        matched = [(900, 1500)]
        n_skips = []
        f = _classify(matched, n_skips)
        assert f["retain_left"] is True
        assert f["retain_right"] is False

    def test_retain_right_only(self):
        # Read starts inside the intron and continues into the 3' exon
        matched = [(1500, 2100)]
        n_skips = []
        f = _classify(matched, n_skips)
        assert f["retain_left"] is False
        assert f["retain_right"] is True

    def test_insufficient_anchor_on_exon_side_no_retain(self):
        # 7 bp on exon side (need 8) -> retain_left False, but the 100 bp
        # inside the intron still count toward intron_coverage_bp.  The read
        # is not "interior" (its alignment starts in the upstream exon).
        matched = [(993, 1100)]
        f = _classify(matched, [])
        assert f is not None
        assert f["retain_left"] is False
        assert f["retain_right"] is False
        assert f["interior"] is False
        assert f["intronic_bp"] == 100

    def test_insufficient_anchor_on_intron_side_no_retain(self):
        # Only 7 bp into the intron (need 8) -> retain_left False, but the
        # 7 intronic bp still contribute to coverage.
        matched = [(900, 1007)]
        f = _classify(matched, [])
        assert f is not None
        assert f["retain_left"] is False
        assert f["intronic_bp"] == 7
        assert f["interior"] is False

    def test_deletion_inside_anchor_window_breaks_retention(self):
        # 2 bp deletion straddling the boundary defeats the strict retain
        # check (only 6 matched bp in the intronic anchor window).  Depth
        # signal still counts the 98 matched bp inside the intron.
        matched = [(900, 1000), (1002, 1100)]
        f = _classify(matched, [])
        assert f is not None
        assert f["retain_left"] is False
        assert f["intronic_bp"] == 98

    def test_deletion_outside_anchor_window_tolerated(self):
        # 5 bp deletion well inside the intron, far from any anchor window.
        # The boundary is still cleanly covered with 8+8 bp.
        matched = [(900, 1500), (1505, 2100)]
        f = _classify(matched, [])
        assert f["retain_left"] is True
        assert f["retain_right"] is True
        # Intronic M bp = (1500-1000) + (2000-1505) = 500 + 495 = 995
        assert f["intronic_bp"] == 995


class TestClassifierMixedAndExclusion:
    def test_read_entirely_inside_intron_marks_interior(self):
        # Fully inside [1000, 2000) with no N -> interior=True, contributes
        # its full length to intron_coverage_bp, but no boundary signals.
        matched = [(1200, 1500)]
        f = _classify(matched, [])
        assert f is not None
        assert f["interior"] is True
        assert f["intronic_bp"] == 300
        assert f["splice_left"] is False
        assert f["splice_right"] is False
        assert f["retain_left"] is False
        assert f["retain_right"] is False

    def test_mixed_splice_left_retain_right(self):
        # CIGAR: 100M starting at 900 to 1000, then N from 1000 to 1800,
        # then M from 1800 to 2100.  The read splices the LEFT boundary
        # (N matches s=1000) but enters the intron region from the 3' side
        # via a M block that includes position 1999 (last intronic) AND 2000
        # (first 3' exon).  Since the M block is [1800, 2100), positions
        # in [e-k, e) = [1992, 2000) are matched (8 bp), and [e, e+k) =
        # [2000, 2008) are matched.  No N op crosses the right boundary.
        # So retain_right = True; splice_left = True.
        matched = [(900, 1000), (1800, 2100)]
        n_skips = [(1000, 1800)]
        f = _classify(matched, n_skips)
        assert f["splice_left"] is True
        assert f["splice_right"] is False
        # The N op ends at 1800, NOT within jitter of 2000.
        assert f["retain_right"] is True

    def test_read_far_away_rejected(self):
        # A read on the other side of the chromosome
        matched = [(10000, 11000)]
        f = _classify(matched, [])
        assert f is None

# ---------------------------------------------------------------------------
# Depth-augmented signal (interior reads + intron_coverage_bp)
# ---------------------------------------------------------------------------
class TestDepthSignal:
    def test_clean_splice_contributes_zero_depth(self):
        # A read that splices the intron exactly contributes 0 intronic_bp.
        matched = [(900, 1000), (2000, 2100)]
        n_skips = [(1000, 2000)]
        f = _classify(matched, n_skips)
        assert f["splice_exact"] is True
        assert f["intronic_bp"] == 0
        assert f["interior"] is False

    def test_partial_overlap_retention_contributes_partial_depth(self):
        # Read covers 500 bp inside the intron via retention at the left
        # boundary; no N; counts those 500 bp toward intron_coverage_bp.
        matched = [(900, 1500)]
        f = _classify(matched, [])
        assert f["retain_left"] is True
        assert f["intronic_bp"] == 500
        assert f["interior"] is False

    def test_full_retention_contributes_full_intron_length(self):
        matched = [(900, 2100)]
        f = _classify(matched, [])
        assert f["retain_left"] is True
        assert f["retain_right"] is True
        # intron length = 1000
        assert f["intronic_bp"] == 1000

    def test_interior_read_with_internal_N_does_not_count(self):
        # Read fully inside the intron but with an N inside it -- e.g. a
        # spurious internal splice.  We refuse to credit this as retention
        # depth.
        matched = [(1100, 1300), (1400, 1700)]
        n_skips = [(1300, 1400)]
        f = _classify(matched, n_skips)
        # n_intersects_intron -> intronic_bp = 0, interior = False
        assert f is None  # no other signal -> no contribution

    def test_read_spliced_at_unrelated_intron_still_contributes_depth(self):
        # Read has a CIGAR `100M_50N_500M` where the N is upstream of our
        # intron entirely.  Inside our intron the read is contiguous M and no
        # N intersects [s, e), so 500 bp of intronic depth count.
        matched = [(800, 900), (950, 1450)]
        n_skips = [(900, 950)]
        s, e = 1000, 2000
        f = classify_read_vs_intron(s, e, matched, n_skips, anchor=8, jitter=3)
        # No N intersects [1000, 2000); intronic_bp = matched bp in [1000, 1450) = 450
        assert f is not None
        assert f["intronic_bp"] == 450

    def test_n_intersecting_intron_zeros_depth(self):
        # Read with an N op that partially overlaps the intron disqualifies
        # all intronic bp contribution from this read.
        matched = [(900, 1100), (1200, 1500)]
        n_skips = [(1100, 1200)]
        f = _classify(matched, n_skips)
        # n_intersects_intron at [1100, 1200) intersects [1000, 2000)
        assert f is not None  # may still have a boundary signal at left
        assert f["intronic_bp"] == 0


# ---------------------------------------------------------------------------
# Aggregation: update_counts wires everything into IntronCounts
# ---------------------------------------------------------------------------
class TestUpdateCounts:
    def test_interior_read_increments_interior_and_bp(self):
        from irfinder_mdl.quantify import IntronCounts, update_counts
        c = IntronCounts()
        flags = {"splice_left": False, "splice_right": False, "splice_exact": False,
                 "retain_left": False, "retain_right": False,
                 "interior": True, "intronic_bp": 300}
        update_counts(c, flags)
        assert c.crossing_reads == 1
        assert c.interior_reads == 1
        assert c.intron_coverage_bp == 300
        assert c.splice_left == c.splice_right == 0
        assert c.retain_left == c.retain_right == 0

    def test_boundary_retain_read_is_not_interior_but_contributes_bp(self):
        from irfinder_mdl.quantify import IntronCounts, update_counts
        c = IntronCounts()
        flags = {"splice_left": False, "splice_right": False, "splice_exact": False,
                 "retain_left": True, "retain_right": False,
                 "interior": False, "intronic_bp": 500}
        update_counts(c, flags)
        assert c.retain_left == 1
        assert c.interior_reads == 0
        assert c.intron_coverage_bp == 500


# ---------------------------------------------------------------------------
# IR ratio semantics: undefined when no splice evidence
# ---------------------------------------------------------------------------
class TestSafeRatio:
    def test_typical_ratio(self):
        # 30 retain, 70 splice -> 30/100 = 0.30
        assert _safe_ratio(30, 70) == "0.300000"

    def test_no_evidence_is_dot(self):
        assert _safe_ratio(0, 0) == "."

    def test_no_splice_is_dot_even_with_retention(self):
        # No splice comparator -- ratio is undefined, not 1.0.
        assert _safe_ratio(50, 0) == "."

    def test_no_retention_is_zero(self):
        # Splice present, retention zero -> 0.0
        assert _safe_ratio(0, 50) == "0.000000"

    def test_works_with_floats(self):
        # Mean depth (float) competes with splice count (int)
        assert _safe_ratio(2.5, 7.5) == "0.250000"

# ---------------------------------------------------------------------------
# GTF parsing roundtrip
# ---------------------------------------------------------------------------
class TestGtfIntronBuild:
    def _write_gtf(self, tmp_path, body: str):
        p = tmp_path / "mini.gtf"
        p.write_text(body)
        return str(p)

    def test_three_exon_transcript_yields_two_introns(self, tmp_path):
        gtf = self._write_gtf(tmp_path, (
            'chr1\tHAVANA\texon\t101\t200\t.\t+\t.\t'
            'gene_id "G1"; transcript_id "T1"; gene_name "GENE1"; '
            'gene_type "protein_coding"; transcript_type "protein_coding";\n'
            'chr1\tHAVANA\texon\t301\t400\t.\t+\t.\t'
            'gene_id "G1"; transcript_id "T1"; gene_name "GENE1";\n'
            'chr1\tHAVANA\texon\t501\t600\t.\t+\t.\t'
            'gene_id "G1"; transcript_id "T1"; gene_name "GENE1";\n'
        ))
        from irfinder_mdl.gtf import build_introns_from_gtf
        introns = build_introns_from_gtf(gtf, log=lambda *a, **kw: None)
        assert len(introns) == 2
        # GTF 1-based inclusive [101, 200] -> BED [100, 200)
        # Intron between exons 1 and 2: [200, 300)
        # Intron between exons 2 and 3: [400, 500)
        assert (introns[0].chrom, introns[0].start, introns[0].end) == ("chr1", 200, 300)
        assert (introns[1].chrom, introns[1].start, introns[1].end) == ("chr1", 400, 500)
        assert introns[0].gene_ids == {"G1"}
        assert introns[0].transcript_ids == {"T1"}

    def test_shared_intron_merged_across_transcripts(self, tmp_path):
        gtf = self._write_gtf(tmp_path, (
            # T1: exons [101, 200] and [301, 400] -> intron [200, 300)
            'chr1\tHAVANA\texon\t101\t200\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n'
            'chr1\tHAVANA\texon\t301\t400\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n'
            # T2: same intron, different exon starts/ends but same junction
            'chr1\tHAVANA\texon\t150\t200\t.\t+\t.\tgene_id "G1"; transcript_id "T2";\n'
            'chr1\tHAVANA\texon\t301\t450\t.\t+\t.\tgene_id "G1"; transcript_id "T2";\n'
        ))
        from irfinder_mdl.gtf import build_introns_from_gtf
        introns = build_introns_from_gtf(gtf, log=lambda *a, **kw: None)
        assert len(introns) == 1
        assert introns[0].transcript_ids == {"T1", "T2"}

    def test_exon_overlap_flag(self, tmp_path):
        # T1: exons [101, 200] and [301, 400] -> intron [200, 300)
        # T2: a single exon at [101, 400] that spans T1's intron -> overlap!
        gtf = self._write_gtf(tmp_path, (
            'chr1\tHAVANA\texon\t101\t200\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n'
            'chr1\tHAVANA\texon\t301\t400\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n'
            'chr1\tHAVANA\texon\t101\t400\t.\t+\t.\tgene_id "G1"; transcript_id "T2";\n'
        ))
        from irfinder_mdl.gtf import build_introns_from_gtf
        introns = build_introns_from_gtf(gtf, log=lambda *a, **kw: None)
        assert len(introns) == 1
        assert introns[0].exon_overlap is True

    def test_clean_intron_no_exon_overlap(self, tmp_path):
        gtf = self._write_gtf(tmp_path, (
            'chr1\tHAVANA\texon\t101\t200\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n'
            'chr1\tHAVANA\texon\t301\t400\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n'
        ))
        from irfinder_mdl.gtf import build_introns_from_gtf
        introns = build_introns_from_gtf(gtf, log=lambda *a, **kw: None)
        assert len(introns) == 1
        assert introns[0].exon_overlap is False
