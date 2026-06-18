"""Unit tests for the CIGAR walker, anchor-window math, and per-(read, intron)
classifier.  These pin the semantics described in irfinder_mdl/quantify.py."""

from __future__ import annotations

import pytest

from irfinder_mdl.quantify import (
    CIGAR_M, CIGAR_I, CIGAR_D, CIGAR_N, CIGAR_S, CIGAR_H,
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

    def test_insufficient_anchor_on_exon_side_rejected(self):
        # Only 7 bp on exon side (993..1000), need 8.  Read end at 1100 also
        # fails the right-boundary anchor.  Result: no boundary supported,
        # classifier returns None.
        matched = [(993, 1100)]
        n_skips = []
        f = _classify(matched, n_skips)
        assert f is None

    def test_insufficient_anchor_on_intron_side_rejected(self):
        # Plenty of exon, only 7 bp into the intron
        matched = [(900, 1007)]
        n_skips = []
        f = _classify(matched, n_skips)
        # crosses_left? not enough anchor inside intron -> None entire read
        assert f is None

    def test_deletion_inside_anchor_window_breaks_retention(self):
        # 2 bp deletion straddling the boundary makes the intronic-side anchor
        # only 6/8 bp.  Strict policy: this read is not a retention witness.
        # Plenty of read context elsewhere; the failure is the anchor itself.
        matched = [(900, 1000), (1002, 1100)]
        f = _classify(matched, [])
        assert f is None  # no boundary crossing reliably supported

    def test_deletion_outside_anchor_window_tolerated(self):
        # 5 bp deletion well inside the intron, far from any anchor window.
        # The boundary is still cleanly covered with 8+8 bp.
        matched = [(900, 1500), (1505, 2100)]
        f = _classify(matched, [])
        assert f["retain_left"] is True
        assert f["retain_right"] is True


class TestClassifierMixedAndExclusion:
    def test_read_entirely_inside_intron_is_rejected(self):
        matched = [(1200, 1500)]
        n_skips = []
        f = _classify(matched, n_skips)
        assert f is None

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
