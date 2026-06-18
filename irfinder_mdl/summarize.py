"""Aggregate a per-intron quant TSV into genome- and chromosome-level metrics.

Tracks two complementary IR signals:

  * **Junction-anchored** (`ir_ratio_junction`): retain reads at the two
    boundaries vs splice reads at the two boundaries.  Conservative; the
    cleanest signal for short introns and any analysis where you do not
    want to credit reads that do not reach a boundary.
  * **Depth-augmented** (`ir_ratio_depth`): IRFinder-S form
    `mean_intron_depth / (mean_intron_depth + max_splice)`.  Includes reads
    that lie fully inside the intron (`interior_reads`) and partial-overlap
    retention reads' intronic bp.  Less biased against long introns.

Both summaries are computed in one pass.
"""

from __future__ import annotations

import gzip
import json
import sys
from dataclasses import asdict, dataclass, field
from typing import IO


@dataclass
class IRStats:
    """Quantiles + thresholded counts for one IR-ratio distribution."""
    n: int = 0
    p05: float | None = None
    p25: float | None = None
    median: float | None = None
    p75: float | None = None
    p90: float | None = None
    p95: float | None = None
    p99: float | None = None
    mean: float | None = None
    n_zero: int = 0
    n_ge_5pct: int = 0
    n_ge_10pct: int = 0
    n_ge_50pct: int = 0


@dataclass
class Summary:
    n_introns: int = 0
    n_with_evidence: int = 0   # ≥1 contributing read (any signal type)

    # Per-boundary event totals
    sum_splice_left: int = 0
    sum_splice_right: int = 0
    sum_splice_exact: int = 0
    sum_retain_left: int = 0
    sum_retain_right: int = 0
    sum_crossing: int = 0

    # Depth-augmented totals
    sum_interior_reads: int = 0
    sum_intron_coverage_bp: int = 0
    sum_intron_length: int = 0  # only for introns counted in this summary

    # Global rates
    global_ir_rate_junction: float = 0.0  # retain / (retain + splice), event-weighted
    global_ir_rate_depth: float = 0.0     # coverage_bp / (coverage_bp + max_splice*intron_length), pooled

    # Per-intron distributions
    junction_stats: IRStats = field(default_factory=IRStats)
    depth_stats: IRStats = field(default_factory=IRStats)


def _quantile(sorted_xs: list[float], q: float) -> float | None:
    if not sorted_xs:
        return None
    pos = max(0, min(len(sorted_xs) - 1, int(q * len(sorted_xs))))
    return sorted_xs[pos]


def _fill_stats(stats: IRStats, values: list[float]) -> None:
    values.sort()
    stats.n = len(values)
    if not values:
        return
    stats.p05    = _quantile(values, 0.05)
    stats.p25    = _quantile(values, 0.25)
    stats.median = _quantile(values, 0.50)
    stats.p75    = _quantile(values, 0.75)
    stats.p90    = _quantile(values, 0.90)
    stats.p95    = _quantile(values, 0.95)
    stats.p99    = _quantile(values, 0.99)
    stats.mean   = sum(values) / len(values)
    stats.n_zero     = sum(1 for x in values if x == 0.0)
    stats.n_ge_5pct  = sum(1 for x in values if x >= 0.05)
    stats.n_ge_10pct = sum(1 for x in values if x >= 0.10)
    stats.n_ge_50pct = sum(1 for x in values if x >= 0.50)


_REQUIRED_COLS = {
    "chrom", "exon_overlap", "intron_length", "crossing_reads",
    "splice_left", "splice_right", "splice_exact",
    "retain_left", "retain_right",
    "interior_reads", "intron_coverage_bp", "intron_mean_depth",
    "ir_ratio_junction", "ir_ratio_depth",
}


def summarize(
    tsv_path: str,
    *,
    min_obs: int = 10,
    only_clean: bool = True,
    by_chrom: bool = False,
) -> dict:
    """Read a quant TSV and return a Summary, plus per-chrom Summaries if
    requested.

    `min_obs` thresholds inclusion in the per-intron distributions.  An intron
    enters the junction-IR distribution when `splice_left + splice_right +
    retain_left + retain_right >= min_obs`; it enters the depth-IR
    distribution when it has *any* of those events OR a non-zero
    `intron_coverage_bp`, normalised so that the equivalent
    `(coverage_bp / intron_length) + max_splice + retain_total >= min_obs`.
    """
    fh: IO[str]
    fh = gzip.open(tsv_path, "rt") if tsv_path.endswith(".gz") else open(tsv_path, "rt")
    try:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {n: i for i, n in enumerate(header)}
        missing = _REQUIRED_COLS - set(idx)
        if missing:
            raise ValueError(f"quant TSV missing columns: {sorted(missing)}")

        agg = Summary()
        ir_j: list[float] = []
        ir_d: list[float] = []

        per_chrom: dict[str, Summary] = {}
        per_chrom_j: dict[str, list[float]] = {}
        per_chrom_d: dict[str, list[float]] = {}

        for line in fh:
            f = line.rstrip("\n").split("\t")
            chrom = f[idx["chrom"]]
            exon_overlap = f[idx["exon_overlap"]] == "1"
            if only_clean and exon_overlap:
                continue

            intron_length = int(f[idx["intron_length"]])
            sL = int(f[idx["splice_left"]])
            sR = int(f[idx["splice_right"]])
            se = int(f[idx["splice_exact"]])
            rL = int(f[idx["retain_left"]])
            rR = int(f[idx["retain_right"]])
            cr = int(f[idx["crossing_reads"]])
            ir_int = int(f[idx["interior_reads"]])
            cov_bp = int(f[idx["intron_coverage_bp"]])
            mean_d = float(f[idx["intron_mean_depth"]])
            ir_j_s = f[idx["ir_ratio_junction"]]
            ir_d_s = f[idx["ir_ratio_depth"]]

            def _ingest(s: Summary, ir_j_list: list[float], ir_d_list: list[float]) -> None:
                s.n_introns += 1
                s.sum_splice_left  += sL
                s.sum_splice_right += sR
                s.sum_splice_exact += se
                s.sum_retain_left  += rL
                s.sum_retain_right += rR
                s.sum_crossing     += cr
                s.sum_interior_reads     += ir_int
                s.sum_intron_coverage_bp += cov_bp
                s.sum_intron_length      += intron_length
                if cr > 0:
                    s.n_with_evidence += 1
                if (sL + sR + rL + rR) >= min_obs and ir_j_s != ".":
                    ir_j_list.append(float(ir_j_s))
                # Depth ratio distribution: include an intron when there is
                # *any* evidence reaching the same min_obs threshold, where
                # interior reads each count as 1.
                if (sL + sR + rL + rR + ir_int) >= min_obs and ir_d_s != ".":
                    ir_d_list.append(float(ir_d_s))

            _ingest(agg, ir_j, ir_d)
            if by_chrom:
                c = per_chrom.get(chrom)
                if c is None:
                    c = Summary()
                    per_chrom[chrom] = c
                    per_chrom_j[chrom] = []
                    per_chrom_d[chrom] = []
                _ingest(c, per_chrom_j[chrom], per_chrom_d[chrom])

        _finalize(agg, ir_j, ir_d)
        out = {"global": asdict(agg)}
        if by_chrom:
            out["by_chrom"] = {}
            for chrom, c in sorted(per_chrom.items()):
                _finalize(c, per_chrom_j[chrom], per_chrom_d[chrom])
                out["by_chrom"][chrom] = asdict(c)
        return out
    finally:
        fh.close()


def _finalize(s: Summary, ir_j: list[float], ir_d: list[float]) -> None:
    # Junction global rate: retain / (retain + splice) over all events
    denom_j = s.sum_splice_left + s.sum_splice_right + s.sum_retain_left + s.sum_retain_right
    if denom_j > 0:
        s.global_ir_rate_junction = (s.sum_retain_left + s.sum_retain_right) / denom_j

    # Depth global rate: pool coverage bp vs splice-events * intron length.
    # This treats each intron's "splice support" as max_splice spread over
    # its length, which is how the per-intron formula would average out.
    # In the pooled aggregate we approximate with `(sum_splice_left +
    # sum_splice_right) / 2` to get a mean splice-event-rate-per-position
    # signal once divided by intron length.  Concretely:
    #   per-intron: depth = cov_bp / L,  rate = depth / (depth + maxS)
    # Pooled equivalent:
    #   sum(cov_bp) / (sum(cov_bp) + sum(maxS * L))
    # We approximate sum(maxS * L) by ((sL + sR) / 2) * mean_L summed = pool.
    # Cheaper and as defensible: just use the symmetric form
    #   sum(cov_bp) / (sum(cov_bp) + (sL + sR) * mean_L / 2)
    # but mean_L over the introns counted varies; using sum_intron_length /
    # n_introns to keep it explicit.
    if s.n_introns > 0 and s.sum_intron_length > 0:
        mean_L = s.sum_intron_length / s.n_introns
        avg_splice = (s.sum_splice_left + s.sum_splice_right) / 2
        denom_d = s.sum_intron_coverage_bp + avg_splice * mean_L
        if denom_d > 0:
            s.global_ir_rate_depth = s.sum_intron_coverage_bp / denom_d

    _fill_stats(s.junction_stats, ir_j)
    _fill_stats(s.depth_stats, ir_d)


def render_text(out: dict, *, min_obs: int) -> str:
    lines = []

    def fmt(x):
        return "." if x is None else f"{x:.4f}"

    def emit(label: str, s: dict) -> None:
        lines.append(f"== {label} ==")
        lines.append(f"  introns:                          {s['n_introns']:>12,}")
        lines.append(f"  ... with any contributing read:   {s['n_with_evidence']:>12,}")
        lines.append("")
        lines.append("  -- junction events --")
        lines.append(f"  splice_left  | splice_right       {s['sum_splice_left']:>10,} | {s['sum_splice_right']:>10,}")
        lines.append(f"  splice_exact (one N op spans)     {s['sum_splice_exact']:>12,}")
        lines.append(f"  retain_left  | retain_right       {s['sum_retain_left']:>10,} | {s['sum_retain_right']:>10,}")
        lines.append("")
        lines.append("  -- interior + depth signal --")
        lines.append(f"  interior reads                    {s['sum_interior_reads']:>12,}")
        lines.append(f"  intron_coverage_bp (M bp inside)  {s['sum_intron_coverage_bp']:>12,}")
        lines.append(f"  sum_intron_length (counted)       {s['sum_intron_length']:>12,}")
        lines.append("")
        lines.append("  -- global rates --")
        lines.append(f"  global IR rate, junction-only:    {s['global_ir_rate_junction']:>12.4%}")
        lines.append(f"  global IR rate, depth-augmented:  {s['global_ir_rate_depth']:>12.4%}")
        lines.append("")

        for name, key in [("junction", "junction_stats"), ("depth", "depth_stats")]:
            stats = s[key]
            lines.append(f"  -- per-intron ir_ratio_{name} (n={stats['n']:,}, >={min_obs} obs) --")
            for q in ("p05", "p25", "median", "p75", "p90", "p95", "p99", "mean"):
                lines.append(f"    {q:<7}                       {fmt(stats[q]):>12}")
            lines.append(f"    n with IR=0:                  {stats['n_zero']:>12,}")
            lines.append(f"    n with IR>=5%:                {stats['n_ge_5pct']:>12,}")
            lines.append(f"    n with IR>=10%:               {stats['n_ge_10pct']:>12,}")
            lines.append(f"    n with IR>=50%:               {stats['n_ge_50pct']:>12,}")
            lines.append("")

    emit("global", out["global"])
    for chrom, s in out.get("by_chrom", {}).items():
        emit(chrom, s)
    return "\n".join(lines)
