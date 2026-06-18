"""Aggregate a per-intron quant TSV into genome- and chromosome-level metrics.

The headline numbers are:

* `global_ir_rate`         event-weighted, retained / (retained + spliced)
  summed over every (intron, boundary) pair that passes the filters;
* per-intron IR distribution -- quantiles of `ir_ratio` over filtered introns
  with sufficient junction support.
"""

from __future__ import annotations

import gzip
import json
import sys
from dataclasses import asdict, dataclass, field
from typing import IO


@dataclass
class Summary:
    n_introns: int = 0
    n_with_evidence: int = 0
    n_with_min_obs: int = 0

    sum_splice_left: int = 0
    sum_splice_right: int = 0
    sum_splice_exact: int = 0
    sum_retain_left: int = 0
    sum_retain_right: int = 0
    sum_crossing: int = 0

    global_ir_rate: float = 0.0
    ir_p05: float | None = None
    ir_p25: float | None = None
    ir_median: float | None = None
    ir_p75: float | None = None
    ir_p90: float | None = None
    ir_p95: float | None = None
    ir_p99: float | None = None
    ir_mean: float | None = None

    n_ir_zero: int = 0
    n_ir_ge_5pct: int = 0
    n_ir_ge_10pct: int = 0
    n_ir_ge_50pct: int = 0


def _quantile(sorted_xs: list[float], q: float) -> float | None:
    if not sorted_xs:
        return None
    pos = max(0, min(len(sorted_xs) - 1, int(q * len(sorted_xs))))
    return sorted_xs[pos]


def summarize(
    tsv_path: str,
    *,
    min_obs: int = 10,
    only_clean: bool = True,
    by_chrom: bool = False,
) -> dict:
    """Read a quant TSV and return a Summary (and optionally per-chrom Summaries)."""
    fh: IO[str]
    fh = gzip.open(tsv_path, "rt") if tsv_path.endswith(".gz") else open(tsv_path, "rt")
    try:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {n: i for i, n in enumerate(header)}
        required = {"chrom", "exon_overlap", "crossing_reads",
                    "splice_left", "splice_right", "splice_exact",
                    "retain_left", "retain_right", "ir_ratio"}
        missing = required - set(idx)
        if missing:
            raise ValueError(f"quant TSV missing columns: {sorted(missing)}")

        agg = Summary()
        per_chrom: dict[str, Summary] = {}
        per_chrom_ir: dict[str, list[float]] = {}
        ir_all: list[float] = []

        for line in fh:
            f = line.rstrip("\n").split("\t")
            chrom = f[idx["chrom"]]
            exon_overlap = f[idx["exon_overlap"]] == "1"
            if only_clean and exon_overlap:
                continue

            sL = int(f[idx["splice_left"]])
            sR = int(f[idx["splice_right"]])
            se = int(f[idx["splice_exact"]])
            rL = int(f[idx["retain_left"]])
            rR = int(f[idx["retain_right"]])
            cr = int(f[idx["crossing_reads"]])
            ir_s = f[idx["ir_ratio"]]

            agg.n_introns += 1
            agg.sum_splice_left  += sL
            agg.sum_splice_right += sR
            agg.sum_splice_exact += se
            agg.sum_retain_left  += rL
            agg.sum_retain_right += rR
            agg.sum_crossing     += cr
            if cr > 0:
                agg.n_with_evidence += 1
            if (sL + sR + rL + rR) >= min_obs and ir_s != ".":
                agg.n_with_min_obs += 1
                ir = float(ir_s)
                ir_all.append(ir)

            if by_chrom:
                c = per_chrom.get(chrom)
                if c is None:
                    c = Summary()
                    per_chrom[chrom] = c
                    per_chrom_ir[chrom] = []
                c.n_introns += 1
                c.sum_splice_left  += sL
                c.sum_splice_right += sR
                c.sum_splice_exact += se
                c.sum_retain_left  += rL
                c.sum_retain_right += rR
                c.sum_crossing     += cr
                if cr > 0:
                    c.n_with_evidence += 1
                if (sL + sR + rL + rR) >= min_obs and ir_s != ".":
                    c.n_with_min_obs += 1
                    per_chrom_ir[chrom].append(float(ir_s))

        _finalize(agg, ir_all)
        out = {"global": asdict(agg)}
        if by_chrom:
            out["by_chrom"] = {}
            for chrom, c in sorted(per_chrom.items()):
                _finalize(c, per_chrom_ir[chrom])
                out["by_chrom"][chrom] = asdict(c)
        return out
    finally:
        fh.close()


def _finalize(s: Summary, ir_values: list[float]) -> None:
    denom = s.sum_splice_left + s.sum_splice_right + s.sum_retain_left + s.sum_retain_right
    if denom > 0:
        s.global_ir_rate = (s.sum_retain_left + s.sum_retain_right) / denom
    ir_values.sort()
    s.ir_p05    = _quantile(ir_values, 0.05)
    s.ir_p25    = _quantile(ir_values, 0.25)
    s.ir_median = _quantile(ir_values, 0.50)
    s.ir_p75    = _quantile(ir_values, 0.75)
    s.ir_p90    = _quantile(ir_values, 0.90)
    s.ir_p95    = _quantile(ir_values, 0.95)
    s.ir_p99    = _quantile(ir_values, 0.99)
    s.ir_mean   = sum(ir_values) / len(ir_values) if ir_values else None
    s.n_ir_zero     = sum(1 for x in ir_values if x == 0.0)
    s.n_ir_ge_5pct  = sum(1 for x in ir_values if x >= 0.05)
    s.n_ir_ge_10pct = sum(1 for x in ir_values if x >= 0.10)
    s.n_ir_ge_50pct = sum(1 for x in ir_values if x >= 0.50)


def render_text(out: dict, *, min_obs: int) -> str:
    lines = []

    def emit(label: str, s: dict) -> None:
        lines.append(f"== {label} ==")
        lines.append(f"  introns:                       {s['n_introns']:>12,}")
        lines.append(f"  ... with crossing reads:       {s['n_with_evidence']:>12,}")
        lines.append(f"  ... with >={min_obs} junction obs:    {s['n_with_min_obs']:>12,}")
        lines.append(f"  splice_left:                   {s['sum_splice_left']:>12,}")
        lines.append(f"  splice_right:                  {s['sum_splice_right']:>12,}")
        lines.append(f"  splice_exact (single N op):    {s['sum_splice_exact']:>12,}")
        lines.append(f"  retain_left:                   {s['sum_retain_left']:>12,}")
        lines.append(f"  retain_right:                  {s['sum_retain_right']:>12,}")
        lines.append(f"  global IR rate:                {s['global_ir_rate']:>12.4%}")
        if s["ir_median"] is not None:
            lines.append("  per-intron IR ratio quantiles:")
            for q, key in [("p05","ir_p05"), ("p25","ir_p25"), ("med","ir_median"),
                           ("p75","ir_p75"), ("p90","ir_p90"), ("p95","ir_p95"),
                           ("p99","ir_p99")]:
                lines.append(f"    {q}:                          {s[key]:>12.4f}")
            lines.append(f"    mean:                         {s['ir_mean']:>12.4f}")
            lines.append(f"    n with IR=0:                  {s['n_ir_zero']:>12,}")
            lines.append(f"    n with IR>=5%:                {s['n_ir_ge_5pct']:>12,}")
            lines.append(f"    n with IR>=10%:               {s['n_ir_ge_10pct']:>12,}")
            lines.append(f"    n with IR>=50%:               {s['n_ir_ge_50pct']:>12,}")
        lines.append("")

    emit("global", out["global"])
    for chrom, s in out.get("by_chrom", {}).items():
        emit(chrom, s)
    return "\n".join(lines)
