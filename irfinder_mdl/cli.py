"""IRfinder-mdl command line.

Subcommands:
  build-introns   GTF -> unique-intron TSV with exon-overlap flag
  quantify        introns TSV + BAM -> per-intron IR counts and ratios
  summarize       per-intron TSV -> global / per-chromosome IR report
"""

from __future__ import annotations

import argparse
import os
import sys

from .gtf import build_introns_from_gtf, read_introns_tsv, write_introns_tsv
from .quantify import QuantParams, quantify_introns, write_quant_tsv
from .version import __version__


def _apply_sample_prefix(output: str, sample_id: str | None) -> str:
    """Prefix the basename of `output` with ``<sample_id>.`` when `sample_id`
    is set.

    No-op when `sample_id` is falsy.  Idempotent: if the basename already
    starts with that prefix it is returned unchanged, so passing both an
    already-prefixed ``--output`` and ``--sample-id`` does not double up.
    Only the basename is touched; any directory component of `output` is
    preserved.
    """
    if not sample_id:
        return output
    d, b = os.path.split(output)
    prefix = f"{sample_id}."
    if b.startswith(prefix):
        return output
    return os.path.join(d, prefix + b)


def _add_build_introns(sub):
    p = sub.add_parser(
        "build-introns",
        help="Derive a unique-intron table from a GTF.",
        description="Stream a GTF and emit one row per unique annotated intron "
                    "interval (0-based half-open), with the set of supporting "
                    "transcripts/genes and a flag marking introns whose interval "
                    "overlaps an exon of any other annotated transcript at the "
                    "same locus.",
    )
    p.add_argument("--gtf", required=True, help="Reference GTF (may be .gz)")
    p.add_argument("--output", "-o", required=True,
                   help="Output introns TSV (.tsv or .tsv.gz)")
    p.add_argument("--sample-id", "-s", default=None,
                   help="If set, prefix the output basename with '<sample-id>.'")
    return p


def _add_quantify(sub):
    p = sub.add_parser(
        "quantify",
        help="Count junction-anchored splice/retention evidence per intron.",
        description="For each annotated intron, count reads from a sorted+indexed "
                    "BAM that splice or retain each of its 5' and 3' boundaries, "
                    "and report per-intron IR ratios.",
    )
    p.add_argument("--bam", required=True, help="Sorted+indexed BAM")
    p.add_argument("--introns", required=True,
                   help="Introns TSV from `build-introns`")
    p.add_argument("--output", "-o", required=True,
                   help="Output per-intron quant TSV (.tsv or .tsv.gz)")
    p.add_argument("--sample-id", "-s", default=None,
                   help="If set, prefix the output basename with '<sample-id>.'")
    p.add_argument("--anchor", type=int, default=8,
                   help="Min matched bp on each anchor side (default: 8)")
    p.add_argument("--jitter", type=int, default=3,
                   help="bp tolerance for splice junction position (default: 3)")
    p.add_argument("--min-mapq", type=int, default=1,
                   help="Discard reads with MAPQ below this (default: 1)")
    p.add_argument("--exclude-flags", type=lambda s: int(s, 0), default=0x900,
                   help="SAM flag mask of reads to exclude "
                        "(default: 0x900 = secondary|supplementary)")
    p.add_argument("--threads", "-t", type=int, default=max(1, os.cpu_count() or 1),
                   help="Parallel chromosome workers (default: all CPUs)")
    p.add_argument("--chrom", action="append", default=None,
                   help="Restrict to this chromosome (may be repeated; "
                        "useful for smoke tests). Takes precedence over the "
                        "primary-chromosome default.")
    p.add_argument("--all-chroms", action="store_true",
                   help="Quantify every chromosome including the mitochondrion "
                        "and unplaced/alt contigs (chrM/MT, GL*, KI*, chrUn_*, "
                        "*_random, *_alt). By default only nuclear primary "
                        "chromosomes (chr1-N, X, Y, with or without the chr "
                        "prefix) are examined.")
    p.add_argument("--skip-exon-overlap", action="store_true",
                   help="Skip introns flagged as overlapping an annotated exon. "
                        "These are alternative-isoform-prone and usually excluded "
                        "from IR analysis.")
    return p

def _add_summarize(sub):
    p = sub.add_parser(
        "summarize",
        help="Aggregate a per-intron quant TSV into a summary report.",
        description="Print or save a JSON / text summary of intron retention "
                    "across the genome and (optionally) per chromosome.",
    )
    p.add_argument("--quant", required=True, help="Output of `quantify`")
    p.add_argument("--min-obs", type=int, default=10,
                   help="Minimum splice+retain count for an intron to "
                        "contribute to the per-intron IR distribution "
                        "(default: 10)")
    p.add_argument("--include-exon-overlap", action="store_true",
                   help="Also include introns whose interval overlaps an "
                        "annotated exon (default: clean introns only)")
    p.add_argument("--by-chrom", action="store_true",
                   help="Also report a per-chromosome breakdown")
    p.add_argument("--json", action="store_true",
                   help="Emit JSON to stdout instead of a text table")
    return p


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="irfinder-mdl",
        description="Junction-anchored intron retention from spliced BAM alignments.",
    )
    parser.add_argument("--version", action="version", version=f"irfinder-mdl {__version__}")
    sub = parser.add_subparsers(dest="cmd", required=True)
    _add_build_introns(sub)
    _add_quantify(sub)
    _add_summarize(sub)
    args = parser.parse_args(argv)

    if args.cmd == "build-introns":
        introns = build_introns_from_gtf(args.gtf)
        output = _apply_sample_prefix(args.output, args.sample_id)
        write_introns_tsv(introns, output)
        print(f"[build-introns] wrote {len(introns):,} introns -> {output}",
              file=sys.stderr)
        return 0

    if args.cmd == "quantify":
        if not os.path.exists(args.bam + ".bai") and not os.path.exists(args.bam.removesuffix(".bam") + ".bai"):
            print(f"[quantify] WARNING: no index found alongside {args.bam}; "
                  "pysam.fetch() will fail. Run `samtools index` first.",
                  file=sys.stderr)
        params = QuantParams(
            anchor=args.anchor,
            jitter=args.jitter,
            min_mapq=args.min_mapq,
            exclude_flags=args.exclude_flags,
        )
        introns = read_introns_tsv(args.introns)
        if args.skip_exon_overlap:
            before = len(introns)
            introns = [i for i in introns if not i.exon_overlap]
            print(f"[quantify] dropped {before - len(introns):,} exon-overlap "
                  f"introns; {len(introns):,} remain", file=sys.stderr)
        # Primary-chromosome default: drop unplaced/alt contigs unless the user
        # asked for everything (--all-chroms) or named explicit chromosomes
        # (--chrom, which already restricts and takes precedence).
        if not args.all_chroms and not args.chrom:
            from .quantify import is_primary_chrom
            before = len(introns)
            introns = [i for i in introns if is_primary_chrom(i.chrom)]
            dropped = before - len(introns)
            if dropped:
                print(f"[quantify] primary-chromosome filter: dropped {dropped:,} "
                      f"introns on unplaced/alt contigs ({len(introns):,} remain); "
                      f"pass --all-chroms to include them", file=sys.stderr)
        counts = quantify_introns(
            args.bam, introns, params,
            threads=args.threads,
            chroms=args.chrom,
        )
        output = _apply_sample_prefix(args.output, args.sample_id)
        write_quant_tsv(introns, counts, output)
        n_evidence = sum(1 for c in counts if c.crossing_reads > 0)
        print(f"[quantify] {n_evidence:,}/{len(introns):,} introns had crossing "
              f"reads; wrote {output}", file=sys.stderr)
        return 0

    if args.cmd == "summarize":
        from .summarize import render_text, summarize
        out = summarize(
            args.quant,
            min_obs=args.min_obs,
            only_clean=not args.include_exon_overlap,
            by_chrom=args.by_chrom,
        )
        if args.json:
            import json as _json
            print(_json.dumps(out, indent=2))
        else:
            print(render_text(out, min_obs=args.min_obs))
        return 0

    parser.error("unknown command")
    return 2


if __name__ == "__main__":
    sys.exit(main())
