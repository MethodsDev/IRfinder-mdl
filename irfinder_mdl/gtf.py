"""GTF parsing and intron derivation.

This module reads a GTF (GENCODE / Ensembl style) and emits a table of
annotated introns suitable for quantification.  Every interval the rest of the
tool sees is 0-based half-open (BAM/BED convention), regardless of GTF's
1-based-inclusive native form.

We also flag introns whose interval overlaps an exon from *any* transcript at
the same locus — these are not "constitutive" introns and an unspliced read
through one of them may simply be the alternative isoform's exon, not retention
of the intron.  Downstream filtering by `exon_overlap == False` recovers the
clean intron set.
"""

from __future__ import annotations

import gzip
import re
import sys
from collections import defaultdict
from dataclasses import dataclass
from typing import IO, Iterable, Iterator

# Only attributes we care about.  Parsing is regex-based and stops at the first
# match so we don't pay for full GTF attribute parsing on every line.
_ATTR_RE = {
    "gene_id":       re.compile(r'gene_id "([^"]+)"'),
    "gene_name":     re.compile(r'gene_name "([^"]+)"'),
    "gene_type":     re.compile(r'gene_(?:type|biotype) "([^"]+)"'),
    "transcript_id": re.compile(r'transcript_id "([^"]+)"'),
    "transcript_type": re.compile(r'transcript_(?:type|biotype) "([^"]+)"'),
}


def _extract(attrs: str, key: str) -> str:
    m = _ATTR_RE[key].search(attrs)
    return m.group(1) if m else ""


@dataclass(frozen=True, slots=True)
class Exon:
    """0-based half-open."""
    chrom: str
    start: int
    end: int
    strand: str
    transcript_id: str
    gene_id: str
    gene_name: str
    gene_type: str
    transcript_type: str


@dataclass(slots=True)
class Intron:
    """0-based half-open.  Aggregated across transcripts sharing the interval."""
    chrom: str
    start: int
    end: int
    strand: str
    gene_ids: set[str]
    gene_names: set[str]
    gene_types: set[str]
    transcript_ids: set[str]
    transcript_types: set[str]
    exon_overlap: bool = False  # set in flag_exon_overlap()


def _open_maybe_gz(path: str) -> IO[str]:
    if path.endswith(".gz"):
        return gzip.open(path, "rt")  # type: ignore[return-value]
    return open(path, "rt")


def iter_exons(gtf_path: str) -> Iterator[Exon]:
    """Stream `exon` rows from a GTF as 0-based half-open Exon records."""
    with _open_maybe_gz(gtf_path) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "exon":
                continue
            chrom = fields[0]
            # GTF: 1-based inclusive [start, end]
            # BED: 0-based half-open [start-1, end)
            start = int(fields[3]) - 1
            end = int(fields[4])
            strand = fields[6]
            attrs = fields[8]
            yield Exon(
                chrom=chrom,
                start=start,
                end=end,
                strand=strand,
                transcript_id=_extract(attrs, "transcript_id"),
                gene_id=_extract(attrs, "gene_id"),
                gene_name=_extract(attrs, "gene_name"),
                gene_type=_extract(attrs, "gene_type"),
                transcript_type=_extract(attrs, "transcript_type"),
            )


def build_introns(exons: Iterable[Exon]) -> list[Intron]:
    """Derive introns from sorted-per-transcript exons.

    Two consecutive exons of the same transcript define one intron whose 0-based
    half-open coordinates are `(prev_exon.end, next_exon.start)`.  Introns
    sharing the same (chrom, start, end, strand) tuple across transcripts are
    merged so each unique splice-junction pair becomes one record.
    """
    by_tx: dict[str, list[Exon]] = defaultdict(list)
    for e in exons:
        by_tx[e.transcript_id].append(e)

    # key: (chrom, start, end, strand) -> Intron
    merged: dict[tuple[str, int, int, str], Intron] = {}

    for tx_id, tx_exons in by_tx.items():
        if len(tx_exons) < 2:
            continue
        tx_exons.sort(key=lambda x: x.start)
        for prev, nxt in zip(tx_exons, tx_exons[1:]):
            if prev.chrom != nxt.chrom or prev.strand != nxt.strand:
                continue  # shouldn't happen in a sane GTF, but be defensive
            i_start, i_end = prev.end, nxt.start
            if i_start >= i_end:
                continue  # touching or overlapping exons in this annotation
            key = (prev.chrom, i_start, i_end, prev.strand)
            rec = merged.get(key)
            if rec is None:
                rec = Intron(
                    chrom=prev.chrom,
                    start=i_start,
                    end=i_end,
                    strand=prev.strand,
                    gene_ids=set(),
                    gene_names=set(),
                    gene_types=set(),
                    transcript_ids=set(),
                    transcript_types=set(),
                )
                merged[key] = rec
            rec.gene_ids.add(prev.gene_id)
            if prev.gene_name:
                rec.gene_names.add(prev.gene_name)
            if prev.gene_type:
                rec.gene_types.add(prev.gene_type)
            rec.transcript_ids.add(tx_id)
            if prev.transcript_type:
                rec.transcript_types.add(prev.transcript_type)

    introns = list(merged.values())
    introns.sort(key=lambda i: (i.chrom, i.start, i.end, i.strand))
    return introns


def flag_exon_overlap(
    introns: list[Intron],
    exons: Iterable[Exon],
) -> None:
    """Mark introns whose interval overlaps *any* exon at the same locus.

    An intron with `exon_overlap = True` is part of an alternative-splicing
    region: a read that fails to splice through it could just be the isoform
    where that interval is an exon, not retention.  These should normally be
    excluded from IR analysis.
    """
    exons_by_chrom: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for e in exons:
        exons_by_chrom[e.chrom].append((e.start, e.end))

    # Per chrom: sort by start and precompute a prefix-max of ends so that
    # "does any exon with start < X have end > Y" becomes O(log n) via one
    # bisect plus one comparison.
    import bisect
    indexed: dict[str, tuple[list[int], list[int]]] = {}
    for chrom, ivs in exons_by_chrom.items():
        ivs.sort()  # by (start, end)
        starts = [a for a, _ in ivs]
        prefix_max_end = []
        running = -1
        for _, end in ivs:
            if end > running:
                running = end
            prefix_max_end.append(running)
        indexed[chrom] = (starts, prefix_max_end)

    for intron in introns:
        ce = indexed.get(intron.chrom)
        if ce is None:
            continue
        starts, prefix_max_end = ce
        # rightmost exon whose start < intron.end
        idx = bisect.bisect_left(starts, intron.end)
        if idx == 0:
            continue
        if prefix_max_end[idx - 1] > intron.start:
            intron.exon_overlap = True


# ---------------------------------------------------------------------------
# Serialization
# ---------------------------------------------------------------------------

INTRONS_TSV_HEADER = [
    "chrom", "start", "end", "strand", "intron_id",
    "gene_ids", "gene_names", "gene_types",
    "transcript_ids", "transcript_types",
    "intron_length", "exon_overlap",
]


def intron_id(intron: Intron) -> str:
    # Compact, parseable, position-keyed identifier (BED-style).
    return f"{intron.chrom}:{intron.start}-{intron.end}:{intron.strand}"


def write_introns_tsv(introns: list[Intron], path: str) -> None:
    """Emit a TSV with one row per unique annotated intron interval."""
    fh: IO[str]
    fh = gzip.open(path, "wt") if path.endswith(".gz") else open(path, "wt")
    try:
        fh.write("\t".join(INTRONS_TSV_HEADER) + "\n")
        for i in introns:
            fh.write("\t".join([
                i.chrom,
                str(i.start),
                str(i.end),
                i.strand,
                intron_id(i),
                ",".join(sorted(i.gene_ids)) if i.gene_ids else "",
                ",".join(sorted(i.gene_names)) if i.gene_names else "",
                ",".join(sorted(i.gene_types)) if i.gene_types else "",
                ",".join(sorted(i.transcript_ids)) if i.transcript_ids else "",
                ",".join(sorted(i.transcript_types)) if i.transcript_types else "",
                str(i.end - i.start),
                "1" if i.exon_overlap else "0",
            ]) + "\n")
    finally:
        fh.close()


def read_introns_tsv(path: str) -> list[Intron]:
    """Inverse of write_introns_tsv.  Returns sorted Intron records."""
    fh: IO[str]
    fh = gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")
    try:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        missing = set(INTRONS_TSV_HEADER) - set(idx)
        if missing:
            raise ValueError(f"introns TSV missing columns: {sorted(missing)}")
        out: list[Intron] = []
        for line in fh:
            f = line.rstrip("\n").split("\t")
            out.append(Intron(
                chrom=f[idx["chrom"]],
                start=int(f[idx["start"]]),
                end=int(f[idx["end"]]),
                strand=f[idx["strand"]],
                gene_ids=set(filter(None, f[idx["gene_ids"]].split(","))),
                gene_names=set(filter(None, f[idx["gene_names"]].split(","))),
                gene_types=set(filter(None, f[idx["gene_types"]].split(","))),
                transcript_ids=set(filter(None, f[idx["transcript_ids"]].split(","))),
                transcript_types=set(filter(None, f[idx["transcript_types"]].split(","))),
                exon_overlap=f[idx["exon_overlap"]] == "1",
            ))
        return out
    finally:
        fh.close()


def build_introns_from_gtf(gtf_path: str, *, log=print) -> list[Intron]:
    """End-to-end helper: GTF path -> sorted list of unique Introns with the
    `exon_overlap` flag populated.  Streams the GTF twice (once for exons used
    to build introns, once for the exon-overlap pass) to keep peak memory low
    on full human annotations."""
    log(f"[gtf] reading exons from {gtf_path}", file=sys.stderr)
    exons1 = list(iter_exons(gtf_path))
    log(f"[gtf] {len(exons1):,} exon rows", file=sys.stderr)
    introns = build_introns(exons1)
    log(f"[gtf] {len(introns):,} unique annotated introns", file=sys.stderr)
    # second pass: just for exon-overlap flag; could re-use exons1 but that
    # forces it to live for the whole call.  Re-streaming is cheap.
    flag_exon_overlap(introns, exons1)
    clean = sum(1 for i in introns if not i.exon_overlap)
    log(f"[gtf] {clean:,} introns with no exon overlap "
        f"({clean / max(1, len(introns)):.1%})", file=sys.stderr)
    return introns
