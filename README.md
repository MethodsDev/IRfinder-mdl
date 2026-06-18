# IRfinder-mdl

> Implemented by Anthropic's Claude Opus 4.7 under direction from
> the Methods Development Lab at the Broad Institute (June 2026).

A minimal, junction-anchored intron-retention counter for spliced-read BAMs.

Given a reference annotation (GTF) and a sorted+indexed BAM of long-read
spliced alignments, for every annotated intron the tool counts, at each of the
two boundaries (5'/donor and 3'/acceptor), reads that:

1. **splice** the boundary — a CIGAR `N` block whose endpoint lies within
   `--jitter` bp of the boundary, supported by at least `--anchor` matched
   bases on the *exonic* side of that junction; or
2. **retain** the boundary — at least `--anchor` matched bases on both the
   *exonic* and *intronic* side of the boundary, with no `N` op crossing the
   anchor window.

The IR ratio at each boundary is then
$$\mathrm{IR}_{\text{side}} = \frac{R_{\text{side}}}{R_{\text{side}} + S_{\text{side}}}$$
and the per-intron summary
$$\mathrm{IR} = \frac{R_L + R_R}{R_L + R_R + S_L + S_R}.$$

Reads that lie entirely within an intron contribute **no** evidence — only
boundary-crossing reads count. This is a deliberate departure from
coverage-based tools like iREAD; for long reads with truncated 5' ends, fully
intronic reads are frequently 3'-fragments of the gene and not retention
witnesses.

The model is the same as IRFinder's per-intron junction counters; the
implementation is in Python, takes any minimap2-spliced BAM, and adds no
intron-coverage signal or CNN filter.

## Install

Python ≥ 3.10. Runtime deps: `pysam`, `intervaltree`.

```bash
git clone git@github.com:MethodsDev/IRfinder-mdl.git
cd IRfinder-mdl
pip install -e .[test]

python -m pytest                          # 34 tests, ~0.1 s
irfinder-mdl --help                       # or: python -m irfinder_mdl --help
```

## Usage

### 1. Build the intron table from a GTF

```bash
python -m irfinder_mdl build-introns \
    --gtf path/to/annotations.gtf \
    --output introns.tsv.gz
```

Output columns:

| Column | Description |
| --- | --- |
| `chrom`, `start`, `end`, `strand` | 0-based half-open BED coordinates |
| `intron_id` | `chrom:start-end:strand` |
| `gene_ids` / `gene_names` / `gene_types` | comma-separated, union over supporting transcripts |
| `transcript_ids` / `transcript_types` | comma-separated, union of every transcript whose splice-junction pair defines this intron |
| `intron_length` | `end - start` |
| `exon_overlap` | `1` if any annotated exon (any strand, any transcript) overlaps `[start, end)`. These introns share a region with an alternative-exon isoform and a non-splicing read through one of them is biologically ambiguous; filter with `--skip-exon-overlap` for clean IR calls. |

### 2. Quantify

```bash
python -m irfinder_mdl quantify \
    --bam path/to/alignments.sorted.bam \
    --introns introns.tsv.gz \
    --output ir.tsv.gz \
    --skip-exon-overlap \
    --threads 16
```

Knobs (all default to IRFinder-S long-read settings):

| Flag | Default | Meaning |
| --- | --- | --- |
| `--anchor`        | `8`     | bp of matched alignment required on each side of every boundary check |
| `--jitter`        | `3`     | bp tolerance for an `N` op aligning to an annotated boundary |
| `--min-mapq`      | `1`     | discard reads with MAPQ below this (excludes multi-mappers at `0`) |
| `--exclude-flags` | `0x900` | SAM flag mask of reads to drop (secondary \| supplementary) |
| `--threads`       | all CPUs| parallel chromosome workers |
| `--chrom`         | —       | restrict to one chromosome (repeatable). Useful for smoke tests. |

### Output columns

| Column | Meaning |
| --- | --- |
| `chrom`, `start`, `end`, `strand`, `intron_id`, `gene_ids`, … | mirror the introns table |
| `intron_length`, `exon_overlap` | mirror the introns table |
| `crossing_reads` | unique reads that triggered ≥1 of the boolean signals below |
| `splice_left`    | reads with an `N` op anchored at the 5' boundary |
| `splice_right`   | reads with an `N` op anchored at the 3' boundary |
| `splice_exact`   | reads where one `N` op covers the full intron (both boundaries) |
| `retain_left`    | reads with matched alignment continuously through the 5' boundary |
| `retain_right`   | reads with matched alignment continuously through the 3' boundary |
| `retain_both`    | reads that retain both boundaries (cleanest IR witness) |
| `mixed`          | reads that splice one boundary and retain the other (alternative splice-site usage at this intron) |
| `ir_ratio_left`  | `retain_left / (retain_left + splice_left)`; `.` if denominator is 0 |
| `ir_ratio_right` | analogous |
| `ir_ratio`       | `(retain_L + retain_R) / (retain_L + retain_R + splice_L + splice_R)` |

`splice_left + splice_right` *over-counts* every `splice_exact` read once
(it contributes to both). That is intentional: each boundary is its own
junction, and the IR ratio treats each as a weighted observation.

## Coordinate convention

GTF is 1-based inclusive on both ends. Internally the tool converts everything
to 0-based half-open (BAM/BED convention). All output coordinates are
0-based half-open. To get a 1-based inclusive coordinate back: `start + 1`,
`end` is unchanged.

## Caveats and design notes

- **Strand.** Read strand is not used. minimap2 `-ax splice -ub` produces
  unstranded alignments, and the `XS` tag (intron strand from splice motif)
  isn't reliable for retention reads (no splice motif present). All read
  evidence is collapsed by genomic interval.
- **MAPQ.** Default `--min-mapq 1` drops MAPQ-0 multi-mappers, which is
  important for paralogous gene families and repetitive UTRs.
- **Exon overlap.** Without `--skip-exon-overlap`, introns whose interval is
  exonic in some other transcript are included with a flag set. Filter
  downstream as needed; for any "rate of intron retention" headline number,
  use the `exon_overlap == 0` subset.
- **Long introns.** If the BAM was aligned with a minimap2 `-G` cap (e.g.
  `1.25 Mb` in this dataset), introns larger than `-G` may show artificially
  high retention because minimap2 couldn't open the splice. Filter on
  `intron_length` if needed.
- **CIGAR semantics.** `M`, `=`, and `X` count as "matched bases". `D`
  advances the reference but creates a gap in matched intervals; a `D` inside
  an anchor window costs you anchor coverage. `I`, `S`, `H`, `P` don't
  advance the reference.
- **Spliced-exact accounting.** A single `N` op whose two endpoints both pass
  the jitter+anchor test at this intron's annotated boundaries increments
  `splice_exact`. Multiple smaller `N` ops crossing the intron region do not
  set `splice_exact` even if each end approximately matches an annotated
  boundary; the intent is to identify reads that report this exact intron as
  one contiguous splice event.

### 3. Summarize

```bash
python -m irfinder_mdl summarize \
    --quant output/ir_genome.tsv.gz \
    --by-chrom
```

Reports global + per-chromosome totals and the per-intron IR-ratio quantile
distribution, restricted to introns with no exon overlap and at least
`--min-obs` (default 10) junction observations.  Use `--json` for a machine-
readable dump, `--include-exon-overlap` to also count alt-isoform-prone
introns.

## Example run: SBX cDNA on GENCODE v47

Numbers below are from a 10 %-subsampled SBX cDNA BAM (5.6 GB, single-end,
median read length ~400 bp, minimap2 `-ax splice -ub -G 1250k`) against
GENCODE v47 on GRCh38, 16 threads on a 16-core AMD EPYC.

| Step | Wall time | Output |
| --- | --- | --- |
| `build-introns` | 28 s    | 516,940 unique introns; 198,297 (38.4 %) with no exon overlap |
| `quantify` (16 threads) | 6 m 36 s | per-intron TSV, 11 MB gzipped |
| `summarize`     | ~1 s    | text or JSON report |

### Observed IR signal

Clean introns (no exon overlap), ≥ 10 junction observations, N = 69,436:

| Statistic | Value |
| --- | --- |
| global IR rate (event-weighted) | **0.46 %** |
| introns with IR = 0    | 49,960 (72.0 %) |
| introns with IR ≥ 5 %  |  5,970 ( 8.6 %) |
| introns with IR ≥ 10 % |  4,125 ( 5.9 %) |
| introns with IR ≥ 50 % |  1,516 ( 2.2 %) |
| per-intron IR p50 / p90 / p95 | 0.000 / 0.036 / 0.135 |
| `splice_exact` reads (one N op spans the whole intron) | 31.1 M |
| `splice_left` + `splice_right` (per-boundary splice events) | 41.2 M / 41.4 M |
| `retain_left` + `retain_right` (per-boundary retention events) | 199 K / 186 K |

These match the published picture for bulk human IR (~5–10 % of introns with
measurable retention; healthy-tissue global event rates well under 1 %).

## Quick start (chr22 only)

```bash
# 1. derive the intron table from a GTF (~30 s on GENCODE v47)
python -m irfinder_mdl build-introns \
    --gtf <REFERENCE>.gtf \
    --output introns.tsv.gz

# 2. quantify on a single chromosome to sanity-check
python -m irfinder_mdl quantify \
    --bam <ALIGNMENTS>.bam \
    --introns introns.tsv.gz \
    --output ir_chr22.tsv.gz \
    --chrom chr22 --threads 1

# 3. report
python -m irfinder_mdl summarize \
    --quant ir_chr22.tsv.gz
```

Drop `--chrom` and bump `--threads` for the full genome.

## Authorship

Designed and implemented by **Claude Opus 4.7** (Anthropic), working from a
specification provided by the Methods Development Lab at the Broad Institute.
The algorithm — junction-anchored splice-vs-retention counting on annotated
introns, intentionally ignoring fully-intronic reads — was the human ask;
the code, tests, and documentation were authored by the model end-to-end in a
single design session.

If you use this tool, please cite the upstream IRFinder-S paper for the
per-intron junction-counting model:

> Lorenzi C, Barriere S, Arnold K, et al. **IRFinder-S: a comprehensive suite
> to discover and explore intron retention.** *Genome Biology* 22, 307 (2021).
> [doi:10.1186/s13059-021-02515-8](https://doi.org/10.1186/s13059-021-02515-8)
