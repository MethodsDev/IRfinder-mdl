# IRfinder-mdl

> Implemented by Anthropic's Claude Opus 4.7 under direction from
> the Methods Development Lab at the Broad Institute (June 2026).

A junction-anchored intron retention counter for spliced-read BAMs, with an
optional intronic-depth signal alongside.

Given a reference annotation (GTF) and a sorted+indexed BAM of long-read
spliced alignments, for every annotated intron the tool counts:

**Junction signal — what does each read say at the boundaries?**

1. **splice** the boundary — a CIGAR `N` block whose endpoint lies within
   `--jitter` bp of the annotated boundary, supported by at least `--anchor`
   matched bases on the *exonic* side of that junction;
2. **retain** the boundary — at least `--anchor` matched bases on both the
   *exonic* and *intronic* side of the boundary, with no `N` op crossing the
   anchor window.

The junction-only IR ratio at each boundary is
$$\mathrm{IR}_{\text{side}} = \frac{R_{\text{side}}}{R_{\text{side}} + S_{\text{side}}}$$
with the per-intron summary
$$\mathrm{IR}_{\text{junction}} = \frac{R_L + R_R}{R_L + R_R + S_L + S_R}.$$

**Depth signal — how much intronic sequence is actually unspliced?**

For each annotated intron the tool also tallies:
- `interior_reads` — reads whose alignment lies fully inside `[s, e)` with
  no `N` op intersecting the intron;
- `intron_coverage_bp` — total matched bp inside `[s, e)` contributed by any
  read whose CIGAR has no `N` op intersecting `[s, e)` (boundary-retain
  reads, interior reads, or reads spliced through unrelated introns elsewhere).

These give the IRFinder-S-style ratio
$$\mathrm{IR}_{\text{depth}} = \frac{D}{D + \max(S_L, S_R)}, \quad D = \frac{\text{intron\_coverage\_bp}}{\text{intron\_length}}.$$

Both ratios are `.` (undefined) when there is no splice evidence to compare
against — the IR ratio is a *relative* measure and a locus with no observed
splicing is uninformative either way.

**Why two ratios?** Long-read RNA-seq frequently produces reads that cross a
boundary non-canonically (alt-3'/5' splice sites near the annotation) — these
show up as boundary-retain reads with little intronic coverage. They look
like IR through the junction lens but are really alternative splicing. The
depth view discriminates: an intron with high `IR_junction` and low
`IR_depth` is likely alternative splicing near a boundary; an intron with
both high is genuine retention. The two views are complementary.

The junction model is the same as IRFinder-S's per-intron junction counters
and the depth model is an adapted form of its `IntronDepth / MaxSplice`
ratio. No CNN filter is applied; reference-mappability filters are left to
the caller.

## Install

Python ≥ 3.10. Runtime deps: `pysam`, `intervaltree`.

```bash
git clone git@github.com:MethodsDev/IRfinder-mdl.git
cd IRfinder-mdl
pip install -e .[test]

python -m pytest                          # 47 tests, ~0.1 s
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
| `crossing_reads` | unique reads that contributed at least one signal at this intron |
| `splice_left`  / `splice_right` / `splice_exact` | per-boundary splice counts; `splice_exact` is reads where one `N` op covers both boundaries |
| `retain_left`  / `retain_right` | per-boundary retention counts |
| `retain_both`  | reads that retain *both* boundaries — cleanest junction IR witness |
| `mixed`        | reads that splice one boundary and retain the other (alternative splice-site usage) |
| `interior_reads` | reads whose alignment lies fully inside the intron with no `N` op intersecting it |
| `intron_coverage_bp` | total matched bp inside the intron from all reads whose CIGAR carries no `N` intersecting the intron |
| `intron_mean_depth` | `intron_coverage_bp / intron_length` |
| `ir_ratio_left`  | `retain_left / (retain_left + splice_left)`; `.` when `splice_left == 0` |
| `ir_ratio_right` | analogous on the 3' side |
| `ir_ratio_junction` | `(R_L + R_R) / (R_L + R_R + S_L + S_R)`; `.` when both sides have no splice events |
| `ir_ratio_depth` | `intron_mean_depth / (intron_mean_depth + max(S_L, S_R))`; `.` when `max(S_L, S_R) == 0` |

`splice_left + splice_right` *over-counts* every `splice_exact` read once
(it contributes to both). That is intentional: each boundary is its own
junction, and the IR ratio treats each as a weighted observation.

All IR ratios return `.` when the splice denominator is zero. The retention
counts are still reported, so a locus with retention reads but no splice
reads is preserved in the raw signal — but its ratio is undefined, not 1.0.
This avoids inflating IR distributions with lowly-expressed loci where stray
interior reads happen to land in an intron.

### 3. Summarize

```bash
python -m irfinder_mdl summarize \
    --quant ir.tsv.gz \
    --by-chrom
```

Reports global + per-chromosome totals and the per-intron IR-ratio quantile
distribution for both `ir_ratio_junction` and `ir_ratio_depth`. Restricted by
default to clean introns (no exon overlap) and to introns with at least
`--min-obs` (default 10) contributing reads. Use `--json` for a machine-
readable dump, `--include-exon-overlap` to also count alt-isoform-prone
introns.


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
- **Interior reads vs the depth signal.** A read whose alignment lies fully
  inside an intron only contributes to `interior_reads` and
  `intron_coverage_bp` when the read's CIGAR carries no `N` op intersecting
  the intron. Any `N` op (including ones that use the intron as a splice
  donor or acceptor) disqualifies that read's intronic bp — we refuse to
  credit retention coverage that overlaps a splice event in the same
  transcript molecule.
- **Two ratios, two stories.** `ir_ratio_junction` measures whether reads
  cross the boundary as exon or as intron; `ir_ratio_depth` measures how
  much intronic sequence is unspliced. Wherever they disagree the
  disagreement is informative — typically alternative splice-site usage
  shows up as high junction IR with near-zero depth IR.
- **Undefined ratios are not zero.** When `splice_left + splice_right == 0`
  (no splice evidence at this locus) every IR ratio is emitted as `.`. A
  locus with retention reads but no splice reads is *not* an IR call -- we
  cannot tell whether the intron is genuinely retained or the gene is
  simply not being processed. Raw counts (`retain_*`, `interior_reads`,
  `intron_coverage_bp`) remain available for downstream analysis.


## Example run: SBX cDNA on GENCODE v47

Numbers below are from a 10 %-subsampled SBX cDNA BAM (5.6 GB, single-end,
median read length ~400 bp, minimap2 `-ax splice -ub -G 1250k`) against
GENCODE v47 on GRCh38, 16 threads on a 16-core AMD EPYC.

| Step | Wall time | Output |
| --- | --- | --- |
| `build-introns` | 28 s    | 516,940 unique introns; 198,297 (38.4 %) with no exon overlap |
| `quantify` (16 threads) | ~7 m | per-intron TSV, ~12 MB gzipped |
| `summarize`     | ~1 s    | text or JSON report |

### Observed IR signal

Clean introns (no exon overlap), `--min-obs` 10:

| | junction | depth |
| --- | ---: | ---: |
| introns scored                | 68,713 | 70,040 |
| **global IR rate**            | **0.46 %** | **0.12 %** |
| introns with IR = 0           | 49,960 (72.7 %) | 18,798 (26.8 %) |
| introns with IR ≥ 5 %         |  5,247 ( 7.6 %) |  5,114 ( 7.3 %) |
| introns with IR ≥ 10 %        |  3,402 ( 4.9 %) |  3,280 ( 4.7 %) |
| introns with IR ≥ 50 %        |    793 ( 1.2 %) |    672 ( 1.0 %) |
| per-intron IR p50 / p90 / p95 | 0.000 / 0.029 / 0.097 | 0.000 / 0.028 / 0.091 |

Total junction events on clean introns:

| Signal | Reads |
| --- | ---: |
| `splice_exact` (one N op spans the whole intron) | 31.1 M |
| `splice_left` / `splice_right`                   | 41.2 M / 41.4 M |
| `retain_left` / `retain_right`                   | 199 K / 186 K |
| `interior_reads` (fully inside intron, no N)     | 548 K |
| `intron_coverage_bp` (M bp inside intron, no N)  | 210 M |

Both ratios agree closely on the IR signal in the bulk of the distribution
(p90 / p95 within ~0.01) but the global *rate* differs: the junction view
credits boundary-crossing reads regardless of how much intronic sequence they
actually deposit, so it picks up alternative-splice-site reads as IR; the
depth view sees those reads as having near-zero intronic coverage and
discounts them. Wherever the two ratios disagree, the discrepancy itself is
the diagnostic.

These are consistent with the published picture for bulk human IR (~5–10 %
of introns with measurable retention; healthy-tissue global rates < 1 %).

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
