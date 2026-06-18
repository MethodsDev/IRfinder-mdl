# IRfinder-mdl WDL workflow

A WDL 1.0 workflow that runs the full IRfinder-mdl pipeline on a single
sorted+indexed BAM. Compatible with Cromwell, miniWDL, and Terra.

## Workflow shape

```
                  ┌───────────────────────────┐
gtf ────(optional)─►  BuildIntrons             ├─► introns.tsv.gz
                  └───────────────────────────┘             │
                                                            ▼
introns_tsv (optional)─────────────────────────► ┌─────────────────┐
                                                 │   Quantify       ├─► ir.tsv.gz
bam, bam_index ───────────────────────────────► └─────────────────┘             │
                                                                                ▼
                                                            ┌───────────────────────┐
                                                            │   Summarize           ├─► .ir_summary.txt
                                                            └───────────────────────┘  .ir_summary.json
```

- **`BuildIntrons`** runs only when no pre-built `introns_tsv` is supplied;
  it derives the unique-intron table from `gtf`.
- **`Quantify`** runs the junction-anchored + depth-augmented per-intron
  counter on the BAM.
- **`Summarize`** emits both text and JSON reports.

Either `introns_tsv` or `gtf` must be provided. Supplying both ignores
`gtf`.

## Inputs

| Input | Type | Default | Description |
| --- | --- | --- | --- |
| `bam`        | `File` | — | Sorted, indexed BAM (minimap2 `-ax splice` for long reads). |
| `bam_index`  | `File` | — | BAI alongside the BAM. |
| `sample_id`  | `String` | `"sample"` | Output prefix. |
| `introns_tsv`| `File?` | unset | Pre-built unique-intron TSV from `irfinder-mdl build-introns`. |
| `gtf`        | `File?` | unset | Reference GTF (required when `introns_tsv` is unset). |
| `docker`     | `String` | `us-central1-docker.pkg.dev/methods-dev-lab/irfinder-mdl/irfinder-mdl:latest` | Container image to run. |
| `anchor`     | `Int` | `8` | bp of matched alignment required on each side of every boundary check. |
| `jitter`     | `Int` | `3` | bp tolerance for an `N` op aligning to an annotated boundary. |
| `min_mapq`   | `Int` | `1` | Discard reads with MAPQ below this. |
| `exclude_flags` | `Int` | `2304` | SAM flag mask of reads to drop (2304 = 0x900 = secondary \| supplementary). |
| `skip_exon_overlap` | `Boolean` | `false` | Drop alt-isoform-prone introns before quantifying. |
| `*_cpu`, `*_memory_gb` | `Int` | see WDL | Per-task resource overrides. |
| `preemptible` | `Int` | `0` | GCP/AWS preemptible retry count. |

## Outputs

| Output | Description |
| --- | --- |
| `introns_table` | The unique-intron TSV used by `Quantify`. |
| `quant_tsv`     | Per-intron counts and ratios. |
| `summary_text`  | Global + per-chromosome IR rates and per-intron quantile distributions, in human-readable text. |
| `summary_json`  | Same content as the text summary, as machine-readable JSON. |

## Quick start

### miniWDL (local, no Cromwell required)

```bash
miniwdl run WDL/irfinder_mdl.wdl \
    bam=path/to/sample.sorted.bam \
    bam_index=path/to/sample.sorted.bam.bai \
    introns_tsv=path/to/introns.tsv.gz \
    sample_id=SID004
```

Outputs land in a timestamped `<wd>/out/` directory; use `--dir` to control
where.

### Cromwell / Terra

Copy `WDL/inputs.template.json` to `inputs.json`, fill in the file paths,
then:

```bash
java -jar cromwell.jar run WDL/irfinder_mdl.wdl --inputs inputs.json
```

For Terra: import `WDL/irfinder_mdl.wdl` directly into a workspace, point
the inputs at workspace data, and submit.

## Notes

- The BAM and its `.bai` are co-located via symlink inside the `Quantify`
  task so `pysam.fetch()` can find the index regardless of where the
  runner localised each file.
- Disk sizing in each task is computed from input file sizes (`size(file,
  "GB")`) with a fixed overhead; override the resource defaults if you
  hit out-of-disk on very large BAMs.
- `Summarize` is run twice in the same task (once for text, once for
  JSON); the work is trivial.
