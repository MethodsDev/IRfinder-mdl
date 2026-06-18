version 1.0

#
# IRfinder-mdl: junction-anchored intron retention from a spliced-read BAM,
# with an optional intronic-depth signal.
#
# Workflow steps:
#   1. (optional) BuildIntrons -- derive a unique-intron TSV from a GTF.
#      Skipped when `introns_tsv` is supplied.
#   2. Quantify -- per-intron splice/retention/interior evidence from the BAM.
#   3. Summarize -- global + per-chromosome IR rates, text and JSON.
#
# Inputs are described in `parameter_meta` below.  Either `introns_tsv` (a
# pre-built unique-intron TSV) or `gtf` is required; if both are supplied,
# `introns_tsv` wins.
#

workflow IRfinderMDL {
    meta {
        author: "Methods Development Lab, The Broad Institute"
        email: "bhaas@broadinstitute.org"
        description: "Run IRfinder-mdl on a single BAM."
        version: "0.0.1"
    }

    parameter_meta {
        bam: "Sorted, indexed BAM of spliced read alignments (minimap2 -ax splice for long reads)."
        bam_index: "BAI index alongside the BAM."
        sample_id: "Output prefix for emitted files."
        introns_tsv: "Optional pre-built unique-intron TSV from `irfinder-mdl build-introns`. If supplied, the GTF is ignored."
        gtf: "Reference GTF (Ensembl/GENCODE).  Required when `introns_tsv` is not supplied."
        docker: "Container image carrying the irfinder-mdl CLI."
        anchor: "bp of matched alignment required on each side of every boundary check."
        jitter: "bp tolerance for an N op aligning to an annotated boundary."
        min_mapq: "Discard reads with MAPQ below this."
        exclude_flags: "SAM flag mask of reads to drop (2304 = 0x900 = secondary | supplementary)."
        skip_exon_overlap: "Drop introns whose interval overlaps any annotated exon before quantifying."
    }

    input {
        # Required: BAM (sorted, indexed)
        File bam
        File bam_index
        String sample_id = "sample"

        # Provide ONE of these:
        File? introns_tsv
        File? gtf

        # Container / runtime
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/irfinder-mdl/irfinder-mdl:latest"

        # Quantify knobs
        Int anchor = 8
        Int jitter = 3
        Int min_mapq = 1
        Int exclude_flags = 2304   # 0x900 = secondary | supplementary
        Boolean skip_exon_overlap = false

        # Resource overrides
        Int build_introns_cpu = 2
        Int build_introns_memory_gb = 8
        Int quantify_cpu = 16
        Int quantify_memory_gb = 16
        Int summarize_cpu = 1
        Int summarize_memory_gb = 4

        # Optional preemptible retries (GCP/AWS backends)
        Int preemptible = 0
    }

    # Build introns from GTF only if no pre-built TSV was supplied.
    # `select_first([gtf])` throws at workflow start when neither input was
    # provided, which gives a clearer error than a downstream task crash.
    if (!defined(introns_tsv)) {
        call BuildIntrons {
            input:
                gtf = select_first([gtf]),
                sample_id = sample_id,
                docker = docker,
                cpu = build_introns_cpu,
                memory_gb = build_introns_memory_gb,
                preemptible = preemptible,
        }
    }

    File introns_resolved = select_first([introns_tsv, BuildIntrons.introns])

    call Quantify {
        input:
            bam = bam,
            bam_index = bam_index,
            introns = introns_resolved,
            sample_id = sample_id,
            anchor = anchor,
            jitter = jitter,
            min_mapq = min_mapq,
            exclude_flags = exclude_flags,
            skip_exon_overlap = skip_exon_overlap,
            cpu = quantify_cpu,
            memory_gb = quantify_memory_gb,
            docker = docker,
            preemptible = preemptible,
    }

    call Summarize {
        input:
            quant_tsv = Quantify.quant_tsv,
            sample_id = sample_id,
            docker = docker,
            cpu = summarize_cpu,
            memory_gb = summarize_memory_gb,
            preemptible = preemptible,
    }

    output {
        File introns_table = introns_resolved
        File quant_tsv     = Quantify.quant_tsv
        File summary_text  = Summarize.summary_text
        File summary_json  = Summarize.summary_json
    }
}


task BuildIntrons {
    input {
        File gtf
        String sample_id
        String docker
        Int cpu
        Int memory_gb
        Int preemptible
    }

    # GTF is ~600 MB for GENCODE human; conservatively allocate 4x for
    # decompression / scratch.  Floor at 10 GB.
    Int disk_gb = ceil(size(gtf, "GB") * 4) + 10

    command <<<
        set -euo pipefail
        irfinder-mdl build-introns \
            --gtf ~{gtf} \
            --output ~{sample_id}.introns.tsv.gz
    >>>

    output {
        File introns = "~{sample_id}.introns.tsv.gz"
    }

    runtime {
        docker: docker
        cpu: cpu
        memory: "~{memory_gb} GB"
        disks: "local-disk ~{disk_gb} HDD"
        preemptible: preemptible
    }
}


task Quantify {
    input {
        File bam
        File bam_index
        File introns
        String sample_id
        Int anchor
        Int jitter
        Int min_mapq
        Int exclude_flags
        Boolean skip_exon_overlap
        Int cpu
        Int memory_gb
        String docker
        Int preemptible
    }

    # Need room for the BAM (input + symlink), the introns table, and the
    # output quant TSV.  The quant TSV is small (~10s of MB).
    Int disk_gb = ceil(size(bam, "GB") * 2 + size(introns, "GB") * 2) + 20

    command <<<
        set -euo pipefail

        # Co-locate the BAM and its index so pysam.fetch finds the .bai
        # automatically, regardless of where the runner localised each file.
        ln -s ~{bam}       input.bam
        ln -s ~{bam_index} input.bam.bai

        irfinder-mdl quantify \
            --bam input.bam \
            --introns ~{introns} \
            --output ~{sample_id}.ir.tsv.gz \
            --anchor ~{anchor} \
            --jitter ~{jitter} \
            --min-mapq ~{min_mapq} \
            --exclude-flags ~{exclude_flags} \
            ~{true="--skip-exon-overlap" false="" skip_exon_overlap} \
            --threads ~{cpu}
    >>>

    output {
        File quant_tsv = "~{sample_id}.ir.tsv.gz"
    }

    runtime {
        docker: docker
        cpu: cpu
        memory: "~{memory_gb} GB"
        disks: "local-disk ~{disk_gb} HDD"
        preemptible: preemptible
    }
}


task Summarize {
    input {
        File quant_tsv
        String sample_id
        String docker
        Int cpu
        Int memory_gb
        Int preemptible
    }

    # Quant TSV is ~12 MB gzipped; small headroom is enough.
    Int disk_gb = ceil(size(quant_tsv, "GB") * 4) + 5

    command <<<
        set -euo pipefail

        irfinder-mdl summarize \
            --quant ~{quant_tsv} \
            --by-chrom \
            > ~{sample_id}.ir_summary.txt

        irfinder-mdl summarize \
            --quant ~{quant_tsv} \
            --by-chrom \
            --json \
            > ~{sample_id}.ir_summary.json
    >>>

    output {
        File summary_text = "~{sample_id}.ir_summary.txt"
        File summary_json = "~{sample_id}.ir_summary.json"
    }

    runtime {
        docker: docker
        cpu: cpu
        memory: "~{memory_gb} GB"
        disks: "local-disk ~{disk_gb} HDD"
        preemptible: preemptible
    }
}
