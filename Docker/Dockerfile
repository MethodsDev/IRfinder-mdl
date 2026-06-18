# syntax=docker/dockerfile:1.7
#
# IRfinder-mdl runtime image.
#
# The image installs the package via `pip install .` so the entrypoint is the
# regular `irfinder-mdl` console script.  `samtools` is bundled for convenience
# (BAM indexing, view, etc.); `pysam` ships its own htslib via wheels so we do
# not need libhts-dev or build tooling in the final image.
#
# Build context MUST be the repository root (the build_and_push.sh script
# enforces this).  `--build-arg VERSION=...` is informational only -- the
# package version is read from VERSION.txt at pip-install time via
# tool.setuptools.dynamic in pyproject.toml.

FROM python:3.12-slim

ARG VERSION=0+unknown

LABEL org.opencontainers.image.title="IRfinder-mdl"
LABEL org.opencontainers.image.description="Junction-anchored intron retention from spliced long-read BAM alignments."
LABEL org.opencontainers.image.source="https://github.com/MethodsDev/IRfinder-mdl"
LABEL org.opencontainers.image.licenses="MIT"
LABEL org.opencontainers.image.version="${VERSION}"
LABEL org.opencontainers.image.authors="Methods Development Lab, The Broad Institute"

# samtools is the only non-Python runtime dep we ship; everything else
# (htslib, libcurl, etc.) is statically embedded in the pysam wheel.
RUN apt-get update \
    && apt-get install -y --no-install-recommends samtools \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/irfinder-mdl
COPY pyproject.toml README.md VERSION.txt ./
COPY irfinder_mdl ./irfinder_mdl

RUN pip install --no-cache-dir .

# Sanity-check at build time so a broken image fails the build, not the user.
RUN irfinder-mdl --version

# Users should mount their data here.
WORKDIR /work

ENTRYPOINT ["irfinder-mdl"]
CMD ["--help"]
