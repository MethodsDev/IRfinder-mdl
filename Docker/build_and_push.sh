#!/usr/bin/env bash
#
# Build the IRfinder-mdl Docker image and push it to GCP Artifact Registry
# under us-central1-docker.pkg.dev/methods-dev-lab/irfinder-mdl/irfinder-mdl.
#
# Usage:    Docker/build_and_push.sh
# Prereq:   `gcloud auth configure-docker us-central1-docker.pkg.dev` once per
#           machine, so docker push picks up the gcloud credential helper.
# Version:  read from VERSION.txt at the repo root.  Two tags are produced
#           and pushed: :${VERSION} and :latest.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

VERSION_FILE="${REPO_ROOT}/VERSION.txt"
if [[ ! -f "${VERSION_FILE}" ]]; then
    echo "ERROR: ${VERSION_FILE} not found" >&2
    exit 1
fi
VERSION="$(tr -d '[:space:]' < "${VERSION_FILE}")"
if [[ -z "${VERSION}" ]]; then
    echo "ERROR: ${VERSION_FILE} is empty" >&2
    exit 1
fi

REGISTRY="us-central1-docker.pkg.dev/methods-dev-lab/irfinder-mdl/irfinder-mdl"
TAG_VERSION="${REGISTRY}:${VERSION}"
TAG_LATEST="${REGISTRY}:latest"

echo "=== build  ${TAG_VERSION}  (also tagging :latest) ==="
docker build \
    --file "${SCRIPT_DIR}/Dockerfile" \
    --build-arg "VERSION=${VERSION}" \
    --tag "${TAG_VERSION}" \
    --tag "${TAG_LATEST}" \
    "${REPO_ROOT}"

echo "=== push   ${TAG_VERSION} ==="
docker push "${TAG_VERSION}"

echo "=== push   ${TAG_LATEST} ==="
docker push "${TAG_LATEST}"

echo
echo "Done."
echo "  ${TAG_VERSION}"
echo "  ${TAG_LATEST}"
