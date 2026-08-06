#!/usr/bin/env bash
# Download, verify, join and extract the newest mouse immune data release.

set -euo pipefail

REPOSITORY="${REPOSITORY:-Ethan-Fu-2000/OpenSpliceAI}"
DESTINATION="${1:-data}"
TAG="${2:-}"

for command in gh sha256sum tar zstd; do
    if ! command -v "${command}" >/dev/null 2>&1; then
        echo "Missing required command: ${command}" >&2
        exit 1
    fi
done

if [[ -z "${TAG}" ]]; then
    TAG="$(
        gh release list \
            --repo "${REPOSITORY}" \
            --limit 100 \
            --json tagName,publishedAt \
            --jq '[.[] | select(.tagName | startswith("mouse-immune-data-"))]
                  | sort_by(.publishedAt)
                  | last
                  | .tagName'
    )"
fi

if [[ -z "${TAG}" || "${TAG}" == "null" ]]; then
    echo "No mouse-immune-data release was found in ${REPOSITORY}." >&2
    exit 1
fi

mkdir -p "${DESTINATION}"
WORK_DIR="$(mktemp -d)"
trap 'rm -rf "${WORK_DIR}"' EXIT

echo "Downloading release ${TAG} from ${REPOSITORY}"
gh release download "${TAG}" \
    --repo "${REPOSITORY}" \
    --dir "${WORK_DIR}" \
    --pattern 'mouse-immune-data.tar.zst.*.part' \
    --pattern 'PARTS_SHA256SUMS' \
    --pattern 'ARCHIVE_SHA256SUM' \
    --pattern 'manifest.json' \
    --clobber

(
    cd "${WORK_DIR}"
    sha256sum --check PARTS_SHA256SUMS
    cat mouse-immune-data.tar.zst.*.part > mouse-immune-data.tar.zst
    sha256sum --check ARCHIVE_SHA256SUM
)

tar --use-compress-program=zstd \
    -xf "${WORK_DIR}/mouse-immune-data.tar.zst" \
    -C "${DESTINATION}"

echo "Extracted ${TAG} to ${DESTINATION}/mouse_immune"
echo "Manifest: ${DESTINATION}/mouse_immune/manifest.json"
