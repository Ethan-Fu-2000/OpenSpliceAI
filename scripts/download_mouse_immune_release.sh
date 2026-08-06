#!/usr/bin/env bash
# 下载、校验、合并并解压最新的小鼠 NCBI+IMGT 数据 Release。

set -euo pipefail

REPOSITORY="${REPOSITORY:-Ethan-Fu-2000/OpenSpliceAI}"
DESTINATION="${1:-data}"
TAG="${2:-}"

for command in gh sha256sum tar zstd; do
    if ! command -v "${command}" >/dev/null 2>&1; then
        echo "缺少命令：${command}" >&2
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
    echo "在 ${REPOSITORY} 中未找到 mouse-immune-data Release。" >&2
    exit 1
fi

mkdir -p "${DESTINATION}"
WORK_DIR="$(mktemp -d)"
trap 'rm -rf "${WORK_DIR}"' EXIT

echo "从 ${REPOSITORY} 下载 ${TAG}"
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

# 只有在完整分卷和总归档校验均通过后，才替换旧数据。
# 这样会同时清除旧版遗留的 gencode/ 目录。
rm -rf "${DESTINATION}/mouse_immune"
tar --use-compress-program=zstd \
    -xf "${WORK_DIR}/mouse-immune-data.tar.zst" \
    -C "${DESTINATION}"

echo "已解压 ${TAG} 到 ${DESTINATION}/mouse_immune"
echo "清单：${DESTINATION}/mouse_immune/manifest.json"
