"""将 NCBI 免疫剪接位点补入 OpenSpliceAI 训练数据的辅助函数。

这个模块解决一个容易漏标的问题：

* 作者原始 ``create-data --biotype protein-coding`` 主要处理普通 mRNA；
* NCBI GFF3 中的 ``J_gene_segment`` / ``C_gene_segment`` 即使存在明确剪接位点，
  也不一定进入作者生成的 protein-coding 基线 HDF5；
* 因此不能只用 NCBI 注释去“排除重复 IMGT 位点”，还必须确认该 NCBI 位点
  是否已经真正存在于基线 HDF5 中。

本模块负责：
1. 从现有 datafile_*.h5 中恢复已有的逐碱基剪接正例；
2. 判断哪些 NCBI 免疫剪接位点仍未进入基线；
3. 将缺失的 NCBI 位点转换成 OpenSpliceAI datafile 记录；
4. 按基因稳定地拆分到 train / validation，永不加入 test。
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
import hashlib
from typing import Dict, Iterable, List, MutableMapping, Sequence, Tuple

from openspliceai.create_data.immune_splice import (
    ACCEPTOR_LABEL,
    DATASET_NAMES,
    DONOR_LABEL,
    SpliceSite,
    context_identity,
    deduplicate_sites,
    normalize_gene_name,
)


@dataclass(frozen=True)
class SiteCoverage:
    """一个待检查位点相对于已有训练标签的覆盖情况。"""

    site: SpliceSite
    covered: bool
    best_identity: float


def _context_around(sequence: str, index: int, radius: int) -> str:
    """返回以 ``index`` 为中心的定长序列；越界位置补 N。"""
    start = index - radius
    end = index + radius + 1
    left = "N" * max(0, -start)
    right = "N" * max(0, end - len(sequence))
    body = sequence[max(0, start):min(len(sequence), end)]
    result = left + body + right
    expected = 2 * radius + 1
    if len(result) < expected:
        result += "N" * (expected - len(result))
    return result


def datafile_positive_sites(
    data: Sequence[Sequence[str]],
    *,
    source: str = "OpenSpliceAI-baseline",
    context_radius: int = 30,
) -> List[SpliceSite]:
    """从 OpenSpliceAI datafile 结构中恢复 donor / acceptor 正例。

    ``data`` 必须按 ``NAME, CHROM, STRAND, TX_START, TX_END, SEQ, LABEL``
    排列。这里不依赖基因名完成去重；同一参考基因组中的同一真实位点，
    其方向统一后的中心序列应完全相同。
    """
    if len(data) != len(DATASET_NAMES):
        raise ValueError("OpenSpliceAI datafile 字段数量不正确")

    sites: List[SpliceSite] = []
    for name, chrom, strand, _tx_start, _tx_end, sequence, labels in zip(*data):
        sequence = str(sequence).upper()
        labels = str(labels)
        if len(sequence) != len(labels):
            raise ValueError(
                f"记录 {name} 的 SEQ 与 LABEL 长度不一致："
                f"{len(sequence)} != {len(labels)}"
            )
        for index, value in enumerate(labels):
            if value not in {"1", "2"}:
                continue
            site_type = "acceptor" if value == "1" else "donor"
            label = ACCEPTOR_LABEL if value == "1" else DONOR_LABEL
            context = _context_around(sequence, index, context_radius)
            sites.append(
                SpliceSite(
                    source=source,
                    accession=str(name),
                    gene=str(name),
                    site_type=site_type,
                    label=label,
                    index=context_radius,
                    oriented_sequence=context,
                    feature_key="datafile-label",
                    seqid=str(chrom),
                    strand=str(strand),
                )
            )
    return deduplicate_sites(sites, context_radius=min(20, context_radius))


def combine_data(*datasets: Sequence[Sequence[str]]) -> List[List[str]]:
    """按字段拼接多个 OpenSpliceAI datafile 数据结构。"""
    result: List[List[str]] = [[] for _ in DATASET_NAMES]
    for data in datasets:
        if len(data) != len(DATASET_NAMES):
            raise ValueError("OpenSpliceAI datafile 字段数量不正确")
        for index, values in enumerate(data):
            result[index].extend(str(value) for value in values)
    return result


def find_sites_missing_from_baseline(
    query_sites: Sequence[SpliceSite],
    baseline_sites: Sequence[SpliceSite],
    *,
    context_radius: int = 20,
    identity_threshold: float = 1.0,
) -> Tuple[List[SpliceSite], List[SiteCoverage]]:
    """找出尚未真正进入作者基线 HDF5 的 NCBI 免疫位点。

    对 NCBI 和作者基线而言，两者都来自同一参考基因组，所以默认要求
    中心上下文完全一致（identity=1.0）。不强制匹配 gene 名称，因为作者
    基线的 ``NAME`` 常是 GFF feature ID，而不是标准 IG/TR gene symbol。
    """
    if not 0 <= identity_threshold <= 1:
        raise ValueError("identity_threshold 必须位于 [0, 1]")

    existing_exact = {
        (site.site_type, site.context(context_radius))
        for site in baseline_sites
    }
    existing_by_type: MutableMapping[str, List[str]] = defaultdict(list)
    if identity_threshold < 1.0:
        for site in baseline_sites:
            existing_by_type[site.site_type].append(site.context(context_radius))

    missing: List[SpliceSite] = []
    coverage: List[SiteCoverage] = []
    for site in deduplicate_sites(query_sites, context_radius=context_radius):
        query_context = site.context(context_radius)
        exact = (site.site_type, query_context) in existing_exact
        if exact:
            best_identity = 1.0
        elif identity_threshold < 1.0:
            best_identity = max(
                (
                    context_identity(query_context, candidate)
                    for candidate in existing_by_type[site.site_type]
                ),
                default=0.0,
            )
        else:
            best_identity = 0.0
        covered = best_identity >= identity_threshold
        coverage.append(SiteCoverage(site, covered, best_identity))
        if not covered:
            missing.append(site)
    return missing, coverage


def _validation_bucket(
    gene: str,
    *,
    validation_ratio: float,
    seed: str,
) -> bool:
    if not 0 <= validation_ratio < 1:
        raise ValueError("validation_ratio 必须位于 [0, 1)")
    normalized = normalize_gene_name(gene) or "UNKNOWN"
    digest = hashlib.sha256(f"{seed}:{normalized}".encode("utf-8")).hexdigest()
    bucket = int(digest[:8], 16) % 10_000
    return bucket < int(validation_ratio * 10_000)


def split_splice_sites_by_gene(
    sites: Sequence[SpliceSite],
    *,
    validation_ratio: float = 0.1,
    seed: str = "openspliceai-imgt-v1",
) -> Tuple[List[SpliceSite], List[SpliceSite]]:
    """按基因稳定拆分 NCBI 免疫位点，避免同一基因横跨 train/validation。"""
    train: List[SpliceSite] = []
    validation: List[SpliceSite] = []
    for site in sites:
        target = validation if _validation_bucket(
            site.gene,
            validation_ratio=validation_ratio,
            seed=seed,
        ) else train
        target.append(site)
    return train, validation


def splice_sites_to_data(
    sites: Sequence[SpliceSite],
    *,
    name_prefix: str = "NCBI-IMMUNE",
) -> List[List[str]]:
    """把 site-centered NCBI 位点转换为 OpenSpliceAI datafile 记录。

    每个位点形成一条独立记录。``oriented_sequence`` 已经按转录方向统一，
    因此 LABEL 可以直接在 ``site.index`` 位置标 donor=2 / acceptor=1。
    """
    data: List[List[str]] = [[] for _ in DATASET_NAMES]
    for ordinal, site in enumerate(sites, start=1):
        sequence = site.oriented_sequence.upper()
        if not 0 <= site.index < len(sequence):
            raise ValueError(
                f"位点索引越界：{site.accession} index={site.index}, "
                f"sequence_length={len(sequence)}"
            )
        labels = ["0"] * len(sequence)
        labels[site.index] = str(site.label)
        position = site.genomic_position if site.genomic_position is not None else "NA"
        gene = site.gene or "unknown"
        accession = site.accession or f"site-{ordinal}"
        name = (
            f"{name_prefix}:{gene}:{site.site_type}:"
            f"{site.seqid or 'unknown'}:{position}:{accession}"
        )
        values = (
            name,
            site.seqid or name_prefix,
            site.strand or "+",
            "1",
            str(len(sequence)),
            sequence,
            "".join(labels),
        )
        for index, value in enumerate(values):
            data[index].append(str(value))
    return data


def coverage_to_report(coverage: Iterable[SiteCoverage]) -> List[Dict[str, object]]:
    """把 NCBI→基线覆盖结果转换为可写入 JSON 的明细。"""
    details: List[Dict[str, object]] = []
    for item in coverage:
        site = item.site
        details.append(
            {
                "source": site.source,
                "accession": site.accession,
                "gene": site.gene,
                "site_type": site.site_type,
                "seqid": site.seqid,
                "genomic_position": site.genomic_position,
                "strand": site.strand,
                "motif": site.motif,
                "already_in_baseline": item.covered,
                "best_baseline_context_identity": round(item.best_identity, 4),
            }
        )
    return details
