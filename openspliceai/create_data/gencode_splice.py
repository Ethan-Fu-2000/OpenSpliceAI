"""读取 GENCODE Mouse Basic 注释并生成严格的小鼠剪接位点。

只接受与当前流程同一参考组装的 GENCODE Mouse M39（GRCm39）Basic GTF。
默认仅保留：

* GENCODE level 1 或 level 2；
* protein_coding 以及功能性 IG/TR 基因类型；
* 至少包含两个 exon 的 transcript；
* GT-AG、GC-AG 或 AT-AC 三类明确剪接对。

GENCODE 使用 ``chr1`` 等名称，NCBI GRCm39 FASTA 通常使用 RefSeq accession。
本模块从 NCBI FASTA 头部建立别名映射，因此不做 liftOver，也不修改坐标。
"""

from __future__ import annotations

from collections import defaultdict
from contextlib import contextmanager
from dataclasses import dataclass
import gzip
from pathlib import Path
import re
from typing import Dict, Iterable, Iterator, List, Mapping, MutableMapping, Sequence, TextIO, Tuple

from openspliceai.create_data.immune_splice import (
    IndexedFasta,
    SpliceSite,
    _extract_context,
    _site_from_context,
)


DEFAULT_GENCODE_RELEASE = "M39"
DEFAULT_GENCODE_ASSEMBLY = "GRCm39"
DEFAULT_TRANSCRIPT_TYPES = frozenset(
    {
        "protein_coding",
        "IG_C_gene",
        "IG_D_gene",
        "IG_J_gene",
        "IG_LV_gene",
        "IG_V_gene",
        "TR_C_gene",
        "TR_D_gene",
        "TR_J_gene",
        "TR_V_gene",
    }
)
CANONICAL_SPLICE_PAIRS = frozenset({("GT", "AG"), ("GC", "AG"), ("AT", "AC")})


@dataclass(frozen=True)
class GencodeTranscript:
    transcript_id: str
    gene: str
    transcript_type: str
    level: int
    seqid: str
    strand: str
    exons: Tuple[Tuple[int, int], ...]


@contextmanager
def open_text(path: Path | str) -> Iterator[TextIO]:
    """以文本方式打开普通文件或 gzip 文件。"""
    path = Path(path)
    if path.suffix == ".gz":
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
            yield handle
        return
    with path.open("rt", encoding="utf-8", errors="replace") as handle:
        yield handle


def parse_gtf_attributes(raw: str) -> Dict[str, List[str]]:
    """解析 GTF 第 9 列，保留重复的 tag 等字段。"""
    result: Dict[str, List[str]] = defaultdict(list)
    for item in raw.strip().strip(";").split(";"):
        item = item.strip()
        if not item:
            continue
        if " " not in item:
            continue
        key, value = item.split(None, 1)
        result[key].append(value.strip().strip('"'))
    return dict(result)


def _first(
    attrs: Mapping[str, Sequence[str]],
    *keys: str,
    default: str = "",
) -> str:
    for key in keys:
        values = attrs.get(key)
        if values:
            return str(values[0])
    return default


def _parse_level(attrs: Mapping[str, Sequence[str]]) -> int | None:
    value = _first(attrs, "level")
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _normalize_chromosome(value: str) -> str:
    value = value.strip()
    if value.lower().startswith("chr"):
        suffix = value[3:]
    else:
        suffix = value
    if suffix.upper() in {"MT", "M"}:
        return "chrM"
    if suffix.upper() in {"X", "Y"}:
        return f"chr{suffix.upper()}"
    if suffix.isdigit():
        return f"chr{int(suffix)}"
    return value


def ncbi_fasta_sequence_aliases(fasta_path: Path | str) -> Dict[str, str]:
    """从 NCBI FASTA 标题建立 ``chrN -> RefSeq accession`` 映射。"""
    aliases: Dict[str, str] = {}
    with Path(fasta_path).open("rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            header = line[1:].strip()
            accession = header.split(None, 1)[0]
            description = header[len(accession):].strip()

            chromosome = re.search(
                r"\bchromosome\s+([0-9]+|X|Y)\b",
                description,
                flags=re.IGNORECASE,
            )
            if chromosome:
                alias = _normalize_chromosome(chromosome.group(1))
                aliases[alias] = accession
                aliases[alias[3:]] = accession
                continue

            if re.search(r"\bmitochondr(?:ion|ial)\b", description, re.IGNORECASE):
                aliases["chrM"] = accession
                aliases["MT"] = accession
                aliases["M"] = accession

    return aliases


def _read_gencode_transcripts(
    gtf_path: Path | str,
    *,
    expected_release: str,
    expected_assembly: str,
    accepted_levels: Sequence[int],
    accepted_transcript_types: Sequence[str],
) -> Tuple[List[GencodeTranscript], Dict[str, object]]:
    metadata: Dict[str, Dict[str, object]] = {}
    exons: MutableMapping[str, List[Tuple[int, int]]] = defaultdict(list)
    header_lines: List[str] = []
    transcript_rows = 0
    exon_rows = 0

    with open_text(gtf_path) as handle:
        for raw_line in handle:
            if raw_line.startswith("#"):
                if len(header_lines) < 50:
                    header_lines.append(raw_line.strip())
                continue
            fields = raw_line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            seqid, _source, feature_type, start_s, end_s, _score, strand, _phase, attr_s = fields
            attrs = parse_gtf_attributes(attr_s)
            transcript_id = _first(attrs, "transcript_id")
            if not transcript_id:
                continue

            if feature_type == "transcript":
                transcript_rows += 1
                level = _parse_level(attrs)
                metadata[transcript_id] = {
                    "gene": _first(attrs, "gene_name", "gene_id", default=transcript_id),
                    "transcript_type": _first(
                        attrs,
                        "transcript_type",
                        "transcript_biotype",
                        "gene_type",
                    ),
                    "level": level,
                    "seqid": seqid,
                    "strand": strand,
                    "tags": tuple(attrs.get("tag", [])),
                }
                continue

            if feature_type != "exon":
                continue
            exon_rows += 1
            try:
                start, end = int(start_s), int(end_s)
            except ValueError:
                continue
            exons[transcript_id].append((start, end))
            if transcript_id not in metadata:
                metadata[transcript_id] = {
                    "gene": _first(attrs, "gene_name", "gene_id", default=transcript_id),
                    "transcript_type": _first(
                        attrs,
                        "transcript_type",
                        "transcript_biotype",
                        "gene_type",
                    ),
                    "level": _parse_level(attrs),
                    "seqid": seqid,
                    "strand": strand,
                    "tags": tuple(attrs.get("tag", [])),
                }

    header_text = "\n".join(header_lines)
    if expected_release and expected_release not in header_text:
        raise ValueError(
            f"GENCODE 注释头部未确认版本 {expected_release}；"
            "拒绝混用其他版本。"
        )
    if expected_assembly and expected_assembly not in header_text:
        raise ValueError(
            f"GENCODE 注释头部未确认组装 {expected_assembly}；"
            "拒绝进行坐标转换或跨版本混用。"
        )

    accepted_level_set = {int(value) for value in accepted_levels}
    accepted_type_set = set(accepted_transcript_types)
    transcripts: List[GencodeTranscript] = []
    excluded_level = 0
    excluded_type = 0
    excluded_single_exon = 0
    excluded_non_basic_tag = 0

    for transcript_id, info in metadata.items():
        level = info["level"]
        transcript_type = str(info["transcript_type"])
        transcript_exons = tuple(exons.get(transcript_id, ()))
        tags = tuple(str(value) for value in info.get("tags", ()))

        if level not in accepted_level_set:
            excluded_level += 1
            continue
        if transcript_type not in accepted_type_set:
            excluded_type += 1
            continue
        if len(transcript_exons) < 2:
            excluded_single_exon += 1
            continue
        if tags and "basic" not in tags:
            excluded_non_basic_tag += 1
            continue

        transcripts.append(
            GencodeTranscript(
                transcript_id=transcript_id,
                gene=str(info["gene"]),
                transcript_type=transcript_type,
                level=int(level),
                seqid=str(info["seqid"]),
                strand=str(info["strand"]),
                exons=transcript_exons,
            )
        )

    report: Dict[str, object] = {
        "header_verified_release": expected_release,
        "header_verified_assembly": expected_assembly,
        "transcript_rows": transcript_rows,
        "exon_rows": exon_rows,
        "accepted_transcripts": len(transcripts),
        "excluded_by_level": excluded_level,
        "excluded_by_transcript_type": excluded_type,
        "excluded_single_exon": excluded_single_exon,
        "excluded_non_basic_tag": excluded_non_basic_tag,
        "accepted_levels": sorted(accepted_level_set),
        "accepted_transcript_types": sorted(accepted_type_set),
    }
    return transcripts, report


def _deduplicate_by_context(
    sites: Iterable[SpliceSite],
    *,
    context_radius: int = 30,
) -> List[SpliceSite]:
    seen = set()
    result: List[SpliceSite] = []
    for site in sites:
        key = (site.site_type, site.context(context_radius))
        if key in seen:
            continue
        seen.add(key)
        result.append(site)
    return result


def read_gencode_basic_splice_sites(
    gtf_path: Path | str,
    genome_fasta: Path | str,
    *,
    context_radius: int = 1000,
    expected_release: str = DEFAULT_GENCODE_RELEASE,
    expected_assembly: str = DEFAULT_GENCODE_ASSEMBLY,
    accepted_levels: Sequence[int] = (1, 2),
    accepted_transcript_types: Sequence[str] = tuple(DEFAULT_TRANSCRIPT_TYPES),
) -> Tuple[List[SpliceSite], Dict[str, object]]:
    """读取 GENCODE Basic 中高置信度、规范 motif 的 exon 边界。"""
    if context_radius < 3:
        raise ValueError("context_radius 至少为 3")

    transcripts, report = _read_gencode_transcripts(
        gtf_path,
        expected_release=expected_release,
        expected_assembly=expected_assembly,
        accepted_levels=accepted_levels,
        accepted_transcript_types=accepted_transcript_types,
    )

    aliases = ncbi_fasta_sequence_aliases(genome_fasta)
    sites: List[SpliceSite] = []
    canonical_pairs = 0
    rejected_noncanonical = 0
    rejected_unmapped_seqid = 0
    rejected_invalid_strand = 0
    mapped_transcripts = 0

    with IndexedFasta(genome_fasta) as fasta:
        for transcript in transcripts:
            mapped_seqid = (
                transcript.seqid
                if transcript.seqid in fasta.entries
                else aliases.get(_normalize_chromosome(transcript.seqid), "")
            )
            if not mapped_seqid or mapped_seqid not in fasta.entries:
                rejected_unmapped_seqid += 1
                continue
            if transcript.strand not in {"+", "-"}:
                rejected_invalid_strand += 1
                continue

            mapped_transcripts += 1
            ordered = sorted(
                transcript.exons,
                key=lambda value: value[0],
                reverse=(transcript.strand == "-"),
            )
            for left, right in zip(ordered, ordered[1:]):
                if transcript.strand == "+":
                    donor_position = left[1]
                    acceptor_position = right[0]
                else:
                    donor_position = left[0]
                    acceptor_position = right[1]

                donor_context = _extract_context(
                    fasta,
                    mapped_seqid,
                    donor_position,
                    transcript.strand,
                    context_radius,
                )
                acceptor_context = _extract_context(
                    fasta,
                    mapped_seqid,
                    acceptor_position,
                    transcript.strand,
                    context_radius,
                )
                donor = _site_from_context(
                    source=f"GENCODE-{expected_release}-Basic",
                    accession=transcript.transcript_id,
                    gene=transcript.gene,
                    site_type="donor",
                    context=donor_context,
                    feature_key="basic-adjacent-exon",
                    seqid=mapped_seqid,
                    genomic_position=donor_position,
                    strand=transcript.strand,
                )
                acceptor = _site_from_context(
                    source=f"GENCODE-{expected_release}-Basic",
                    accession=transcript.transcript_id,
                    gene=transcript.gene,
                    site_type="acceptor",
                    context=acceptor_context,
                    feature_key="basic-adjacent-exon",
                    seqid=mapped_seqid,
                    genomic_position=acceptor_position,
                    strand=transcript.strand,
                )
                if (donor.motif, acceptor.motif) not in CANONICAL_SPLICE_PAIRS:
                    rejected_noncanonical += 1
                    continue
                canonical_pairs += 1
                sites.extend((donor, acceptor))

    unique_sites = _deduplicate_by_context(sites, context_radius=min(30, context_radius))
    report.update(
        {
            "mapped_transcripts": mapped_transcripts,
            "rejected_unmapped_seqid": rejected_unmapped_seqid,
            "rejected_invalid_strand": rejected_invalid_strand,
            "canonical_intron_pairs": canonical_pairs,
            "rejected_noncanonical_intron_pairs": rejected_noncanonical,
            "site_count_before_deduplication": len(sites),
            "unique_site_count": len(unique_sites),
            "coordinate_policy": (
                "GENCODE M39 GRCm39 coordinates mapped only by chromosome aliases "
                "to the NCBI GRCm39 FASTA; no liftOver"
            ),
        }
    )
    return unique_sites, report
