"""将 IMGT 明确剪接位点严格定位到 NCBI GRCm39。"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Sequence, Tuple

from openspliceai.create_data.immune_splice import ImgtRecord, SpliceSite


def _build_aligner(genome_fasta: Path | str):
    try:
        import mappy as mp
    except ImportError as error:
        raise RuntimeError(
            "严格匹配 IMGT 到 GRCm39 需要 mappy；请执行 pip install -e ."
        ) from error

    genome_fasta = Path(genome_fasta)
    index_path = Path(str(genome_fasta) + ".mmi")
    if index_path.exists():
        aligner = mp.Aligner(str(index_path), preset="sr", best_n=10)
    else:
        aligner = mp.Aligner(
            str(genome_fasta),
            preset="sr",
            best_n=10,
            fn_idx_out=str(index_path),
        )
    if not aligner:
        raise RuntimeError(f"无法建立 GRCm39 索引：{genome_fasta}")
    return aligner


def _exact_hits(aligner, query: str) -> List[object]:
    hits = []
    for hit in aligner.map(query):
        q_st = int(getattr(hit, "q_st", -1))
        q_en = int(getattr(hit, "q_en", -1))
        nm = int(getattr(hit, "NM", -1))
        mlen = int(getattr(hit, "mlen", q_en - q_st))
        blen = int(getattr(hit, "blen", q_en - q_st))
        if (
            q_st == 0
            and q_en == len(query)
            and nm == 0
            and mlen == len(query)
            and blen == len(query)
        ):
            hits.append(hit)
    unique = {}
    for hit in hits:
        key = (
            str(getattr(hit, "ctg", "")),
            int(getattr(hit, "r_st", -1)),
            int(getattr(hit, "r_en", -1)),
            int(getattr(hit, "strand", 0)),
        )
        unique[key] = hit
    return list(unique.values())


def filter_imgt_records_by_grcm39(
    records: Sequence[ImgtRecord],
    genome_fasta: Path | str,
    *,
    context_radius: int = 30,
    require_unique: bool = True,
    aligner=None,
) -> Tuple[List[ImgtRecord], Dict[str, object]]:
    """只保留中心上下文在 GRCm39 全长零错配定位的 IMGT 位点。"""
    if context_radius < 10:
        raise ValueError("context_radius 至少为 10")
    if aligner is None:
        aligner = _build_aligner(genome_fasta)

    kept_records: List[ImgtRecord] = []
    details = []
    accepted = 0
    rejected_no_match = 0
    rejected_multiple = 0
    rejected_ambiguous = 0
    total = 0

    for record in records:
        kept_sites = []
        for site in record.sites:
            total += 1
            query = site.context(context_radius).upper()
            hits = []
            status = ""
            if "N" in query:
                rejected_ambiguous += 1
                status = "rejected_ambiguous_context"
            else:
                hits = _exact_hits(aligner, query)
                if not hits:
                    rejected_no_match += 1
                    status = "rejected_no_exact_grcm39_match"
                elif require_unique and len(hits) != 1:
                    rejected_multiple += 1
                    status = "rejected_multiple_exact_grcm39_matches"
                else:
                    accepted += 1
                    status = "accepted_exact_grcm39_match"
                    hit = hits[0]
                    hit_strand = int(getattr(hit, "strand", 1))
                    r_st = int(getattr(hit, "r_st", -1))
                    r_en = int(getattr(hit, "r_en", -1))
                    genomic_position = (
                        r_st + context_radius + 1
                        if hit_strand == 1
                        else r_en - context_radius
                    )
                    kept_sites.append(
                        SpliceSite(
                            source=site.source,
                            accession=site.accession,
                            gene=site.gene,
                            site_type=site.site_type,
                            label=site.label,
                            index=site.index,
                            oriented_sequence=site.oriented_sequence,
                            feature_key=site.feature_key,
                            functionality=site.functionality,
                            seqid=str(getattr(hit, "ctg", "")),
                            genomic_position=genomic_position,
                            strand="+" if hit_strand == 1 else "-",
                        )
                    )
            details.append(
                {
                    "accession": site.accession,
                    "gene": site.gene,
                    "site_type": site.site_type,
                    "motif": site.motif,
                    "status": status,
                    "exact_match_count": len(hits),
                }
            )

        if kept_sites:
            kept_records.append(
                ImgtRecord(
                    accession=record.accession,
                    species=record.species,
                    sequence=record.sequence,
                    gene=record.gene,
                    functionality=record.functionality,
                    strand=record.strand,
                    sites=kept_sites,
                )
            )

    return kept_records, {
        "policy": (
            "site-centered context must map full-length to GRCm39 with NM=0; "
            "unique location required by default"
        ),
        "context_radius": context_radius,
        "input_record_count": len(records),
        "input_site_count": total,
        "accepted_record_count": len(kept_records),
        "accepted_exact_unique_site_count": accepted,
        "rejected_ambiguous_context": rejected_ambiguous,
        "rejected_no_exact_grcm39_match": rejected_no_match,
        "rejected_multiple_exact_grcm39_matches": rejected_multiple,
        "sites": details,
    }


def remove_imgt_sites_on_chromosomes(
    records: Sequence[ImgtRecord],
    excluded_chromosomes: Sequence[str],
) -> Tuple[List[ImgtRecord], int]:
    """排除已定位到测试染色体的 IMGT 位点。"""
    excluded = set(str(value) for value in excluded_chromosomes)
    kept_records = []
    removed = 0
    for record in records:
        sites = []
        for site in record.sites:
            if site.seqid in excluded:
                removed += 1
            else:
                sites.append(site)
        if sites:
            kept_records.append(
                ImgtRecord(
                    accession=record.accession,
                    species=record.species,
                    sequence=record.sequence,
                    gene=record.gene,
                    functionality=record.functionality,
                    strand=record.strand,
                    sites=sites,
                )
            )
    return kept_records, removed
