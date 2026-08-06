from types import SimpleNamespace

from openspliceai.create_data.imgt_genome_match import (
    filter_imgt_records_by_grcm39,
)
from openspliceai.create_data.immune_splice import ImgtRecord, SpliceSite


class FakeAligner:
    def __init__(self, hits):
        self.hits = hits

    def map(self, query):
        return iter(self.hits)


def _record():
    sequence = "A" * 30 + "AGT" + "C" * 30
    site = SpliceSite(
        source="IMGT",
        accession="A1",
        gene="IGHJ1*01",
        site_type="donor",
        label=2,
        index=30,
        oriented_sequence=sequence,
        feature_key="DONOR-SPLICE",
        functionality="functional",
    )
    return ImgtRecord(
        accession="A1",
        species="Mus musculus",
        sequence=sequence,
        gene="IGHJ1*01",
        functionality="functional",
        strand="+",
        sites=[site],
    )


def test_imgt_exact_unique_grcm39_match_is_kept():
    query_length = 61
    hit = SimpleNamespace(
        q_st=0,
        q_en=query_length,
        NM=0,
        mlen=query_length,
        blen=query_length,
        ctg="NC_000078.7",
        r_st=100,
        r_en=161,
        strand=1,
    )
    kept, report = filter_imgt_records_by_grcm39(
        [_record()],
        "unused.fa",
        context_radius=30,
        aligner=FakeAligner([hit]),
    )
    assert len(kept) == 1
    assert kept[0].sites[0].seqid == "NC_000078.7"
    assert kept[0].sites[0].genomic_position == 131
    assert report["accepted_exact_unique_site_count"] == 1


def test_imgt_multiple_exact_matches_are_rejected():
    query_length = 61
    hits = [
        SimpleNamespace(
            q_st=0,
            q_en=query_length,
            NM=0,
            mlen=query_length,
            blen=query_length,
            ctg="NC_000078.7",
            r_st=start,
            r_en=start + query_length,
            strand=1,
        )
        for start in (100, 500)
    ]
    kept, report = filter_imgt_records_by_grcm39(
        [_record()],
        "unused.fa",
        context_radius=30,
        aligner=FakeAligner(hits),
    )
    assert kept == []
    assert report["rejected_multiple_exact_grcm39_matches"] == 1
