from pathlib import Path

from openspliceai.create_data.immune_splice import (
    ACCEPTOR_LABEL,
    DONOR_LABEL,
    find_missing_imgt_sites,
    imgt_records_to_data,
    read_imgt_splice_records,
    read_ncbi_immune_splice_sites,
)


def test_imgt_explicit_j_donor_and_c_acceptor(tmp_path: Path):
    flat = tmp_path / "imgt.dat"
    flat.write_text(
        """ID   TEST1
AC   TEST1;
DE   Mus musculus functional IGHJ1*01 and IGHM*01 sequence
OS   Mus musculus (house mouse).
FT   J-REGION        1..8
FT                   /allele="IGHJ1*01"
FT                   /functionality="functional"
FT   DONOR-SPLICE    8..10
FT   ACCEPTOR-SPLICE 16..20
FT   C-REGION        19..28
FT                   /allele="IGHM*01"
SQ   Sequence 30 BP;
     aaaaaaaagttttttcagccgggggggggg 30
//
""",
        encoding="utf-8",
    )
    records = read_imgt_splice_records(flat)
    assert len(records) == 1
    sites = records[0].sites
    assert [(site.site_type, site.index, site.label) for site in sites] == [
        ("donor", 7, DONOR_LABEL),
        ("acceptor", 18, ACCEPTOR_LABEL),
    ]
    assert sites[0].gene == "IGHJ1*01"
    assert sites[1].gene == "IGHM*01"
    assert sites[0].motif == "GT"
    assert sites[1].motif == "AG"


def test_ncbi_adjacent_exons_do_not_invent_j_to_c_pair(tmp_path: Path):
    fasta = tmp_path / "mouse.fa"
    fasta.write_text(
        ">chr1\n" + "A" * 9 + "AGT" + "A" * 17 + "AG" + "C" * 30 + "\n"
    )
    gff = tmp_path / "mouse.gff3"
    gff.write_text(
        """##gff-version 3
chr1\tRefSeq\tgene\t10\t20\t.\t+\t.\tID=gene-j;gene=Ighj1;gene_biotype=J_segment
chr1\tRefSeq\tJ_gene_segment\t10\t20\t.\t+\t.\tID=jseg;Parent=gene-j
chr1\tRefSeq\texon\t10\t20\t.\t+\t.\tID=jex;Parent=jseg
chr1\tRefSeq\tgene\t30\t60\t.\t+\t.\tID=gene-c;gene=Ighm;gene_biotype=C_region
chr1\tRefSeq\tC_gene_segment\t30\t60\t.\t+\t.\tID=cseg;Parent=gene-c
chr1\tRefSeq\texon\t30\t40\t.\t+\t.\tID=cex1;Parent=cseg
chr1\tRefSeq\texon\t50\t60\t.\t+\t.\tID=cex2;Parent=cseg
""",
        encoding="utf-8",
    )
    sites = read_ncbi_immune_splice_sites(gff, fasta, context_radius=5)
    assert [
        (site.gene, site.site_type, site.genomic_position)
        for site in sites
    ] == [
        ("Ighm", "donor", 40),
        ("Ighm", "acceptor", 50),
    ]


def test_ncbi_context_prevents_duplicate_imgt_site(tmp_path: Path):
    flat = tmp_path / "imgt.dat"
    flat.write_text(
        """ID   TEST2
AC   TEST2;
DE   Mus musculus functional IGHM*01 sequence
OS   Mus musculus.
FT   DONOR-SPLICE    6..8
FT                   /allele="IGHM*01"
FT                   /functionality="functional"
SQ   Sequence 13 BP;
     aaaaagtaaaaaa 13
//
""",
        encoding="utf-8",
    )
    records = read_imgt_splice_records(flat)
    site = records[0].sites[0]
    from openspliceai.create_data.immune_splice import SpliceSite

    ncbi = [
        SpliceSite(
            source="NCBI-inferred",
            accession="x",
            gene="Ighm",
            site_type="donor",
            label=2,
            index=site.index,
            oriented_sequence=site.oriented_sequence,
            feature_key="adjacent-exon",
        )
    ]
    missing, report = find_missing_imgt_sites(records, ncbi, context_radius=4)
    assert missing == []
    assert report["covered_by_ncbi"] == 1
    assert report["supplemented_from_imgt"] == 0


def test_imgt_records_convert_to_open_splice_data(tmp_path: Path):
    flat = tmp_path / "imgt.dat"
    flat.write_text(
        """ID   TEST3
AC   TEST3;
DE   Mus musculus functional IGHJ1*01 sequence
OS   Mus musculus.
FT   DONOR-SPLICE    4..6
FT                   /allele="IGHJ1*01"
FT                   /functionality="functional"
SQ   Sequence 9 BP;
     aaagtaaaa 9
//
""",
        encoding="utf-8",
    )
    records = read_imgt_splice_records(flat)
    data = imgt_records_to_data(records)
    assert data[0][0].startswith("IMGT:TEST3")
    assert data[5][0] == "AAAGTAAAA"
    assert data[6][0] == "000200000"


def test_gene_grouped_split_keeps_same_gene_together(tmp_path: Path):
    from openspliceai.create_data.immune_splice import (
        ImgtRecord,
        SpliceSite,
        split_imgt_records_by_gene,
    )

    sequence = "AAAGTAAAA"
    records = []
    for accession in ("A", "B"):
        site = SpliceSite(
            source="IMGT",
            accession=accession,
            gene="IGHJ1*01",
            site_type="donor",
            label=2,
            index=2,
            oriented_sequence=sequence,
            feature_key="DONOR-SPLICE",
        )
        records.append(
            ImgtRecord(
                accession,
                "Mus musculus",
                sequence,
                "IGHJ1*01",
                "functional",
                "+",
                [site],
            )
        )
    train, validation = split_imgt_records_by_gene(
        records,
        validation_ratio=0.5,
    )
    assert not (train and validation)
    assert len(train) + len(validation) == 2


def test_h5_merge_round_trip(tmp_path: Path):
    from openspliceai.create_data.immune_splice import (
        merge_data,
        read_datafile_h5,
        write_datafile_h5,
    )

    base = [["base"] for _ in range(7)]
    supplement = [["extra"] for _ in range(7)]
    path = tmp_path / "datafile_train.h5"
    write_datafile_h5(path, merge_data(base, supplement))
    loaded = read_datafile_h5(path)
    assert loaded == [["base", "extra"] for _ in range(7)]
