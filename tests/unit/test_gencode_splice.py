from pathlib import Path

from openspliceai.create_data.gencode_splice import (
    read_gencode_basic_splice_sites,
)


def _write_fasta(path: Path, *, canonical: bool = True) -> None:
    sequence = list("A" * 80)
    if canonical:
        sequence[20:22] = list("GT")
        sequence[27:29] = list("AG")
    else:
        sequence[20:22] = list("TT")
        sequence[27:29] = list("CC")
    path.write_text(
        ">NC_000067.7 Mus musculus chromosome 1, GRCm39\n"
        + "".join(sequence)
        + "\n",
        encoding="utf-8",
    )


def _write_gtf(path: Path, *, level: int = 2) -> None:
    path.write_text(
        f"""##description: evidence-based annotation of the mouse genome (GRCm39), version M39
chr1\tHAVANA\ttranscript\t10\t40\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; gene_name "Gene1"; gene_type "protein_coding"; transcript_type "protein_coding"; level "{level}"; tag "basic";
chr1\tHAVANA\texon\t10\t20\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; gene_name "Gene1"; transcript_type "protein_coding"; level "{level}"; tag "basic";
chr1\tHAVANA\texon\t30\t40\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; gene_name "Gene1"; transcript_type "protein_coding"; level "{level}"; tag "basic";
""",
        encoding="utf-8",
    )


def test_gencode_m39_basic_canonical_pair(tmp_path: Path):
    fasta = tmp_path / "mouse.fa"
    gtf = tmp_path / "gencode.gtf"
    _write_fasta(fasta)
    _write_gtf(gtf)

    sites, report = read_gencode_basic_splice_sites(
        gtf,
        fasta,
        context_radius=5,
    )

    assert [(site.site_type, site.genomic_position) for site in sites] == [
        ("donor", 20),
        ("acceptor", 30),
    ]
    assert [site.motif for site in sites] == ["GT", "AG"]
    assert report["header_verified_release"] == "M39"
    assert report["header_verified_assembly"] == "GRCm39"
    assert report["canonical_intron_pairs"] == 1


def test_gencode_level3_is_rejected(tmp_path: Path):
    fasta = tmp_path / "mouse.fa"
    gtf = tmp_path / "gencode.gtf"
    _write_fasta(fasta)
    _write_gtf(gtf, level=3)

    sites, report = read_gencode_basic_splice_sites(
        gtf,
        fasta,
        context_radius=5,
    )

    assert sites == []
    assert report["excluded_by_level"] == 1


def test_gencode_noncanonical_pair_is_rejected(tmp_path: Path):
    fasta = tmp_path / "mouse.fa"
    gtf = tmp_path / "gencode.gtf"
    _write_fasta(fasta, canonical=False)
    _write_gtf(gtf)

    sites, report = read_gencode_basic_splice_sites(
        gtf,
        fasta,
        context_radius=5,
    )

    assert sites == []
    assert report["rejected_noncanonical_intron_pairs"] == 1
