from openspliceai.create_data.immune_splice import SpliceSite
from openspliceai.create_data.immune_training import (
    combine_data,
    datafile_positive_sites,
    find_sites_missing_from_baseline,
    splice_sites_to_data,
    split_splice_sites_by_gene,
)


def _site(sequence: str = "AAAAGTAAAA", gene: str = "Ighj1") -> SpliceSite:
    return SpliceSite(
        source="NCBI-explicit",
        accession="j1",
        gene=gene,
        site_type="donor",
        label=2,
        index=3,
        oriented_sequence=sequence,
        feature_key="donor_splice_site",
        seqid="chr12",
        genomic_position=100,
        strand="+",
    )


def test_datafile_positive_sites_recovers_existing_label():
    data = [
        ["gene1"],
        ["chr12"],
        ["+"],
        ["1"],
        ["9"],
        ["AAAAGTAAA"],
        ["000200000"],
    ]
    sites = datafile_positive_sites(data, context_radius=3)
    assert len(sites) == 1
    assert sites[0].site_type == "donor"
    assert sites[0].context(3) == "AAAAGTA"


def test_ncbi_site_already_in_baseline_is_not_added():
    query = _site(sequence="AAAAGTAAAA")
    baseline = [
        SpliceSite(
            source="baseline",
            accession="x",
            gene="different-feature-id",
            site_type="donor",
            label=2,
            index=3,
            oriented_sequence="AAAAGTAAAA",
            feature_key="datafile-label",
        )
    ]
    missing, coverage = find_sites_missing_from_baseline(
        [query],
        baseline,
        context_radius=3,
    )
    assert missing == []
    assert coverage[0].covered is True


def test_ncbi_site_missing_from_baseline_is_added():
    missing, coverage = find_sites_missing_from_baseline(
        [_site()],
        [],
        context_radius=3,
    )
    assert len(missing) == 1
    assert coverage[0].covered is False


def test_splice_sites_to_data_places_label_at_center():
    data = splice_sites_to_data([_site()])
    assert len(data[0]) == 1
    assert data[5][0] == "AAAAGTAAAA"
    assert data[6][0] == "000200000"
    assert data[0][0].startswith("NCBI-IMMUNE:Ighj1:donor")


def test_split_sites_keeps_same_gene_together():
    first = _site(gene="Ighj1")
    second = SpliceSite(
        source=first.source,
        accession="j2",
        gene="Ighj1",
        site_type=first.site_type,
        label=first.label,
        index=first.index,
        oriented_sequence="CCCCGTCCCC",
        feature_key=first.feature_key,
        seqid=first.seqid,
        genomic_position=200,
        strand=first.strand,
    )
    train, validation = split_splice_sites_by_gene(
        [first, second],
        validation_ratio=0.5,
    )
    assert not (train and validation)
    assert len(train) + len(validation) == 2


def test_combine_data_appends_records_fieldwise():
    first = [[f"a{i}"] for i in range(7)]
    second = [[f"b{i}"] for i in range(7)]
    combined = combine_data(first, second)
    assert combined == [[f"a{i}", f"b{i}"] for i in range(7)]
