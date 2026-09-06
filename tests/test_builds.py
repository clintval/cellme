import pytest

from cellme.builds import ContigStyle
from cellme.builds import GenomeBuild
from cellme.builds import contig_order
from cellme.builds import contigs_for
from cellme.builds import styled_contig


def test_canonical_members_have_two_builds() -> None:
    assert list(GenomeBuild) == [GenomeBuild.hg38, GenomeBuild.hg19]


def test_members_are_strings() -> None:
    assert GenomeBuild.hg38 == "hg38"
    assert isinstance(GenomeBuild.hg19, str)


@pytest.mark.parametrize(
    ("value", "expected"),
    [
        ("hg38", GenomeBuild.hg38),
        ("HG38", GenomeBuild.hg38),
        ("hg19", GenomeBuild.hg19),
        ("Hg19", GenomeBuild.hg19),
    ],
)
def test_canonical_names_resolve(value: str, expected: GenomeBuild) -> None:
    assert GenomeBuild(value) is expected


@pytest.mark.parametrize(
    ("alias", "expected"),
    [
        ("GRCh38", GenomeBuild.hg38),
        ("grch38", GenomeBuild.hg38),
        ("GRCh37", GenomeBuild.hg19),
        ("GRCH37", GenomeBuild.hg19),
    ],
)
def test_grch_aliases_resolve(alias: str, expected: GenomeBuild) -> None:
    assert GenomeBuild(alias) is expected


def test_grch_aliases_are_not_canonical_members() -> None:
    # The GRCh names resolve but are hidden: they never appear as CLI choices.
    assert "GRCh38" not in [build.name for build in GenomeBuild]


def test_unknown_build_raises() -> None:
    with pytest.raises(ValueError):
        GenomeBuild("hg17")


def test_ucsc_names() -> None:
    assert GenomeBuild.hg38.ucsc_name == "hg38"
    assert GenomeBuild.hg19.ucsc_name == "hg19"


def test_grch_names() -> None:
    assert GenomeBuild.hg38.grch_name == "GRCh38"
    assert GenomeBuild.hg19.grch_name == "GRCh37"


def test_each_build_has_all_primary_contigs() -> None:
    for build in GenomeBuild:
        names = [name for name, _length in contigs_for(build)]
        assert names == [*(str(index) for index in range(1, 23)), "X", "Y", "MT"]


def test_builds_disagree_on_chromosome_lengths() -> None:
    hg19 = dict(contigs_for(GenomeBuild.hg19))
    hg38 = dict(contigs_for(GenomeBuild.hg38))
    assert hg19["1"] != hg38["1"]
    assert hg19["MT"] == hg38["MT"]


def test_contig_order_is_karyotypic() -> None:
    order = contig_order(GenomeBuild.hg38)
    assert order["1"] < order["2"] < order["22"] < order["X"] < order["MT"]


@pytest.mark.parametrize(
    ("name", "expected"),
    [("1", "chr1"), ("22", "chr22"), ("X", "chrX"), ("Y", "chrY")],
)
def test_styled_contig_ucsc_prefixes_primary_contigs(name: str, expected: str) -> None:
    assert styled_contig(name, ContigStyle.ucsc) == expected


def test_styled_contig_ucsc_mitochondrion_is_chr_m() -> None:
    assert styled_contig("MT", ContigStyle.ucsc) == "chrM"


@pytest.mark.parametrize("name", ["1", "22", "X", "Y", "MT"])
def test_styled_contig_ensembl_is_identity(name: str) -> None:
    assert styled_contig(name, ContigStyle.ensembl) == name
