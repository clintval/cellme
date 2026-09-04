import pytest

from cellme.builds import GenomeBuild
from cellme.builds import contig_order
from cellme.builds import contigs_for


def test_canonical_members_have_two_builds() -> None:
    assert list(GenomeBuild) == [GenomeBuild.GRCh38, GenomeBuild.GRCh37]


def test_members_are_strings() -> None:
    assert GenomeBuild.GRCh38 == "GRCh38"
    assert isinstance(GenomeBuild.GRCh37, str)


@pytest.mark.parametrize(
    ("alias", "expected"),
    [
        ("hg38", GenomeBuild.GRCh38),
        ("HG38", GenomeBuild.GRCh38),
        ("hg19", GenomeBuild.GRCh37),
        ("Hg19", GenomeBuild.GRCh37),
    ],
)
def test_hg_aliases_resolve(alias: str, expected: GenomeBuild) -> None:
    assert GenomeBuild(alias) is expected


def test_unknown_build_raises() -> None:
    with pytest.raises(ValueError):
        GenomeBuild("hg17")


def test_ucsc_names() -> None:
    assert GenomeBuild.GRCh38.ucsc_name == "hg38"
    assert GenomeBuild.GRCh37.ucsc_name == "hg19"


def test_each_build_has_all_primary_contigs() -> None:
    for build in GenomeBuild:
        names = [name for name, _length in contigs_for(build)]
        assert names == [*(str(index) for index in range(1, 23)), "X", "Y", "MT"]


def test_builds_disagree_on_chromosome_lengths() -> None:
    grch37 = dict(contigs_for(GenomeBuild.GRCh37))
    grch38 = dict(contigs_for(GenomeBuild.GRCh38))
    assert grch37["1"] != grch38["1"]
    assert grch37["MT"] == grch38["MT"]


def test_contig_order_is_karyotypic() -> None:
    order = contig_order(GenomeBuild.GRCh38)
    assert order["1"] < order["2"] < order["22"] < order["X"] < order["MT"]
