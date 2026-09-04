import pytest
from pytest import MonkeyPatch

from cellme.cbioportal import AmbiguousCellLineError
from cellme.cbioportal import CellLineNotFoundError
from cellme.cbioportal import cell_line_name
from cellme.cbioportal import fetch_mutations
from cellme.cbioportal import fetch_sample_ids
from cellme.cbioportal import humanize_cell_line
from cellme.cbioportal import normalize_cell_line
from cellme.cbioportal import resolve_sample
from cellme.cbioportal import suggest_cell_lines

SAMPLES = [
    "MOLT4_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE",
    "MOLT3_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE",
    "MOLT3_MOLT4",
    "JURKAT_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE",
]

PTEN_RECORD = {
    "gene": {"hugoGeneSymbol": "PTEN", "entrezGeneId": 5728},
    "entrezGeneId": 5728,
    "chr": "10",
    "startPosition": 89717770,
    "endPosition": 89717770,
    "referenceAllele": "A",
    "variantAllele": "-",
    "proteinChange": "K267Rfs*9",
    "mutationType": "Frame_Shift_Del",
    "variantType": "DEL",
    "ncbiBuild": "GRCh37",
    "refseqMrnaId": "NM_000314.4",
    "proteinPosStart": 265,
    "proteinPosEnd": 267,
}


class _FakeResponse:
    """A minimal stand-in for a requests Response carrying a canned payload."""

    def __init__(self, payload: object) -> None:
        self._payload = payload

    def raise_for_status(self) -> None:
        """Pretend the request succeeded."""

    def json(self) -> object:
        """Return the canned payload."""
        return self._payload


@pytest.mark.parametrize("query", ["MOLT-4", "molt4", "MOLT4", "  molt 4 "])
def test_resolve_matches_leading_token(query: str) -> None:
    assert resolve_sample(query, SAMPLES) == "MOLT4_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE"


def test_resolve_prefers_whole_sample_id() -> None:
    assert resolve_sample("MOLT3_MOLT4", SAMPLES) == "MOLT3_MOLT4"


def test_resolve_not_found_raises() -> None:
    with pytest.raises(CellLineNotFoundError):
        resolve_sample("HELA", SAMPLES)


def test_resolve_ambiguous_raises() -> None:
    with pytest.raises(AmbiguousCellLineError):
        resolve_sample("MOLT4", ["MOLT4_A", "MOLT4_B"])


def test_cell_line_name_takes_leading_token() -> None:
    assert cell_line_name("MOLT4_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE") == "MOLT4"


def test_normalize_strips_punctuation_and_uppercases() -> None:
    assert normalize_cell_line("MOLT-4") == "MOLT4"
    assert normalize_cell_line("nci h660") == "NCIH660"


def test_fetch_sample_ids(monkeypatch: MonkeyPatch) -> None:
    def fake_get(url: str, **_kwargs: object) -> _FakeResponse:
        assert url.endswith("/studies/ccle_broad_2019/samples")
        return _FakeResponse([{"sampleId": "MOLT4_X"}, {"sampleId": "JURKAT_Y"}])

    monkeypatch.setattr("cellme.cbioportal.requests.get", fake_get)
    assert fetch_sample_ids("ccle_broad_2019") == ["MOLT4_X", "JURKAT_Y"]


def test_fetch_mutations_parses_and_cleans(monkeypatch: MonkeyPatch) -> None:
    captured: dict[str, object] = {}

    def fake_post(url: str, **kwargs: object) -> _FakeResponse:
        captured["url"] = url
        captured["json"] = kwargs.get("json")
        return _FakeResponse([PTEN_RECORD, {**PTEN_RECORD, "refseqMrnaId": "NA"}])

    monkeypatch.setattr("cellme.cbioportal.requests.post", fake_post)
    mutations = fetch_mutations("ccle_broad_2019", "MOLT4_X")

    assert captured["json"] == {"sampleIds": ["MOLT4_X"]}
    assert str(captured["url"]).endswith("/ccle_broad_2019_mutations/mutations/fetch")
    assert mutations[0].gene == "PTEN"
    assert mutations[0].entrez_gene_id == 5728
    assert mutations[0].variant_allele == "-"
    assert mutations[0].refseq_mrna_id == "NM_000314.4"
    assert mutations[1].refseq_mrna_id is None


def test_humanize_cell_line_hyphenates_trailing_digits() -> None:
    assert humanize_cell_line("MOLT4") == "MOLT-4"
    assert humanize_cell_line("MOLT16") == "MOLT-16"
    assert humanize_cell_line("JURKAT") == "JURKAT"


def test_suggest_cell_lines_ranks_close_names() -> None:
    assert suggest_cell_lines("JURKKAT", SAMPLES) == ["JURKAT"]


def test_suggest_cell_lines_is_empty_below_cutoff() -> None:
    assert suggest_cell_lines("ZQXWVPYK", SAMPLES) == []


def test_resolve_typo_raises_with_suggestions() -> None:
    with pytest.raises(CellLineNotFoundError) as excinfo:
        resolve_sample("MOLT", SAMPLES)
    message = str(excinfo.value)
    assert "no CCLE cell line matched 'MOLT'" in message
    assert "Did you mean:" in message
    assert "MOLT-4" in message


def test_resolve_nonsense_reports_no_close_matches() -> None:
    with pytest.raises(CellLineNotFoundError) as excinfo:
        resolve_sample("ZQXWVPYK", SAMPLES)
    message = str(excinfo.value)
    assert "No close matches." in message
    assert "Did you mean" not in message


def test_resolve_exact_query_is_unchanged() -> None:
    assert resolve_sample("MOLT-4", SAMPLES) == "MOLT4_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE"


def test_fetch_mutations_maps_numeric_sex_chromosomes(monkeypatch: MonkeyPatch) -> None:
    def fake_post(_url: str, **_kwargs: object) -> _FakeResponse:
        return _FakeResponse([{**PTEN_RECORD, "chr": "23"}, {**PTEN_RECORD, "chr": "24"}])

    monkeypatch.setattr("cellme.cbioportal.requests.post", fake_post)
    mutations = fetch_mutations("ccle_broad_2019", "SAMPLE")
    assert mutations[0].chromosome == "X"
    assert mutations[1].chromosome == "Y"
