"""Client for the cBioPortal REST API and its CCLE mutation study."""

import re
from collections.abc import Iterable
from dataclasses import dataclass
from typing import Any

import requests
from rapidfuzz import fuzz
from rapidfuzz import process

CBIOPORTAL_API: str = "https://www.cbioportal.org/api"
"""Base URL of the public cBioPortal REST API."""

DEFAULT_STUDY: str = "ccle_broad_2019"
"""The Cancer Cell Line Encyclopedia study used as the mutation source."""

_REQUEST_TIMEOUT: float = 60.0
"""Seconds to wait on any single cBioPortal request before giving up."""

_MISSING_VALUES: frozenset[str] = frozenset({"", ".", "NA", "N/A", "NULL"})
"""String placeholders that cBioPortal uses to mean a value is absent."""

_SUGGESTION_LIMIT: int = 5
"""How many close cell-line names to offer when a query does not resolve."""

_SUGGESTION_SCORE_CUTOFF: float = 70.0
"""The minimum rapidfuzz similarity, out of 100, for a name to be suggested."""


class CellLineError(ValueError):
    """Base class for failures to resolve a cell line query to a sample."""


class CellLineNotFoundError(CellLineError):
    """Raised when a query matches no cell line in the study."""


class AmbiguousCellLineError(CellLineError):
    """Raised when a query matches more than one cell line in the study."""


@dataclass(frozen=True)
class Mutation:
    """
    A single somatic mutation call for one cell line in the CCLE study.

    The fields mirror the subset of the cBioPortal mutation record needed to
    build a VCF record. Insertions and deletions follow the MAF convention in
    which the absent allele is encoded as a single dash.
    """

    gene: str
    entrez_gene_id: int | None
    chromosome: str
    start_position: int
    end_position: int
    reference_allele: str
    variant_allele: str
    protein_change: str | None
    variant_class: str
    variant_type: str
    ncbi_build: str
    refseq_mrna_id: str | None
    protein_pos_start: int | None
    protein_pos_end: int | None


def _clean(value: object) -> str | None:
    """
    Normalize a cBioPortal string field, mapping placeholders to None.

    Args:
        value: The raw field value from the API response.

    Returns:
        The stripped string, or None when the value is missing or a placeholder.
    """
    if value is None:
        return None
    text = str(value).strip()
    if text.upper() in _MISSING_VALUES:
        return None
    return text


def normalize_cell_line(name: str) -> str:
    """
    Reduce a cell line name to an uppercase alphanumeric matching key.

    This lets a user pass ``MOLT-4``, ``MOLT4``, or ``molt 4`` and have them all
    compare equal to the leading token of a CCLE sample identifier.

    Args:
        name: A cell line name or CCLE sample identifier token.

    Returns:
        The name uppercased with every non-alphanumeric character removed.
    """
    return re.sub(r"[^0-9A-Za-z]", "", name).upper()


def cell_line_name(sample_id: str) -> str:
    """
    Return the leading cell line token of a CCLE sample identifier.

    CCLE sample identifiers join a cell line name to its tissue of origin with
    underscores, for example ``MOLT4_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE``.

    Args:
        sample_id: A CCLE sample identifier.

    Returns:
        The substring before the first underscore.
    """
    return sample_id.split("_", 1)[0]


def humanize_cell_line(token: str) -> str:
    """
    Render a CCLE cell-line token in its conventional hyphenated form.

    A token of a letter prefix followed by a trailing digit run is split with a
    hyphen, so ``MOLT4`` becomes ``MOLT-4``. This is a light heuristic for
    readability; tokens that do not fit that simple shape are returned unchanged.

    Args:
        token: The leading cell line token of a CCLE sample identifier.

    Returns:
        The token with a hyphen inserted before a trailing digit run when one
        applies, otherwise the token unchanged.
    """
    match = re.fullmatch(r"([A-Za-z]+)(\d+)", token)
    if match is None:
        return token
    return f"{match.group(1)}-{match.group(2)}"


def suggest_cell_lines(
    query: str,
    sample_ids: Iterable[str],
    *,
    limit: int = _SUGGESTION_LIMIT,
    score_cutoff: float = _SUGGESTION_SCORE_CUTOFF,
) -> list[str]:
    """
    Return the closest CCLE cell-line names to a query by fuzzy similarity.

    The query is compared against the unique leading cell-line token of every
    sample, on their normalized alphanumeric form, using a Levenshtein-based
    similarity. Only names scoring at or above the cutoff are kept, so a query
    with no near neighbors yields an empty list.

    Args:
        query: The unresolved cell line identifier supplied by the user.
        sample_ids: All sample identifiers available in the study.
        limit: The maximum number of names to return.
        score_cutoff: The minimum similarity, out of 100, to include a name.

    Returns:
        Up to ``limit`` human-friendly cell-line names, best match first.
    """
    display_names: list[str] = []
    normalized_names: list[str] = []
    seen: set[str] = set()
    for sample_id in sample_ids:
        token = cell_line_name(sample_id)
        key = normalize_cell_line(token)
        if key in seen:
            continue
        seen.add(key)
        display_names.append(token)
        normalized_names.append(key)

    matches = process.extract(
        normalize_cell_line(query),
        normalized_names,
        scorer=fuzz.ratio,
        limit=limit,
        score_cutoff=score_cutoff,
    )
    return [humanize_cell_line(display_names[index]) for _name, _score, index in matches]


def resolve_sample(query: str, sample_ids: Iterable[str]) -> str:
    """
    Resolve a cell line query to a single CCLE sample identifier.

    A query first tries to match a whole sample identifier, then falls back to
    matching the leading cell line token, both compared on their normalized
    alphanumeric form.

    Args:
        query: The cell line identifier supplied by the user.
        sample_ids: All sample identifiers available in the study.

    Returns:
        The single matching sample identifier.

    Raises:
        CellLineNotFoundError: If no sample matches the query.
        AmbiguousCellLineError: If more than one sample matches the query.
    """
    wanted = normalize_cell_line(query)
    identifiers = list(sample_ids)

    whole = [s for s in identifiers if normalize_cell_line(s) == wanted]
    if len(whole) == 1:
        return whole[0]

    by_token = [s for s in identifiers if normalize_cell_line(cell_line_name(s)) == wanted]
    if len(by_token) == 1:
        return by_token[0]
    if not by_token:
        suggestions = suggest_cell_lines(query, identifiers)
        if suggestions:
            hint = "Did you mean: " + ", ".join(suggestions) + "?"
        else:
            hint = "No close matches."
        raise CellLineNotFoundError(f"no CCLE cell line matched {query!r}. {hint}")
    raise AmbiguousCellLineError(
        f"Query {query!r} matched multiple cell lines: {', '.join(sorted(by_token))}."
    )


def _mutation_from_record(record: dict[str, Any]) -> Mutation:
    """
    Build a Mutation from one cBioPortal mutation record.

    Args:
        record: A decoded JSON object from the mutations fetch endpoint.

    Returns:
        The parsed Mutation.
    """
    gene = record.get("gene") or {}
    symbol = gene.get("hugoGeneSymbol") or record.get("hugoGeneSymbol") or ""
    entrez = record.get("entrezGeneId")
    if entrez is None:
        entrez = gene.get("entrezGeneId")
    return Mutation(
        gene=symbol,
        entrez_gene_id=int(entrez) if entrez is not None else None,
        chromosome=str(record["chr"]),
        start_position=int(record["startPosition"]),
        end_position=int(record["endPosition"]),
        reference_allele=str(record.get("referenceAllele") or ""),
        variant_allele=str(record.get("variantAllele") or ""),
        protein_change=_clean(record.get("proteinChange")),
        variant_class=str(record.get("mutationType") or ""),
        variant_type=str(record.get("variantType") or ""),
        ncbi_build=str(record.get("ncbiBuild") or ""),
        refseq_mrna_id=_clean(record.get("refseqMrnaId")),
        protein_pos_start=record.get("proteinPosStart"),
        protein_pos_end=record.get("proteinPosEnd"),
    )


def fetch_sample_ids(
    study: str = DEFAULT_STUDY,
    *,
    session: requests.Session | None = None,
    base_url: str = CBIOPORTAL_API,
) -> list[str]:
    """
    Fetch every sample identifier in a cBioPortal study.

    Args:
        study: The cBioPortal study identifier.
        session: An optional requests session to reuse a connection.
        base_url: The base URL of the cBioPortal API.

    Returns:
        The sample identifiers reported by the study.
    """
    http = session or requests
    response = http.get(f"{base_url}/studies/{study}/samples", timeout=_REQUEST_TIMEOUT)
    response.raise_for_status()
    return [str(sample["sampleId"]) for sample in response.json()]


def fetch_mutations(
    study: str,
    sample_id: str,
    *,
    session: requests.Session | None = None,
    base_url: str = CBIOPORTAL_API,
) -> list[Mutation]:
    """
    Fetch all mutations for a single sample from a study's mutation profile.

    Args:
        study: The cBioPortal study identifier.
        sample_id: The sample identifier to fetch mutations for.
        session: An optional requests session to reuse a connection.
        base_url: The base URL of the cBioPortal API.

    Returns:
        The parsed mutations for the sample.
    """
    http = session or requests
    profile = f"{study}_mutations"
    url = f"{base_url}/molecular-profiles/{profile}/mutations/fetch"
    response = http.post(
        url,
        params={"projection": "DETAILED"},
        json={"sampleIds": [sample_id]},
        timeout=_REQUEST_TIMEOUT,
    )
    response.raise_for_status()
    return [_mutation_from_record(record) for record in response.json()]
