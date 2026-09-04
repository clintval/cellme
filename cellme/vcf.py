"""Conversion of CCLE mutations into a sorted, well-described truth-track VCF."""

from collections.abc import Callable
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

import pysam
from liftover import get_lifter

from cellme.builds import GenomeBuild
from cellme.builds import contig_order
from cellme.builds import contigs_for
from cellme.cbioportal import Mutation

LiftPosition = Callable[[str, int], int | None]
"""A callable mapping a source (contig, 1-based position) to a lifted position."""

AnchorBase = Callable[[str, int], str]
"""A callable returning the single reference base at a (contig, 1-based position)."""

_DASH: str = "-"
"""The MAF sentinel for the absent allele of an insertion or deletion."""

InfoValue = str | int | bool
"""The value types cellme writes into a VCF INFO field."""


@dataclass(frozen=True)
class InfoField:
    """A single INFO field definition used for both the header and records."""

    key: str
    number: str
    type: str
    description: str


INFO_FIELDS: tuple[InfoField, ...] = (
    InfoField("GENE", "1", "String", "HUGO gene symbol"),
    InfoField(
        "PROTEIN_CHANGE",
        "1",
        "String",
        "Protein-level change in CCLE short form (HGVS p.-like), e.g. R306*",
    ),
    InfoField(
        "VARIANT_CLASS",
        "1",
        "String",
        "Variant classification (MAF-style), e.g. Missense_Mutation",
    ),
    InfoField(
        "VARIANT_TYPE",
        "1",
        "String",
        "Sequence alteration type reported by CCLE: SNP, DNP, INS, or DEL",
    ),
    InfoField("CELL_LINE", "1", "String", "Cell line resolved from the query"),
    InfoField("SAMPLE_ID", "1", "String", "CCLE sample identifier for the cell line"),
    InfoField("ENTREZ", "1", "Integer", "NCBI Entrez gene identifier"),
    InfoField("REFSEQ", "1", "String", "RefSeq mRNA accession for the annotated transcript"),
    InfoField("PROTEIN_POS", "1", "String", "Affected protein position or range, 1-based"),
    InfoField(
        "SOURCE",
        "1",
        "String",
        "Source database and study, e.g. cBioPortal CCLE ccle_broad_2019",
    ),
    InfoField(
        "ORIGINAL_BUILD",
        "1",
        "String",
        "Genome build of the source coordinates before any liftover",
    ),
    InfoField(
        "ORIGINAL_LOCUS",
        "1",
        "String",
        "Original locus before liftover in UCSC position format (chrom:pos) on ORIGINAL_BUILD",
    ),
    InfoField(
        "LIFTED",
        "0",
        "Flag",
        "Coordinate was lifted from ORIGINAL_BUILD to the output reference build",
    ),
    InfoField(
        "ANCHOR",
        "1",
        "String",
        "Provenance of the indel anchor base: reference (from --reference) or placeholder (N)",
    ),
)
"""The complete, ordered INFO schema emitted by cellme."""


@dataclass(frozen=True)
class TrackContext:
    """The per-run context shared by every record in a truth-track VCF."""

    cell_line: str
    sample_id: str
    study: str
    source_build: GenomeBuild
    target_build: GenomeBuild


@dataclass(frozen=True)
class VcfRecord:
    """A single VCF record with its 1-based position and resolved alleles."""

    contig: str
    position: int
    identifier: str | None
    reference_allele: str
    alternate_allele: str
    info: dict[str, InfoValue]


def _identifier(mutation: Mutation) -> str | None:
    """
    Build a stable, meaningful VCF identifier for a mutation.

    Args:
        mutation: The source mutation.

    Returns:
        A ``gene:proteinChange`` token when available, a ``gene:chrom:pos`` token
        when only the gene is known, otherwise None.
    """
    if mutation.gene and mutation.protein_change:
        return f"{mutation.gene}:{mutation.protein_change}"
    if mutation.gene:
        return f"{mutation.gene}:{mutation.chromosome}:{mutation.start_position}"
    return None


def _format_protein_position(mutation: Mutation) -> str | None:
    """
    Format the affected protein position or range for the PROTEIN_POS field.

    Args:
        mutation: The source mutation.

    Returns:
        A single position, a hyphenated range, or None when unknown.
    """
    start = mutation.protein_pos_start
    end = mutation.protein_pos_end
    if start is None:
        return None
    if end is None or end == start:
        return str(start)
    return f"{start}-{end}"


def _original_locus(mutation: Mutation) -> str:
    """
    Format the pre-liftover locus in UCSC position syntax.

    CCLE contigs arrive without a ``chr`` prefix, so one is added. A variant
    spanning a single base is rendered as ``chrom:pos``; a wider variant is
    rendered as the range ``chrom:start-end`` using the original coordinates.

    Args:
        mutation: The source mutation, on the source build.

    Returns:
        The UCSC-style locus string, for example ``chr17:7577022``.
    """
    if mutation.start_position == mutation.end_position:
        return f"chr{mutation.chromosome}:{mutation.start_position}"
    return f"chr{mutation.chromosome}:{mutation.start_position}-{mutation.end_position}"


def _build_info(
    mutation: Mutation,
    context: TrackContext,
    *,
    lifted: bool,
    anchor_source: str | None,
) -> dict[str, InfoValue]:
    """
    Assemble the INFO mapping for one record from its mutation and context.

    Args:
        mutation: The source mutation.
        context: The shared per-run context.
        lifted: Whether the coordinate was lifted to the target build.
        anchor_source: The indel anchor provenance, or None for substitutions.

    Returns:
        A mapping of INFO keys to values, omitting fields with no value.
    """
    info: dict[str, InfoValue] = {}
    if mutation.gene:
        info["GENE"] = mutation.gene
    if mutation.protein_change:
        info["PROTEIN_CHANGE"] = mutation.protein_change
    if mutation.variant_class:
        info["VARIANT_CLASS"] = mutation.variant_class
    if mutation.variant_type:
        info["VARIANT_TYPE"] = mutation.variant_type
    info["CELL_LINE"] = context.cell_line
    info["SAMPLE_ID"] = context.sample_id
    if mutation.entrez_gene_id is not None:
        info["ENTREZ"] = mutation.entrez_gene_id
    if mutation.refseq_mrna_id:
        info["REFSEQ"] = mutation.refseq_mrna_id
    protein_position = _format_protein_position(mutation)
    if protein_position is not None:
        info["PROTEIN_POS"] = protein_position
    info["SOURCE"] = f"cBioPortal CCLE {context.study}"
    info["ORIGINAL_BUILD"] = str(context.source_build)
    if lifted:
        info["ORIGINAL_LOCUS"] = _original_locus(mutation)
        info["LIFTED"] = True
    if anchor_source is not None:
        info["ANCHOR"] = anchor_source
    return info


def build_record(
    mutation: Mutation,
    context: TrackContext,
    *,
    lift_position: LiftPosition | None,
    anchor_base: AnchorBase | None,
) -> VcfRecord | None:
    """
    Convert one CCLE mutation into a VCF record on the target build.

    Substitutions carry their CCLE alleles through unchanged. Insertions and
    deletions are left-anchored: a deletion sits one base to the left of the
    deleted sequence and an insertion at the base preceding the insertion point.
    The anchoring position is computed in the source build, then lifted, so the
    anchor base is read from the target reference when one is supplied.

    Args:
        mutation: The source mutation to convert.
        context: The shared per-run context.
        lift_position: A callable that lifts a source position to the target
            build, or None when the source and target builds are the same.
        anchor_base: A callable returning the target reference base for indel
            anchoring, or None to use a placeholder ``N`` anchor.

    Returns:
        The converted record, or None when the mutation is malformed or its
        coordinate cannot be lifted to the target build.
    """
    is_deletion = mutation.variant_allele == _DASH
    is_insertion = mutation.reference_allele == _DASH
    if is_deletion and is_insertion:
        return None

    if is_deletion:
        source_position = mutation.start_position - 1
        if source_position < 1:
            return None
    else:
        source_position = mutation.start_position

    if lift_position is None:
        position = source_position
    else:
        lifted_position = lift_position(mutation.chromosome, source_position)
        if lifted_position is None:
            return None
        position = lifted_position

    anchor_source: str | None = None
    if is_deletion or is_insertion:
        if anchor_base is None:
            anchor = "N"
            anchor_source = "placeholder"
        else:
            anchor = anchor_base(mutation.chromosome, position)
            anchor_source = "reference"
        if is_deletion:
            reference_allele = anchor + mutation.reference_allele
            alternate_allele = anchor
        else:
            reference_allele = anchor
            alternate_allele = anchor + mutation.variant_allele
    else:
        reference_allele = mutation.reference_allele
        alternate_allele = mutation.variant_allele
        if not reference_allele or not alternate_allele:
            return None

    info = _build_info(
        mutation,
        context,
        lifted=lift_position is not None,
        anchor_source=anchor_source,
    )
    return VcfRecord(
        contig=mutation.chromosome,
        position=position,
        identifier=_identifier(mutation),
        reference_allele=reference_allele.upper(),
        alternate_allele=alternate_allele.upper(),
        info=info,
    )


def build_records(
    mutations: Iterable[Mutation],
    context: TrackContext,
    *,
    lift_position: LiftPosition | None,
    anchor_base: AnchorBase | None,
) -> tuple[list[VcfRecord], int]:
    """
    Convert many mutations, sorted into the target build's contig order.

    Args:
        mutations: The source mutations.
        context: The shared per-run context.
        lift_position: A callable that lifts a source position to the target
            build, or None when the source and target builds are the same.
        anchor_base: A callable returning the target reference base for indel
            anchoring, or None to use a placeholder ``N`` anchor.

    Returns:
        A tuple of the sorted records and the count of mutations that were
        dropped because they were malformed, could not be lifted, or sit on a
        contig that is not part of the target build.
    """
    order = contig_order(context.target_build)
    records: list[VcfRecord] = []
    dropped = 0
    for mutation in mutations:
        record = build_record(
            mutation, context, lift_position=lift_position, anchor_base=anchor_base
        )
        if record is None or record.contig not in order:
            dropped += 1
            continue
        records.append(record)
    records.sort(key=lambda record: (order[record.contig], record.position))
    return records, dropped


def build_header(context: TrackContext, version: str) -> "pysam.VariantHeader":
    """
    Build a VCF header for the target build with the full INFO schema.

    Args:
        context: The shared per-run context.
        version: The cellme version string, written to the ``source`` line.

    Returns:
        A pysam variant header populated with meta, contig, and INFO lines.
    """
    header = pysam.VariantHeader()
    header.add_line(f"##source=cellme {version}")
    header.add_line(f"##reference={context.target_build}")
    header.add_line(f"##cellme_cellLine={context.cell_line}")
    header.add_line(f"##cellme_sampleId={context.sample_id}")
    header.add_line(f"##cellme_sourceStudy=cBioPortal CCLE {context.study}")
    header.add_line(f"##cellme_sourceBuild={context.source_build}")
    for name, length in contigs_for(context.target_build):
        header.add_line(f"##contig=<ID={name},length={length}>")
    for field in INFO_FIELDS:
        header.add_line(
            f"##INFO=<ID={field.key},Number={field.number},"
            f'Type={field.type},Description="{field.description}">'
        )
    return header


def _to_pysam_record(header: "pysam.VariantHeader", record: VcfRecord) -> "pysam.VariantRecord":
    """
    Materialize a VcfRecord as a pysam record bound to a header.

    Args:
        header: The header the record will be written under.
        record: The record to materialize.

    Returns:
        A pysam variant record ready to be written.
    """
    start = record.position - 1
    stop = start + len(record.reference_allele)
    return header.new_record(
        contig=record.contig,
        start=start,
        stop=stop,
        alleles=(record.reference_allele, record.alternate_allele),
        id=record.identifier,
        info=dict(record.info),
    )


def write_vcf(
    records: Iterable[VcfRecord],
    header: "pysam.VariantHeader",
    output: Path | None,
) -> None:
    """
    Write records to a VCF at a path, or to standard output when no path given.

    Args:
        records: The records to write, already in sorted order.
        header: The header to write them under.
        output: The destination path, or None to write to standard output.
    """
    destination = str(output) if output is not None else "-"
    with pysam.VariantFile(destination, "w", header=header) as out:
        for record in records:
            out.write(_to_pysam_record(header, record))


def make_lifter(source: GenomeBuild, target: GenomeBuild) -> LiftPosition | None:
    """
    Create a position lifter between two builds, or None when they match.

    The returned callable drops coordinates that fail to lift or that lift to a
    different contig, so that only confidently mapped positions are emitted.

    Args:
        source: The build the source coordinates are on.
        target: The build to lift coordinates onto.

    Returns:
        A lifter callable, or None when ``source`` and ``target`` are equal.
    """
    if source == target:
        return None
    lifter = get_lifter(source.ucsc_name, target.ucsc_name)

    def lift(chromosome: str, position: int) -> int | None:
        result = lifter.query(chromosome, position)
        if not result:
            return None
        lifted_chromosome, lifted_position, _ = result[0]
        if lifted_chromosome.removeprefix("chr") != chromosome.removeprefix("chr"):
            return None
        return int(lifted_position)

    return lift


def make_anchor_base(reference: Path | None) -> AnchorBase | None:
    """
    Open a reference FASTA and return a base lookup for indel anchoring.

    The FASTA is matched by contig name with and without a ``chr`` prefix, so a
    reference using either naming convention works.

    Args:
        reference: The path to a FASTA for the target build, or None.

    Returns:
        A callable returning the reference base at a 1-based position, or None
        when no reference was supplied.
    """
    if reference is None:
        return None
    fasta = pysam.FastaFile(str(reference))
    contigs = set(fasta.references)

    def anchor(chromosome: str, position: int) -> str:
        for name in (chromosome, f"chr{chromosome}"):
            if name in contigs:
                base = fasta.fetch(name, position - 1, position)
                return base.upper() if base else "N"
        return "N"

    return anchor
