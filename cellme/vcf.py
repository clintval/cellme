"""Conversion of CCLE mutations into a sorted, well-described truth-track VCF."""

import logging
from collections.abc import Callable
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path
from typing import ClassVar

import pysam
from liftover import get_lifter

from cellme.builds import GenomeBuild
from cellme.builds import contig_order
from cellme.builds import contigs_for
from cellme.cbioportal import Mutation

logger = logging.getLogger("cellme")


@dataclass(frozen=True)
class LiftedCoordinate:
    """A coordinate lifted to the target build, with the strand of its chain block."""

    position: int
    """The 1-based position on the target build."""

    negative_strand: bool
    """True when the source maps to the minus strand, so alleles must be reoriented."""


LiftPosition = Callable[[str, int], LiftedCoordinate | None]
"""A callable mapping a source (contig, 1-based position) to a lifted coordinate."""

AnchorBase = Callable[[str, int], str]
"""A callable returning the single reference base at a (contig, 1-based position)."""

_COMPLEMENTS: dict[str, str] = {
    "A": "T",
    "C": "G",
    "G": "C",
    "T": "A",
    "U": "A",
    "M": "K",
    "K": "M",
    "R": "Y",
    "Y": "R",
    "W": "W",
    "S": "S",
    "B": "V",
    "V": "B",
    "H": "D",
    "D": "H",
    "N": "N",
    "a": "t",
    "c": "g",
    "g": "c",
    "t": "a",
    "u": "a",
    "m": "k",
    "k": "m",
    "r": "y",
    "y": "r",
    "w": "w",
    "s": "s",
    "b": "v",
    "v": "b",
    "h": "d",
    "d": "h",
    "n": "n",
}
"""Complement of every IUPAC nucleotide code, in both cases, including the N placeholder."""

_INVALID_BASES: str = "".join(c for c in (chr(o) for o in range(256)) if c not in _COMPLEMENTS)
"""Every byte that is not an IUPAC nucleotide code, dropped during complementing."""

_COMPLEMENTS_TABLE: dict[int, int | None] = str.maketrans(
    "".join(_COMPLEMENTS.keys()), "".join(_COMPLEMENTS.values()), _INVALID_BASES
)
"""Translation table mapping each base to its complement for fast reverse-complementing."""


def _reverse_complement(bases: str) -> str:
    """
    Return the reverse complement of a DNA sequence over the IUPAC alphabet.

    Alleles lifted across a minus-strand chain block are read on the opposite
    strand, so both their order and their bases must be flipped to sit on the
    target build's forward strand. Degenerate IUPAC codes are complemented too.

    Args:
        bases: A DNA sequence over the IUPAC nucleotide codes.

    Returns:
        The reverse complement of the sequence.

    Raises:
        KeyError: When the sequence contains a character that is not an IUPAC
            nucleotide code, which would otherwise be silently dropped.
    """
    reverse_complemented = bases.translate(_COMPLEMENTS_TABLE)[::-1]
    if len(reverse_complemented) != len(bases):
        bad_bases = "".join({base for base in bases if base not in _COMPLEMENTS})
        raise KeyError(f"Invalid bases found: {bad_bases}")
    return reverse_complemented


ReferenceLookup = Callable[[str, int, int], str]
"""A callable returning the reference bases over a 1-based, inclusive ``[start, end]`` span."""

_DASH: str = "-"
"""The MAF sentinel for the absent allele of an insertion or deletion."""

InfoValue = str | int | bool
"""The value types cellme writes into a VCF INFO field."""


class TruthTrackError(Exception):
    """
    Base class for failures raised while building truth-track VCF records.

    A subclass may record the command-line flag that downgrades the failure from
    a raised error to a dropped-and-warned record, so the command line can point
    the user at the matching opt-out. It is left empty when dropping is already
    the default and no such opt-out flag exists.
    """

    opt_out_flag: ClassVar[str] = ""
    """The command-line flag that turns this error into a dropped, warned record."""


class LiftoverError(TruthTrackError):
    """
    Raised when a mutation's coordinate cannot be lifted to the target build.

    Liftover failures are lenient by default and only surface as a raised error
    when the run opts into strict liftover with ``--raise-on-liftover-fails``, so
    there is no opt-out flag to advertise.
    """


class ReferenceMismatchError(TruthTrackError):
    """Raised when a record's REF allele disagrees with the target reference."""

    opt_out_flag: ClassVar[str] = "--skip-ref-mismatch"


def _describe(mutation: Mutation) -> str:
    """
    Render a short, human-readable label for a mutation for use in messages.

    Args:
        mutation: The source mutation.

    Returns:
        A ``gene proteinChange`` label when both are known, the gene alone when
        only it is known, otherwise the literal ``variant``.
    """
    if mutation.gene and mutation.protein_change:
        return f"{mutation.gene} {mutation.protein_change}"
    if mutation.gene:
        return mutation.gene
    return "variant"


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
    info["ORIGINAL_BUILD"] = context.source_build.grch_name
    if lifted:
        info["ORIGINAL_LOCUS"] = _original_locus(mutation)
        info["LIFTED"] = True
    if anchor_source is not None:
        info["ANCHOR"] = anchor_source
    return info


def _lift(
    lift_position: LiftPosition | None,
    mutation: Mutation,
    context: TrackContext,
    position: int,
) -> LiftedCoordinate:
    """
    Lift one source position to the target build, or raise :class:`LiftoverError`.

    Args:
        lift_position: The lifter, or None when the source and target builds match.
        mutation: The source mutation, used only to describe a failure.
        context: The shared per-run context, used to name the builds in a failure.
        position: The 1-based source position to lift.

    Returns:
        The lifted coordinate. When no lifter is supplied the position passes
        through unchanged on the forward strand.

    Raises:
        LiftoverError: When the position cannot be mapped to the target build.
    """
    if lift_position is None:
        return LiftedCoordinate(position=position, negative_strand=False)
    lifted = lift_position(mutation.chromosome, position)
    if lifted is None:
        raise LiftoverError(
            f"Could not lift {_describe(mutation)} at {_original_locus(mutation)} "
            f"from {context.source_build.grch_name} to {context.target_build.grch_name}."
        )
    return lifted


def _resolve_anchor(
    anchor_base: AnchorBase | None, chromosome: str, position: int
) -> tuple[str, str]:
    """
    Resolve the indel anchor base and its provenance at a target position.

    Args:
        anchor_base: A callable returning the target reference base, or None.
        chromosome: The contig of the anchor position.
        position: The 1-based anchor position on the target build.

    Returns:
        A tuple of the anchor base and its provenance, ``reference`` when read
        from the target FASTA or ``placeholder`` when it falls back to ``N``.
    """
    if anchor_base is None:
        return "N", "placeholder"
    return anchor_base(chromosome, position), "reference"


def _lift_span(
    lift_position: LiftPosition | None,
    mutation: Mutation,
    context: TrackContext,
) -> tuple[int, int, bool]:
    """
    Lift both ends of a variant's source span onto the target build.

    Args:
        lift_position: The lifter, or None when the source and target builds match.
        mutation: The source mutation whose start and end positions are lifted.
        context: The shared per-run context, used to name the builds in a failure.

    Returns:
        A tuple of the lower and upper 1-based target positions and whether the
        span sits on the minus strand.

    Raises:
        LiftoverError: When either end fails to lift or the two ends disagree on
            strand, which means the span straddles a chain break.
    """
    start = _lift(lift_position, mutation, context, mutation.start_position)
    end = _lift(lift_position, mutation, context, mutation.end_position)
    if start.negative_strand != end.negative_strand:
        raise LiftoverError(
            f"Could not lift {_describe(mutation)} at {_original_locus(mutation)} "
            f"from {context.source_build.grch_name} to {context.target_build.grch_name}: "
            f"the variant straddles a strand boundary between builds."
        )
    low, high = sorted((start.position, end.position))
    return low, high, start.negative_strand


def build_record(
    mutation: Mutation,
    context: TrackContext,
    *,
    lift_position: LiftPosition | None,
    anchor_base: AnchorBase | None,
) -> VcfRecord | None:
    """
    Convert one CCLE mutation into a VCF record on the target build.

    Coordinates are lifted with strand awareness: a source that maps to a
    minus-strand chain block is placed at the lifted position of its left-most
    base on the target build, and its alleles are reverse-complemented so the
    emitted REF sits on the target build's forward strand. Insertions and
    deletions are left-anchored, so a deletion carries the reference base to the
    left of the deleted sequence and an insertion the base preceding the insertion
    point; the anchor base is read from the target reference when one is supplied.

    Args:
        mutation: The source mutation to convert.
        context: The shared per-run context.
        lift_position: A callable that lifts a source position to the target
            build, or None when the source and target builds are the same.
        anchor_base: A callable returning the target reference base for indel
            anchoring, or None to use a placeholder ``N`` anchor.

    Returns:
        The converted record, or None when the mutation is malformed.

    Raises:
        LiftoverError: When ``lift_position`` cannot map the coordinate to the
            target build. The message names the variant and its source locus.
    """
    is_deletion = mutation.variant_allele == _DASH
    is_insertion = mutation.reference_allele == _DASH
    if is_deletion and is_insertion:
        return None

    reference_allele = mutation.reference_allele.upper()
    variant_allele = mutation.variant_allele.upper()

    anchor_source: str | None = None
    if is_deletion:
        low, _high, negative = _lift_span(lift_position, mutation, context)
        position = low - 1
        if position < 1:
            return None
        deleted = _reverse_complement(reference_allele) if negative else reference_allele
        anchor, anchor_source = _resolve_anchor(anchor_base, mutation.chromosome, position)
        reference_allele = anchor + deleted
        variant_allele = anchor
    elif is_insertion:
        lifted = _lift(lift_position, mutation, context, mutation.start_position)
        position = lifted.position - 1 if lifted.negative_strand else lifted.position
        if position < 1:
            return None
        inserted = _reverse_complement(variant_allele) if lifted.negative_strand else variant_allele
        anchor, anchor_source = _resolve_anchor(anchor_base, mutation.chromosome, position)
        reference_allele = anchor
        variant_allele = anchor + inserted
    else:
        if not reference_allele or not variant_allele:
            return None
        low, _high, negative = _lift_span(lift_position, mutation, context)
        position = low
        if negative:
            reference_allele = _reverse_complement(reference_allele)
            variant_allele = _reverse_complement(variant_allele)

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
        reference_allele=reference_allele,
        alternate_allele=variant_allele,
        info=info,
    )


def validate_reference_allele(
    record: VcfRecord,
    mutation: Mutation,
    context: TrackContext,
    reference_lookup: ReferenceLookup,
) -> None:
    """
    Check that a record's REF allele matches the target reference sequence.

    The reference bases spanning the record's REF allele are read from the target
    FASTA and compared to the emitted REF. A single-base REF (a substitution or an
    insertion anchor) checks one base; a multi-base REF (a deletion's anchor base
    plus its deleted bases) checks the whole span.

    Args:
        record: The built record, already positioned on the target build.
        mutation: The source mutation, used only to describe a mismatch.
        context: The shared per-run context, used to name the target build.
        reference_lookup: A callable returning the target reference bases over a
            1-based, inclusive span.

    Raises:
        ReferenceMismatchError: When the REF allele disagrees with the reference.
            The message names the variant, the expected REF, and the actual bases.
    """
    expected = record.reference_allele
    end = record.position + len(expected) - 1
    actual = reference_lookup(record.contig, record.position, end)
    if actual != expected:
        raise ReferenceMismatchError(
            f"REF allele for {_describe(mutation)} at {record.contig}:{record.position} "
            f"on {context.target_build.grch_name} does not match the reference: "
            f"expected REF {expected!r} but the reference has {actual!r}."
        )


def build_records(
    mutations: Iterable[Mutation],
    context: TrackContext,
    *,
    lift_position: LiftPosition | None,
    anchor_base: AnchorBase | None,
    reference_lookup: ReferenceLookup | None = None,
    raise_on_liftover_fails: bool = False,
    skip_ref_mismatch: bool = False,
) -> tuple[list[VcfRecord], int]:
    """
    Convert many mutations, sorted into the target build's contig order.

    Liftover failures are lenient by default: a mutation that cannot be lifted to
    the target build is dropped with a warning, so a truth track still builds from
    the coordinates that do lift. Reference-base mismatches are strict by default:
    an emitted record whose REF disagrees with the supplied reference raises rather
    than being silently dropped, since a mismatch signals a miscalled record. Each
    default has an opt-in that flips it.

    Args:
        mutations: The source mutations.
        context: The shared per-run context.
        lift_position: A callable that lifts a source position to the target
            build, or None when the source and target builds are the same.
        anchor_base: A callable returning the target reference base for indel
            anchoring, or None to use a placeholder ``N`` anchor.
        reference_lookup: A callable returning target reference bases over a span,
            used to validate each record's REF allele, or None to skip that check.
        raise_on_liftover_fails: Raise :class:`LiftoverError` on a liftover failure
            instead of dropping the mutation with a warning.
        skip_ref_mismatch: Drop and warn on a reference-base mismatch instead of
            raising :class:`ReferenceMismatchError`.

    Returns:
        A tuple of the sorted records and the count of mutations that were
        dropped because they were malformed, sit on a contig that is not part of
        the target build, could not be lifted, or failed reference validation.

    Raises:
        LiftoverError: When a coordinate cannot be lifted and
            ``raise_on_liftover_fails`` is True.
        ReferenceMismatchError: When a record's REF disagrees with the reference
            and ``skip_ref_mismatch`` is False.
    """
    order = contig_order(context.target_build)
    records: list[VcfRecord] = []
    dropped = 0
    for mutation in mutations:
        try:
            record = build_record(
                mutation, context, lift_position=lift_position, anchor_base=anchor_base
            )
        except LiftoverError as error:
            if raise_on_liftover_fails:
                raise
            logger.warning(f"Dropping mutation: {error}")
            dropped += 1
            continue
        if record is None or record.contig not in order:
            dropped += 1
            continue
        if reference_lookup is not None:
            try:
                validate_reference_allele(record, mutation, context, reference_lookup)
            except ReferenceMismatchError as error:
                if not skip_ref_mismatch:
                    raise
                logger.warning(f"Dropping mutation: {error}")
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
    header.add_line(f"##reference={context.target_build.grch_name}")
    header.add_line(f"##cellme_cellLine={context.cell_line}")
    header.add_line(f"##cellme_sampleId={context.sample_id}")
    header.add_line(f"##cellme_sourceStudy=cBioPortal CCLE {context.study}")
    header.add_line(f"##cellme_sourceBuild={context.source_build.grch_name}")
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

    When the output path ends in ``.gz`` the VCF is block-gzip (BGZF) compressed
    and a tabix ``.tbi`` index is written alongside it, so the result can be used
    directly as a random-access truth track. A plain path is written uncompressed
    and unindexed, as is standard output.

    Args:
        records: The records to write, already in sorted order.
        header: The header to write them under.
        output: The destination path, or None to write to standard output.
    """
    if output is None:
        _write_records("-", header, records, compressed=False)
        return
    destination = str(output)
    compressed = destination.endswith(".gz")
    _write_records(destination, header, records, compressed=compressed)
    if compressed:
        pysam.tabix_index(destination, preset="vcf", force=True)


def _write_records(
    destination: str,
    header: "pysam.VariantHeader",
    records: Iterable[VcfRecord],
    *,
    compressed: bool,
) -> None:
    """
    Write records to a destination, BGZF-compressed when ``compressed`` is set.

    Args:
        destination: The pysam destination, a path or ``-`` for standard output.
        header: The header to write the records under.
        records: The records to write, already in sorted order.
        compressed: Whether to write a block-gzip compressed VCF.
    """
    mode = "wz" if compressed else "w"
    with pysam.VariantFile(destination, mode, header=header) as out:
        for record in records:
            out.write(_to_pysam_record(header, record))


def make_lifter(source: GenomeBuild, target: GenomeBuild) -> LiftPosition | None:
    """
    Create a position lifter between two builds, or None when they match.

    The lifter works in 1-based coordinates and reports the strand of each match,
    so the caller can reverse-complement alleles that map to a minus-strand chain
    block. The returned callable returns None for a coordinate that fails to lift
    or that lifts to a different contig (including an alt contig), so only
    confidently mapped positions carry through; the caller decides how a None is
    handled.

    Args:
        source: The build the source coordinates are on.
        target: The build to lift coordinates onto.

    Returns:
        A lifter callable, or None when ``source`` and ``target`` are equal.
    """
    if source == target:
        return None
    lifter = get_lifter(source.ucsc_name, target.ucsc_name, one_based=True)

    def lift(chromosome: str, position: int) -> LiftedCoordinate | None:
        result = lifter.query(chromosome, position)
        if not result:
            return None
        lifted_chromosome, lifted_position, strand = result[0]
        if lifted_chromosome.removeprefix("chr") != chromosome.removeprefix("chr"):
            return None
        return LiftedCoordinate(position=int(lifted_position), negative_strand=strand == "-")

    return lift


def make_reference_lookup(reference: Path | None) -> ReferenceLookup | None:
    """
    Open a reference FASTA and return a lookup over 1-based, inclusive spans.

    The FASTA is matched by contig name with and without a ``chr`` prefix, so a
    reference using either naming convention works. A contig absent from the
    FASTA, or a span past its end, yields an empty string.

    Args:
        reference: The path to a FASTA for the target build, or None.

    Returns:
        A callable returning the uppercased reference bases over a span, or None
        when no reference was supplied.
    """
    if reference is None:
        return None
    fasta = pysam.FastaFile(str(reference))
    contigs = set(fasta.references)

    def lookup(chromosome: str, start: int, end: int) -> str:
        for name in (chromosome, f"chr{chromosome}"):
            if name in contigs:
                return fasta.fetch(name, start - 1, end).upper()
        return ""

    return lookup


def make_anchor_base(reference_lookup: ReferenceLookup | None) -> AnchorBase | None:
    """
    Adapt a reference lookup into a single-base anchor lookup for indels.

    Args:
        reference_lookup: A span lookup from :func:`make_reference_lookup`, or None.

    Returns:
        A callable returning the reference base at a 1-based position, falling
        back to a placeholder ``N`` when the base is unavailable, or None when no
        reference lookup was supplied.
    """
    if reference_lookup is None:
        return None

    def anchor(chromosome: str, position: int) -> str:
        base = reference_lookup(chromosome, position, position)
        return base if base else "N"

    return anchor
