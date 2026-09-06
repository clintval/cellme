import logging
from collections.abc import Callable
from pathlib import Path
from typing import Any

import pysam
import pytest

from cellme.builds import GenomeBuild
from cellme.cbioportal import Mutation
from cellme.vcf import INFO_FIELDS
from cellme.vcf import LiftedCoordinate
from cellme.vcf import LiftoverError
from cellme.vcf import ReferenceLookup
from cellme.vcf import ReferenceMismatchError
from cellme.vcf import TrackContext
from cellme.vcf import VcfRecord
from cellme.vcf import build_header
from cellme.vcf import build_record
from cellme.vcf import build_records
from cellme.vcf import make_anchor_base
from cellme.vcf import make_reference_lookup
from cellme.vcf import validate_reference_allele
from cellme.vcf import write_vcf

CONTEXT = TrackContext(
    cell_line="MOLT4",
    sample_id="MOLT4_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE",
    study="ccle_broad_2019",
    source_build=GenomeBuild.hg19,
    target_build=GenomeBuild.hg38,
)

SAME_BUILD_CONTEXT = TrackContext(
    cell_line="MOLT4",
    sample_id="MOLT4_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE",
    study="ccle_broad_2019",
    source_build=GenomeBuild.hg19,
    target_build=GenomeBuild.hg19,
)


def make_mutation(**overrides: Any) -> Mutation:
    defaults: dict[str, Any] = {
        "gene": "TP53",
        "entrez_gene_id": 7157,
        "chromosome": "17",
        "start_position": 7577022,
        "end_position": 7577022,
        "reference_allele": "G",
        "variant_allele": "A",
        "protein_change": "R306*",
        "variant_class": "Nonsense_Mutation",
        "variant_type": "SNP",
        "ncbi_build": "GRCh37",
        "refseq_mrna_id": "NM_001126112.2",
        "protein_pos_start": 306,
        "protein_pos_end": 306,
    }
    defaults.update(overrides)
    return Mutation(**defaults)


def add_offset(_chromosome: str, position: int) -> LiftedCoordinate:
    return LiftedCoordinate(position=position + 100000, negative_strand=False)


def minus_lifter(mirror: int) -> Callable[[str, int], LiftedCoordinate]:
    """
    Return a minus-strand lifter mapping source position ``p`` to ``mirror - p``.

    A minus-strand chain block reverses coordinates, so consecutive source
    positions map to consecutive decreasing target positions. Choosing ``mirror``
    as ``source_position + target_position`` places a chosen source at a chosen
    target while preserving that reversal for its neighbors.
    """

    def lift(_chromosome: str, position: int) -> LiftedCoordinate:
        return LiftedCoordinate(position=mirror - position, negative_strand=True)

    return lift


def never_lifts(_chromosome: str, _position: int) -> None:
    return None


def anchor_is_c(_chromosome: str, _position: int) -> str:
    return "C"


def test_substitution_carries_alleles_through_and_flags_lift() -> None:
    record = build_record(make_mutation(), CONTEXT, lift_position=add_offset, anchor_base=None)
    assert record is not None
    assert record.contig == "17"
    assert record.position == 7677022
    assert record.reference_allele == "G"
    assert record.alternate_allele == "A"
    assert record.identifier == "TP53:R306*"
    assert record.info["LIFTED"] is True
    assert record.info["ORIGINAL_BUILD"] == "GRCh37"
    assert "ANCHOR" not in record.info


def test_same_build_does_not_lift_or_flag() -> None:
    record = build_record(make_mutation(), SAME_BUILD_CONTEXT, lift_position=None, anchor_base=None)
    assert record is not None
    assert record.position == 7577022
    assert "LIFTED" not in record.info


def test_deletion_left_anchors_with_placeholder() -> None:
    mutation = make_mutation(
        gene="PTEN",
        chromosome="10",
        start_position=89717770,
        end_position=89717770,
        reference_allele="A",
        variant_allele="-",
        protein_change="K267Rfs*9",
        variant_class="Frame_Shift_Del",
        variant_type="DEL",
        protein_pos_start=265,
        protein_pos_end=267,
    )
    record = build_record(mutation, CONTEXT, lift_position=add_offset, anchor_base=None)
    assert record is not None
    assert record.position == 89717769 + 100000
    assert record.reference_allele == "NA"
    assert record.alternate_allele == "N"
    assert record.info["ANCHOR"] == "placeholder"
    assert record.info["PROTEIN_POS"] == "265-267"


def test_deletion_uses_reference_anchor_base() -> None:
    mutation = make_mutation(
        reference_allele="A",
        variant_allele="-",
        variant_class="Frame_Shift_Del",
        variant_type="DEL",
    )
    record = build_record(mutation, CONTEXT, lift_position=add_offset, anchor_base=anchor_is_c)
    assert record is not None
    assert record.reference_allele == "CA"
    assert record.alternate_allele == "C"
    assert record.info["ANCHOR"] == "reference"


def test_insertion_left_anchors_with_placeholder() -> None:
    mutation = make_mutation(
        reference_allele="-",
        variant_allele="T",
        variant_class="Frame_Shift_Ins",
        variant_type="INS",
    )
    record = build_record(mutation, CONTEXT, lift_position=add_offset, anchor_base=None)
    assert record is not None
    assert record.reference_allele == "N"
    assert record.alternate_allele == "NT"
    assert record.info["ANCHOR"] == "placeholder"


def test_minus_strand_substitution_lifts_and_reverse_complements() -> None:
    # AVPR1B R364C: hg19 chr1:206230957 C>T lifts to hg38 chr1:206110374 G>A (CrossMap oracle).
    mutation = make_mutation(
        gene="AVPR1B",
        chromosome="1",
        start_position=206230957,
        end_position=206230957,
        reference_allele="C",
        variant_allele="T",
        protein_change="R364C",
        variant_class="Missense_Mutation",
        variant_type="SNP",
    )
    mirror = 206230957 + 206110374
    record = build_record(mutation, CONTEXT, lift_position=minus_lifter(mirror), anchor_base=None)
    assert record is not None
    assert record.position == 206110374
    assert record.reference_allele == "G"
    assert record.alternate_allele == "A"
    assert record.info["LIFTED"] is True


def test_minus_strand_substitution_passes_reference_validation() -> None:
    # The reverse-complemented REF (G) must match the hg38 base at the lifted position.
    mutation = make_mutation(
        gene="AVPR1B",
        chromosome="1",
        start_position=206230957,
        end_position=206230957,
        reference_allele="C",
        variant_allele="T",
        protein_change="R364C",
        variant_type="SNP",
    )
    mirror = 206230957 + 206110374
    record = build_record(mutation, CONTEXT, lift_position=minus_lifter(mirror), anchor_base=None)
    assert record is not None
    reference_lookup = make_fake_lookup("G", 206110374)
    validate_reference_allele(record, mutation, CONTEXT, reference_lookup)


def test_minus_strand_deletion_left_anchors_and_reverse_complements() -> None:
    # Deleting hg19 chr1:206230957 (C) lands at hg38 chr1:206110373 with REF=CG, ALT=C,
    # matching the hg38 reference bases CG at chr1:206110373-206110374.
    mutation = make_mutation(
        gene="AVPR1B",
        chromosome="1",
        start_position=206230957,
        end_position=206230957,
        reference_allele="C",
        variant_allele="-",
        variant_class="Frame_Shift_Del",
        variant_type="DEL",
    )
    mirror = 206230957 + 206110374

    def anchor_is_reference_c(_chromosome: str, _position: int) -> str:
        return "C"

    record = build_record(
        mutation, CONTEXT, lift_position=minus_lifter(mirror), anchor_base=anchor_is_reference_c
    )
    assert record is not None
    assert record.position == 206110373
    assert record.reference_allele == "CG"
    assert record.alternate_allele == "C"
    assert record.info["ANCHOR"] == "reference"
    reference_lookup = make_fake_lookup("CG", 206110373)
    validate_reference_allele(record, mutation, CONTEXT, reference_lookup)


def test_minus_strand_insertion_left_anchors_and_reverse_complements() -> None:
    # RBP3 insertion of G between hg19 chr10:48390168-48390169 lands at hg38 chr10:47349193
    # with REF=G, ALT=GC (revcomp of the inserted G), matching the hg38 base G there.
    mutation = make_mutation(
        gene="RBP3",
        chromosome="10",
        start_position=48390168,
        end_position=48390169,
        reference_allele="-",
        variant_allele="G",
        variant_class="Frame_Shift_Ins",
        variant_type="INS",
    )
    mirror = 48390168 + 47349194

    def anchor_is_reference_g(_chromosome: str, _position: int) -> str:
        return "G"

    record = build_record(
        mutation, CONTEXT, lift_position=minus_lifter(mirror), anchor_base=anchor_is_reference_g
    )
    assert record is not None
    assert record.position == 47349193
    assert record.reference_allele == "G"
    assert record.alternate_allele == "GC"
    assert record.info["ANCHOR"] == "reference"
    reference_lookup = make_fake_lookup("G", 47349193)
    validate_reference_allele(record, mutation, CONTEXT, reference_lookup)


def test_unliftable_mutation_raises_and_names_the_variant() -> None:
    with pytest.raises(LiftoverError) as excinfo:
        build_record(make_mutation(), CONTEXT, lift_position=never_lifts, anchor_base=None)
    message = str(excinfo.value)
    assert "TP53 R306*" in message
    assert "chr17:7577022" in message
    assert "GRCh37" in message
    assert "GRCh38" in message


def test_identifier_falls_back_to_locus_without_protein_change() -> None:
    record = build_record(
        make_mutation(protein_change=None),
        SAME_BUILD_CONTEXT,
        lift_position=None,
        anchor_base=None,
    )
    assert record is not None
    assert record.identifier == "TP53:17:7577022"


def test_build_records_sorts_and_counts_drops() -> None:
    mutations = [
        make_mutation(chromosome="17", start_position=7577022),
        make_mutation(chromosome="1", start_position=1000),
        make_mutation(chromosome="unmappable", start_position=5),
    ]

    def lift(chromosome: str, position: int) -> LiftedCoordinate | None:
        if chromosome == "unmappable":
            return None
        return LiftedCoordinate(position=position, negative_strand=False)

    records, dropped = build_records(mutations, CONTEXT, lift_position=lift, anchor_base=None)
    assert dropped == 1
    assert [record.contig for record in records] == ["1", "17"]


def test_lifted_record_carries_original_locus_in_ucsc_format() -> None:
    record = build_record(make_mutation(), CONTEXT, lift_position=add_offset, anchor_base=None)
    assert record is not None
    assert record.info["ORIGINAL_LOCUS"] == "chr17:7577022"


def test_passthrough_record_has_no_original_locus() -> None:
    record = build_record(make_mutation(), SAME_BUILD_CONTEXT, lift_position=None, anchor_base=None)
    assert record is not None
    assert "ORIGINAL_LOCUS" not in record.info


def test_multi_base_variant_original_locus_is_a_range() -> None:
    mutation = make_mutation(
        reference_allele="AT",
        variant_allele="-",
        start_position=100,
        end_position=101,
        variant_class="Frame_Shift_Del",
        variant_type="DEL",
    )
    record = build_record(mutation, CONTEXT, lift_position=add_offset, anchor_base=None)
    assert record is not None
    assert record.info["ORIGINAL_LOCUS"] == "chr17:100-101"


def test_build_records_drops_unknown_contigs() -> None:
    mutations = [
        make_mutation(chromosome="17"),
        make_mutation(chromosome="GL000209"),
    ]
    records, dropped = build_records(
        mutations, SAME_BUILD_CONTEXT, lift_position=None, anchor_base=None
    )
    assert dropped == 1
    assert [record.contig for record in records] == ["17"]


def test_build_records_drops_liftover_failure_by_default(
    caplog: pytest.LogCaptureFixture,
) -> None:
    mutations = [make_mutation()]
    with caplog.at_level(logging.WARNING, logger="cellme"):
        records, dropped = build_records(
            mutations,
            CONTEXT,
            lift_position=never_lifts,
            anchor_base=None,
        )
    assert records == []
    assert dropped == 1
    assert "TP53 R306*" in caplog.text


def test_build_records_raises_on_liftover_failure_with_flag() -> None:
    mutations = [make_mutation(chromosome="1", start_position=1000)]
    with pytest.raises(LiftoverError):
        build_records(
            mutations,
            CONTEXT,
            lift_position=never_lifts,
            anchor_base=None,
            raise_on_liftover_fails=True,
        )


def make_fake_lookup(sequence: str, start: int) -> ReferenceLookup:
    """Return a reference lookup where `sequence` begins at 1-based `start` on any contig."""

    def lookup(_chromosome: str, span_start: int, span_end: int) -> str:
        return sequence[span_start - start : span_end - start + 1]

    return lookup


def test_validate_reference_allele_accepts_matching_snv() -> None:
    record = build_record(make_mutation(), CONTEXT, lift_position=add_offset, anchor_base=None)
    assert record is not None
    reference_lookup = make_fake_lookup("G", record.position)
    validate_reference_allele(record, make_mutation(), CONTEXT, reference_lookup)


def test_validate_reference_allele_raises_on_mismatched_snv() -> None:
    record = build_record(make_mutation(), CONTEXT, lift_position=add_offset, anchor_base=None)
    assert record is not None
    reference_lookup = make_fake_lookup("C", record.position)
    with pytest.raises(ReferenceMismatchError) as excinfo:
        validate_reference_allele(record, make_mutation(), CONTEXT, reference_lookup)
    message = str(excinfo.value)
    assert "TP53 R306*" in message
    assert "'G'" in message
    assert "'C'" in message
    assert "GRCh38" in message


def test_validate_reference_allele_accepts_matching_deletion() -> None:
    mutation = make_mutation(
        reference_allele="A",
        variant_allele="-",
        variant_class="Frame_Shift_Del",
        variant_type="DEL",
    )
    record = build_record(mutation, CONTEXT, lift_position=add_offset, anchor_base=anchor_is_c)
    assert record is not None
    assert record.reference_allele == "CA"
    reference_lookup = make_fake_lookup("CA", record.position)
    validate_reference_allele(record, mutation, CONTEXT, reference_lookup)


def test_validate_reference_allele_raises_on_mismatched_deletion() -> None:
    mutation = make_mutation(
        reference_allele="A",
        variant_allele="-",
        variant_class="Frame_Shift_Del",
        variant_type="DEL",
    )
    record = build_record(mutation, CONTEXT, lift_position=add_offset, anchor_base=anchor_is_c)
    assert record is not None
    reference_lookup = make_fake_lookup("CT", record.position)
    with pytest.raises(ReferenceMismatchError) as excinfo:
        validate_reference_allele(record, mutation, CONTEXT, reference_lookup)
    assert "'CA'" in str(excinfo.value)
    assert "'CT'" in str(excinfo.value)


def test_build_records_raises_on_reference_mismatch_by_default() -> None:
    mutations = [make_mutation(start_position=7577022, end_position=7577022)]
    reference_lookup = make_fake_lookup("C", 7677022)
    with pytest.raises(ReferenceMismatchError):
        build_records(
            mutations,
            CONTEXT,
            lift_position=add_offset,
            anchor_base=None,
            reference_lookup=reference_lookup,
        )


def test_build_records_skips_reference_mismatch_with_flag(
    caplog: pytest.LogCaptureFixture,
) -> None:
    mutations = [make_mutation(start_position=7577022, end_position=7577022)]
    reference_lookup = make_fake_lookup("C", 7677022)
    with caplog.at_level(logging.WARNING, logger="cellme"):
        records, dropped = build_records(
            mutations,
            CONTEXT,
            lift_position=add_offset,
            anchor_base=None,
            reference_lookup=reference_lookup,
            skip_ref_mismatch=True,
        )
    assert records == []
    assert dropped == 1
    assert "TP53 R306*" in caplog.text


def test_build_records_happy_path_with_matching_reference_is_unchanged() -> None:
    mutations = [make_mutation(start_position=7577022, end_position=7577022)]
    reference_lookup = make_fake_lookup("G", 7677022)
    records, dropped = build_records(
        mutations,
        CONTEXT,
        lift_position=add_offset,
        anchor_base=None,
        reference_lookup=reference_lookup,
    )
    assert dropped == 0
    assert [record.reference_allele for record in records] == ["G"]


def test_make_reference_lookup_reads_span_from_fasta(tmp_path: Path) -> None:
    fasta = tmp_path / "ref.fa"
    fasta.write_text(">17\nACGTACGTACGT\n")
    pysam.faidx(str(fasta))
    reference_lookup = make_reference_lookup(fasta)
    assert reference_lookup is not None
    assert reference_lookup("17", 1, 4) == "ACGT"
    assert reference_lookup("17", 5, 5) == "A"
    assert reference_lookup("unplaced", 1, 4) == ""
    anchor_base = make_anchor_base(reference_lookup)
    assert anchor_base is not None
    assert anchor_base("17", 3) == "G"


def test_make_reference_lookup_matches_chr_prefixed_contigs(tmp_path: Path) -> None:
    fasta = tmp_path / "ref.fa"
    fasta.write_text(">chr17\nACGTACGTACGT\n")
    pysam.faidx(str(fasta))
    reference_lookup = make_reference_lookup(fasta)
    assert reference_lookup is not None
    assert reference_lookup("17", 1, 4) == "ACGT"


def test_make_reference_lookup_is_none_without_reference() -> None:
    assert make_reference_lookup(None) is None
    assert make_anchor_base(None) is None


def test_header_declares_reference_contigs_and_full_info_schema() -> None:
    header_text = str(build_header(CONTEXT, "0.1.0"))
    assert "##reference=GRCh38" in header_text
    assert "##source=cellme 0.1.0" in header_text
    assert "##contig=<ID=17,length=83257441>" in header_text
    for field in INFO_FIELDS:
        assert f"##INFO=<ID={field.key}," in header_text


def _one_record_and_header() -> tuple[list[VcfRecord], "pysam.VariantHeader"]:
    records, _dropped = build_records(
        [make_mutation()], SAME_BUILD_CONTEXT, lift_position=None, anchor_base=None
    )
    return records, build_header(SAME_BUILD_CONTEXT, "0.1.0")


def test_write_vcf_plain_path_is_uncompressed_and_unindexed(tmp_path: Path) -> None:
    records, header = _one_record_and_header()
    output = tmp_path / "molt4.vcf"

    write_vcf(records, header, output)

    assert output.read_bytes()[:2] != b"\x1f\x8b"  # not gzip-compressed
    assert not (tmp_path / "molt4.vcf.tbi").exists()
    with pysam.VariantFile(str(output)) as vcf:
        assert [record.info["GENE"] for record in vcf] == ["TP53"]


def test_write_vcf_gz_path_is_bgzipped_and_tabix_indexed(tmp_path: Path) -> None:
    records, header = _one_record_and_header()
    output = tmp_path / "molt4.vcf.gz"

    write_vcf(records, header, output)

    # BGZF is gzip with the extra-field flag (FLG.FEXTRA) set: magic 1f 8b 08 04.
    assert output.read_bytes()[:4] == b"\x1f\x8b\x08\x04"
    assert (tmp_path / "molt4.vcf.gz.tbi").exists()
    with pysam.VariantFile(str(output)) as vcf:
        assert [record.info["GENE"] for record in vcf] == ["TP53"]
