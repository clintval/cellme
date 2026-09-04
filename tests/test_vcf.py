from typing import Any

from cellme.builds import GenomeBuild
from cellme.cbioportal import Mutation
from cellme.vcf import INFO_FIELDS
from cellme.vcf import TrackContext
from cellme.vcf import build_header
from cellme.vcf import build_record
from cellme.vcf import build_records

CONTEXT = TrackContext(
    cell_line="MOLT4",
    sample_id="MOLT4_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE",
    study="ccle_broad_2019",
    source_build=GenomeBuild.GRCh37,
    target_build=GenomeBuild.GRCh38,
)

SAME_BUILD_CONTEXT = TrackContext(
    cell_line="MOLT4",
    sample_id="MOLT4_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE",
    study="ccle_broad_2019",
    source_build=GenomeBuild.GRCh37,
    target_build=GenomeBuild.GRCh37,
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


def add_offset(_chromosome: str, position: int) -> int:
    return position + 100000


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


def test_unliftable_mutation_is_dropped() -> None:
    assert (
        build_record(make_mutation(), CONTEXT, lift_position=never_lifts, anchor_base=None) is None
    )


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

    def lift(chromosome: str, position: int) -> int | None:
        if chromosome == "unmappable":
            return None
        return position

    records, dropped = build_records(mutations, CONTEXT, lift_position=lift, anchor_base=None)
    assert dropped == 1
    assert [record.contig for record in records] == ["1", "17"]


def test_header_declares_reference_contigs_and_full_info_schema() -> None:
    header_text = str(build_header(CONTEXT, "0.1.0"))
    assert "##reference=GRCh38" in header_text
    assert "##source=cellme 0.1.0" in header_text
    assert "##contig=<ID=17,length=83257441>" in header_text
    for field in INFO_FIELDS:
        assert f"##INFO=<ID={field.key}," in header_text
