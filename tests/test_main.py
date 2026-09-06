from inspect import isfunction
from pathlib import Path
from typing import Callable

import pysam
import pytest
from defopt import signature
from pytest import CaptureFixture
from pytest import MonkeyPatch

from cellme import __version__
from cellme import main
from cellme.cbioportal import Mutation

SAMPLE_ID = "MOLT4_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE"


def _snv_on_contig_17(reference_allele: str) -> Mutation:
    """Build a single hg19 SNV on contig 17 at position 5, for CLI reference checks."""
    return Mutation(
        gene="TP53",
        entrez_gene_id=7157,
        chromosome="17",
        start_position=5,
        end_position=5,
        reference_allele=reference_allele,
        variant_allele="T",
        protein_change="R306*",
        variant_class="Nonsense_Mutation",
        variant_type="SNP",
        ncbi_build="GRCh37",
        refseq_mrna_id="NM_001126112.2",
        protein_pos_start=306,
        protein_pos_end=306,
    )


def _write_reference(tmp_path: Path) -> Path:
    """Write and index a tiny FASTA whose contig 17 has base 'A' at 1-based position 5."""
    fasta = tmp_path / "ref.fa"
    fasta.write_text(">17\nACGTACGTACGT\n")
    pysam.faidx(str(fasta))
    return fasta


def _patch_cbioportal(monkeypatch: MonkeyPatch, mutation: Mutation) -> None:
    """Stub the cBioPortal fetches so the CLI runs offline against one sample and mutation."""
    monkeypatch.setattr("cellme.tools.truth_track.fetch_sample_ids", lambda *_a, **_k: [SAMPLE_ID])
    monkeypatch.setattr("cellme.tools.truth_track.fetch_mutations", lambda *_a, **_k: [mutation])


@pytest.mark.parametrize("tool", main._tools)
def test_tools_are_defined(tool: Callable[..., None]) -> None:
    """Test that all command line tools passed to defopt are defined functions."""
    assert isfunction(tool)


@pytest.mark.parametrize("tool", main._tools)
def test_tools_have_valid_docstrings(tool: Callable[..., None]) -> None:
    """Test that all command line tools have a valid defopt docstring."""
    try:
        signature(tool)
    except TypeError:
        raise AssertionError(f"defopt could not parse docstring for {tool.__name__}") from None


def test_cli_version_flag_prints_version_and_exits(
    monkeypatch: MonkeyPatch, capsys: CaptureFixture[str]
) -> None:
    """`cellme --version` should print the package version and exit cleanly."""
    monkeypatch.setattr("sys.argv", ["cellme", "--version"])
    with pytest.raises(SystemExit) as exc_info:
        main.run()
    assert exc_info.value.code == 0
    assert capsys.readouterr().out.strip() == __version__


def test_cli_help_lists_the_validation_flags(
    monkeypatch: MonkeyPatch, capsys: CaptureFixture[str]
) -> None:
    """`cellme --help` should advertise the liftover and reference validation flags."""
    monkeypatch.setattr("sys.argv", ["cellme", "--help"])
    with pytest.raises(SystemExit):
        main.run()
    help_text = capsys.readouterr().out
    assert "--raise-on-liftover-fails" in help_text
    assert "--skip-ref-mismatch" in help_text


def test_cli_help_lists_the_contig_style_flag(
    monkeypatch: MonkeyPatch, capsys: CaptureFixture[str]
) -> None:
    """`cellme --help` should advertise the contig-style opt-out and its choices."""
    monkeypatch.setattr("sys.argv", ["cellme", "--help"])
    with pytest.raises(SystemExit):
        main.run()
    help_text = capsys.readouterr().out
    assert "--contig-style" in help_text
    assert "{ucsc,ensembl}" in help_text


def test_cli_writes_chr_prefixed_contigs_by_default(
    monkeypatch: MonkeyPatch, tmp_path: Path
) -> None:
    """Without --contig-style the emitted CHROM is UCSC chr-prefixed."""
    reference = _write_reference(tmp_path)
    _patch_cbioportal(monkeypatch, _snv_on_contig_17("A"))
    output = tmp_path / "out.vcf"
    monkeypatch.setattr(
        "sys.argv",
        [
            "cellme",
            "MOLT4",
            "--build",
            "hg19",
            "--reference",
            str(reference),
            "--output",
            str(output),
        ],
    )
    main.run()
    with pysam.VariantFile(str(output)) as vcf:
        assert [record.contig for record in vcf] == ["chr17"]
        assert "##contig=<ID=chr17," in str(vcf.header)


def test_cli_ensembl_contig_style_writes_unprefixed_contigs(
    monkeypatch: MonkeyPatch, tmp_path: Path
) -> None:
    """--contig-style ensembl opts back into the unprefixed Ensembl names."""
    reference = _write_reference(tmp_path)
    _patch_cbioportal(monkeypatch, _snv_on_contig_17("A"))
    output = tmp_path / "out.vcf"
    monkeypatch.setattr(
        "sys.argv",
        [
            "cellme",
            "MOLT4",
            "--build",
            "hg19",
            "--reference",
            str(reference),
            "--output",
            str(output),
            "--contig-style",
            "ensembl",
        ],
    )
    main.run()
    with pysam.VariantFile(str(output)) as vcf:
        assert [record.contig for record in vcf] == ["17"]
        assert "##contig=<ID=17," in str(vcf.header)


def test_cli_aborts_on_reference_mismatch_by_default(
    monkeypatch: MonkeyPatch, capsys: CaptureFixture[str], tmp_path: Path
) -> None:
    """With --reference, a REF that disagrees with the FASTA should abort with a hint."""
    reference = _write_reference(tmp_path)
    _patch_cbioportal(monkeypatch, _snv_on_contig_17("G"))
    output = tmp_path / "out.vcf"
    monkeypatch.setattr(
        "sys.argv",
        [
            "cellme",
            "MOLT4",
            "--build",
            "hg19",
            "--reference",
            str(reference),
            "--output",
            str(output),
        ],
    )
    with pytest.raises(SystemExit) as exc_info:
        main.run()
    assert exc_info.value.code == 1
    stderr = capsys.readouterr().err
    assert "does not match the reference" in stderr
    assert "--skip-ref-mismatch" in stderr
    assert not output.exists()


def test_cli_skips_reference_mismatch_with_flag(monkeypatch: MonkeyPatch, tmp_path: Path) -> None:
    """With --skip-ref-mismatch the run completes, dropping the mismatched record."""
    reference = _write_reference(tmp_path)
    _patch_cbioportal(monkeypatch, _snv_on_contig_17("G"))
    output = tmp_path / "out.vcf"
    monkeypatch.setattr(
        "sys.argv",
        [
            "cellme",
            "MOLT4",
            "--build",
            "hg19",
            "--reference",
            str(reference),
            "--output",
            str(output),
            "--skip-ref-mismatch",
        ],
    )
    main.run()
    with pysam.VariantFile(str(output)) as vcf:
        assert list(vcf) == []


def test_cli_writes_records_when_reference_matches(
    monkeypatch: MonkeyPatch, tmp_path: Path
) -> None:
    """A REF matching the FASTA passes validation and is written to the VCF."""
    reference = _write_reference(tmp_path)
    _patch_cbioportal(monkeypatch, _snv_on_contig_17("A"))
    output = tmp_path / "out.vcf"
    monkeypatch.setattr(
        "sys.argv",
        [
            "cellme",
            "MOLT4",
            "--build",
            "hg19",
            "--reference",
            str(reference),
            "--output",
            str(output),
        ],
    )
    main.run()
    with pysam.VariantFile(str(output)) as vcf:
        assert [record.info["GENE"] for record in vcf] == ["TP53"]
