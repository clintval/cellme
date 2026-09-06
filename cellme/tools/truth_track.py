"""The cellme command: resolve a cell line and write its truth-track VCF."""

import logging
from pathlib import Path

import requests

from cellme import __version__
from cellme.builds import GenomeBuild
from cellme.cbioportal import DEFAULT_STUDY
from cellme.cbioportal import cell_line_name
from cellme.cbioportal import fetch_mutations
from cellme.cbioportal import fetch_sample_ids
from cellme.cbioportal import resolve_sample
from cellme.vcf import TrackContext
from cellme.vcf import build_header
from cellme.vcf import build_records
from cellme.vcf import make_anchor_base
from cellme.vcf import make_lifter
from cellme.vcf import make_reference_lookup
from cellme.vcf import write_vcf

logger = logging.getLogger("cellme")

CCLE_BUILD: GenomeBuild = GenomeBuild.hg19
"""CCLE reports its coordinates against GRCh37 (hg19), the source build for lifting."""


def truth_track(
    query: str,
    *,
    build: GenomeBuild = GenomeBuild.hg38,
    reference: Path | None = None,
    output: Path | None = None,
    study: str = DEFAULT_STUDY,
    raise_on_liftover_fails: bool = False,
    skip_ref_mismatch: bool = False,
) -> None:
    """
    Write a truth-track VCF of a human cell line's known CCLE mutations.

    The query is resolved to a Cancer Cell Line Encyclopedia sample, its somatic
    mutations are fetched from cBioPortal, and each is written as a VCF record on
    the requested genome build. CCLE coordinates are hg19 (GRCh37); when the
    target build is hg38 (GRCh38) every coordinate is lifted with a UCSC chain
    file, reverse-complementing the alleles of any variant that maps to a
    minus-strand chain block. When the output path ends in ``.gz`` the VCF is BGZF
    compressed and a tabix ``.tbi`` index is written alongside it.

    Liftover failures are lenient by default: a coordinate that cannot be lifted
    to the target build is dropped with a warning, so the truth track still builds
    from the coordinates that do lift. Pass ``--raise-on-liftover-fails`` to abort
    the run on the first coordinate that fails to lift instead. Reference
    validation is strict by default: when ``reference`` is supplied, every emitted
    record's REF allele is checked against the reference sequence and a mismatch
    aborts the run, since a mismatch signals a miscalled record. Pass
    ``--skip-ref-mismatch`` to drop such records with a warning instead.

    Args:
        query: Cell line identifier, e.g. MOLT-4, MOLT4, or a full CCLE sample id.
        build: Target genome build for the emitted VCF: hg38 or hg19 (the aliases
            GRCh38 and GRCh37 are also accepted).
        reference: Reference FASTA for the target build. It is used to place
            spec-compliant anchor bases on insertions and deletions and, when
            supplied, to validate that each record's REF allele matches the
            reference. Without it, indel anchors use a placeholder N and are
            marked ANCHOR=placeholder and no REF validation is performed.
        output: Output VCF path. Writes to standard output when omitted.
        study: cBioPortal study identifier to query.
        raise_on_liftover_fails: Abort the run on a variant that cannot be lifted
            to the target build instead of dropping it with a warning.
        skip_ref_mismatch: Drop and warn on a record whose REF allele does not
            match the reference instead of aborting the run. Has no effect unless
            ``reference`` is supplied.
    """
    with requests.Session() as session:
        session.headers["User-Agent"] = f"cellme/{__version__}"
        sample_ids = fetch_sample_ids(study, session=session)
        sample_id = resolve_sample(query, sample_ids)
        logger.info(f"Resolved {query!r} to CCLE sample {sample_id}")
        mutations = fetch_mutations(study, sample_id, session=session)
    logger.info(f"Fetched {len(mutations)} mutations for {sample_id}")

    context = TrackContext(
        cell_line=cell_line_name(sample_id),
        sample_id=sample_id,
        study=study,
        source_build=CCLE_BUILD,
        target_build=build,
    )
    lift_position = make_lifter(CCLE_BUILD, build)
    reference_lookup = make_reference_lookup(reference)
    anchor_base = make_anchor_base(reference_lookup)
    records, dropped = build_records(
        mutations,
        context,
        lift_position=lift_position,
        anchor_base=anchor_base,
        reference_lookup=reference_lookup,
        raise_on_liftover_fails=raise_on_liftover_fails,
        skip_ref_mismatch=skip_ref_mismatch,
    )
    if dropped:
        logger.warning(f"Dropped {dropped} of {len(mutations)} mutations while building {build}")
    header = build_header(context, __version__)
    write_vcf(records, header, output)
    logger.info(f"Wrote {len(records)} records for {context.cell_line} on {build}")
