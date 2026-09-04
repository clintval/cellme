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
from cellme.vcf import write_vcf

logger = logging.getLogger("cellme")

CCLE_BUILD: GenomeBuild = GenomeBuild.GRCh37
"""CCLE reports its coordinates against GRCh37, the source build for lifting."""


def truth_track(
    query: str,
    *,
    build: GenomeBuild = GenomeBuild.GRCh38,
    reference: Path | None = None,
    output: Path | None = None,
    study: str = DEFAULT_STUDY,
) -> None:
    """
    Write a truth-track VCF of a human cell line's known CCLE mutations.

    The query is resolved to a Cancer Cell Line Encyclopedia sample, its somatic
    mutations are fetched from cBioPortal, and each is written as a VCF record on
    the requested genome build. CCLE coordinates are GRCh37; when the target
    build is GRCh38 every coordinate is lifted with a UCSC chain file.

    Args:
        query: Cell line identifier, e.g. MOLT-4, MOLT4, or a full CCLE sample id.
        build: Target genome build for the emitted VCF.
        reference: Reference FASTA for the target build, used only to place
            spec-compliant anchor bases on insertions and deletions. Without it,
            indel anchors use a placeholder N and are marked ANCHOR=placeholder.
        output: Output VCF path. Writes to standard output when omitted.
        study: cBioPortal study identifier to query.
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
    anchor_base = make_anchor_base(reference)
    records, dropped = build_records(
        mutations, context, lift_position=lift_position, anchor_base=anchor_base
    )
    if dropped:
        logger.warning(f"Dropped {dropped} mutations that could not be lifted to {build}")
    header = build_header(context, __version__)
    write_vcf(records, header, output)
    logger.info(f"Wrote {len(records)} records for {context.cell_line} on {build}")
