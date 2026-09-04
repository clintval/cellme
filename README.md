# cellme

[![CI](https://github.com/clintval/cellme/actions/workflows/python_package.yml/badge.svg?branch=main)](https://github.com/clintval/cellme/actions/workflows/python_package.yml?query=branch%3Amain)
[![Python Versions](https://img.shields.io/badge/python-3.12_|_3.13_|_3.14-blue)](https://github.com/clintval/cellme)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![MyPy Checked](http://www.mypy-lang.org/static/mypy_badge.svg)](http://mypy-lang.org/)
[![uv](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/uv/main/assets/badge/v0.json)](https://github.com/astral-sh/uv)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://docs.astral.sh/ruff/)

Convert a human cell line identifier into a truth-track VCF of its known mutations.

> [!NOTE]
> cellme currently supports human cell lines only.

## Introduction

A truth track is a curated list of the variants a sample is expected to carry.
cellme builds one for a cell line by resolving a name such as `MOLT-4` to a [Cancer Cell Line Encyclopedia](https://sites.broadinstitute.org/ccle/) sample, fetching that sample's somatic mutations from the [cBioPortal](https://www.cbioportal.org/) REST API, and writing them as a sorted VCF.
The mutation source is the `ccle_broad_2019` study, whose coordinates are reported against GRCh37.
When the requested build is GRCh38 every coordinate is lifted with a UCSC chain file using [`liftover`](https://github.com/jeremymcrae/liftover), and the VCF is written with [`pysam`](https://github.com/pysam-developers/pysam).

## Recommended Installation

Install the Python package and dependency management tool [`uv`](https://docs.astral.sh/uv/getting-started/installation/) using the official documentation.

Install the dependencies of the project with:

```console
uv sync --locked
```

To check successful installation, run:

```console
uv run cellme --help
```

## Usage

Write a GRCh38 truth track for MOLT-4 to a file:

```console
uv run cellme "MOLT-4" --build GRCh38 --output MOLT-4.GRCh38.vcf
```

The query is matched case-insensitively on the leading cell line token, so `MOLT-4`, `MOLT4`, and the full `MOLT4_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE` sample id all resolve to the same sample.
Omit `--output` to stream the VCF to standard output.
Pass `--build GRCh37` to emit the CCLE coordinates unchanged, with no liftover.

The resulting VCF carries a full header and one record per mutation, including the classic MOLT-4 markers TP53 R306* (nonsense) and PTEN K267Rfs*9 (frameshift deletion):

```vcf
##fileformat=VCFv4.2
##source=cellme 0.1.0
##reference=GRCh38
##cellme_cellLine=MOLT4
##cellme_sampleId=MOLT4_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE
##cellme_sourceStudy=cBioPortal CCLE ccle_broad_2019
##cellme_sourceBuild=GRCh37
##contig=<ID=1,length=248956422>
...
##contig=<ID=17,length=83257441>
...
##INFO=<ID=GENE,Number=1,Type=String,Description="HUGO gene symbol">
##INFO=<ID=PROTEIN_CHANGE,Number=1,Type=String,Description="Protein-level change in CCLE short form (HGVS p.-like), e.g. R306*">
##INFO=<ID=VARIANT_CLASS,Number=1,Type=String,Description="Variant classification (MAF-style), e.g. Missense_Mutation">
##INFO=<ID=VARIANT_TYPE,Number=1,Type=String,Description="Sequence alteration type reported by CCLE: SNP, DNP, INS, or DEL">
##INFO=<ID=CELL_LINE,Number=1,Type=String,Description="Cell line resolved from the query">
##INFO=<ID=SAMPLE_ID,Number=1,Type=String,Description="CCLE sample identifier for the cell line">
##INFO=<ID=ENTREZ,Number=1,Type=Integer,Description="NCBI Entrez gene identifier">
##INFO=<ID=REFSEQ,Number=1,Type=String,Description="RefSeq mRNA accession for the annotated transcript">
##INFO=<ID=PROTEIN_POS,Number=1,Type=String,Description="Affected protein position or range, 1-based">
##INFO=<ID=SOURCE,Number=1,Type=String,Description="Source database and study, e.g. cBioPortal CCLE ccle_broad_2019">
##INFO=<ID=ORIGINAL_BUILD,Number=1,Type=String,Description="Genome build of the source coordinates before any liftover">
##INFO=<ID=LIFTED,Number=0,Type=Flag,Description="Coordinate was lifted from ORIGINAL_BUILD to the output reference build">
##INFO=<ID=ANCHOR,Number=1,Type=String,Description="Provenance of the indel anchor base: reference (from --reference) or placeholder (N)">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	6197724	RPL22:K15Rfs*5	NT	N	.	.	GENE=RPL22;PROTEIN_CHANGE=K15Rfs*5;VARIANT_CLASS=Frame_Shift_Del;VARIANT_TYPE=DEL;CELL_LINE=MOLT4;SAMPLE_ID=MOLT4_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE;ENTREZ=6146;REFSEQ=NM_000983.3;PROTEIN_POS=15;SOURCE=cBioPortal CCLE ccle_broad_2019;ORIGINAL_BUILD=GRCh37;LIFTED;ANCHOR=placeholder
10	87958012	PTEN:K267Rfs*9	NA	N	.	.	GENE=PTEN;PROTEIN_CHANGE=K267Rfs*9;VARIANT_CLASS=Frame_Shift_Del;VARIANT_TYPE=DEL;CELL_LINE=MOLT4;SAMPLE_ID=MOLT4_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE;ENTREZ=5728;REFSEQ=NM_000314.4;PROTEIN_POS=265-267;SOURCE=cBioPortal CCLE ccle_broad_2019;ORIGINAL_BUILD=GRCh37;LIFTED;ANCHOR=placeholder
17	7673704	TP53:R306*	G	A	.	.	GENE=TP53;PROTEIN_CHANGE=R306*;VARIANT_CLASS=Nonsense_Mutation;VARIANT_TYPE=SNP;CELL_LINE=MOLT4;SAMPLE_ID=MOLT4_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE;ENTREZ=7157;REFSEQ=NM_001126112.2;PROTEIN_POS=306;SOURCE=cBioPortal CCLE ccle_broad_2019;ORIGINAL_BUILD=GRCh37;LIFTED
```

## VCF INFO Fields

Every record is annotated with a compact, self-describing INFO block.

| Key              | Type    | Description                                                                 |
| ---------------- | ------- | --------------------------------------------------------------------------- |
| `GENE`           | String  | HUGO gene symbol.                                                           |
| `PROTEIN_CHANGE` | String  | Protein-level change in CCLE short form (HGVS p.-like), e.g. `R306*`.        |
| `VARIANT_CLASS`  | String  | Variant classification (MAF-style), e.g. `Missense_Mutation`.               |
| `VARIANT_TYPE`   | String  | Sequence alteration type reported by CCLE: `SNP`, `DNP`, `INS`, or `DEL`.    |
| `CELL_LINE`      | String  | Cell line resolved from the query.                                          |
| `SAMPLE_ID`      | String  | CCLE sample identifier for the cell line.                                   |
| `ENTREZ`         | Integer | NCBI Entrez gene identifier.                                                |
| `REFSEQ`         | String  | RefSeq mRNA accession for the annotated transcript.                        |
| `PROTEIN_POS`    | String  | Affected protein position or range, 1-based.                                |
| `SOURCE`         | String  | Source database and study, e.g. `cBioPortal CCLE ccle_broad_2019`.          |
| `ORIGINAL_BUILD` | String  | Genome build of the source coordinates before any liftover.                 |
| `LIFTED`         | Flag    | Coordinate was lifted from `ORIGINAL_BUILD` to the output reference build.   |
| `ANCHOR`         | String  | Provenance of the indel anchor base: `reference` or `placeholder`.          |

The `ID` column is set to a stable `gene:proteinChange` token where one is available, for example `TP53:R306*`, and falls back to `gene:chrom:pos` otherwise.

## Notes and Limitations

The `ccle_broad_2019` mutation profile is emitted in full, so the truth track contains every mutation CCLE reports for the cell line rather than a hand-picked driver subset.

Insertions and deletions are written in the left-anchored VCF convention.
Supply a reference FASTA for the target build with `--reference` to fill the anchor base from the reference; without one, the anchor base is a placeholder `N` and the record is marked `ANCHOR=placeholder`.

Substitution alleles are carried through from CCLE as reported on GRCh37 and are not re-derived from the target reference after lifting.
Coordinates that fail to lift, or that lift to a different contig, are dropped and counted in a warning.

The first GRCh38 run downloads the UCSC `hg19ToHg38` chain file into `~/.liftover`, so that run needs network access to UCSC.

## Development and Testing

See the [contributing guide](./CONTRIBUTING.md) for more information.
