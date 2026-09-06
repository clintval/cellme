# cellme

[![PyPI version](https://img.shields.io/pypi/v/cellme)](https://pypi.org/project/cellme/)
[![CI](https://github.com/clintval/cellme/actions/workflows/python_package.yml/badge.svg?branch=main)](https://github.com/clintval/cellme/actions/workflows/python_package.yml?query=branch%3Amain)
[![Python Versions](https://img.shields.io/badge/python-3.12_|_3.13_|_3.14-blue)](https://github.com/clintval/cellme)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![MyPy Checked](http://www.mypy-lang.org/static/mypy_badge.svg)](http://mypy-lang.org/)
[![uv](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/uv/main/assets/badge/v0.json)](https://github.com/astral-sh/uv)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://docs.astral.sh/ruff/)

Convert a human cell line identifier into a truth-track VCF of its known mutations.

## Introduction

cellme builds a truth/known VCF for a cell line by resolving a name such as `MOLT-4` to a [Cancer Cell Line Encyclopedia](https://sites.broadinstitute.org/ccle/) sample, fetching that sample's somatic mutations from the [cBioPortal](https://www.cbioportal.org/) REST API, and writing them as a sorted VCF.
The mutation source is the `ccle_broad_2019` study.
CCLE reports its coordinates against hg19 (GRCh37); by default cellme lifts them to hg38 (GRCh38), and `--build hg19` keeps the native coordinates.

## Installation

Install the release from PyPI with `pip`:

```console
pip install cellme
```

## Usage

Write an hg38 truth track for the cell line MOLT-4 to a file:

```console
cellme "MOLT-4" --build hg38 --reference /ref/hg38.fa --output MOLT-4.hg38.vcf.gz
```

> [!IMPORTANT]
> Always pass `--reference` with a FASTA for the build you target, as shown above.
> It ensures the correct anchor bases are set for insertions and deletions.
> Without it, indel anchor bases fall back to a placeholder `N` and those records are marked `ANCHOR=placeholder`.
> When `--reference` is given, cellme also validates that every record's REF allele matches the reference at its target position (see [Validation](#validation)).

The build keys are `hg38` and `hg19` (the Ensembl names `GRCh38` and `GRCh37` are accepted as aliases).

The query is matched case-insensitively on the leading cell line token, so `MOLT-4`, `MOLT4`, and the full `MOLT4_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE` sample ID all resolve to the same sample.

## Validation

cellme is strict by default so a truth track is never silently partial or inconsistent.
Both checks below abort the run with an actionable error, and each has an opt-out that drops the offending variant with a warning instead.

- **Liftover failures.**
  CCLE reports hg19 (GRCh37) coordinates; lifting them to the target build can fail.
  By default a failure aborts the run, naming the offending variant and its source locus, so a partial truth track never goes unnoticed.
  Pass `--skip-liftover-fails` to drop unliftable variants with a warning instead.
- **REF-base checks.**
  Whenever `--reference` is given, cellme verifies that each emitted record's REF allele matches the reference sequence at its (post-liftover) target position: a single base for substitutions, and the anchor plus deleted bases for deletions.
  By default a mismatch aborts the run, naming the variant, the expected REF, and the actual reference base(s).
  Pass `--skip-ref-mismatch` to drop mismatched variants with a warning instead.

## VCF INFO Fields

Every record is annotated with compact, self-describing INFO fields.

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
| `ORIGINAL_LOCUS` | String  | Original locus before liftover in UCSC position format, e.g. `chr17:7577022`. |
| `LIFTED`         | Flag    | Coordinate was lifted from `ORIGINAL_BUILD` to the output reference build.   |
| `ANCHOR`         | String  | Provenance of the indel anchor base: `reference` or `placeholder`.          |

The `ID` column is set to a stable `gene:proteinChange` token where one is available, for example `TP53:R306*`, and falls back to `gene:chrom:pos` otherwise.

## Development and Testing

See the [contributing guide](./CONTRIBUTING.md) for more information.
