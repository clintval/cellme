"""Human genome build definitions and their reference contig tables."""

from enum import StrEnum


class GenomeBuild(StrEnum):
    """
    A human reference genome build supported by cellme.

    Only the two builds relevant to CCLE truth tracks are modeled: GRCh37, the
    build that CCLE reports its coordinates against, and GRCh38, the current
    human reference. The ``hg19`` and ``hg38`` aliases resolve to GRCh37 and
    GRCh38 respectively when the enum is constructed from a string.
    """

    GRCh38 = "GRCh38"
    GRCh37 = "GRCh37"

    @classmethod
    def _missing_(cls, value: object) -> "GenomeBuild | None":
        """
        Resolve the ``hg19`` and ``hg38`` aliases case-insensitively.

        Args:
            value: The value that did not match a canonical member name.

        Returns:
            The matching build for a known alias, otherwise None.
        """
        aliases = {"hg38": cls.GRCh38, "hg19": cls.GRCh37}
        if isinstance(value, str):
            return aliases.get(value.lower())
        return None

    @property
    def ucsc_name(self) -> str:
        """The UCSC-style assembly name used by liftover chain files."""
        return {GenomeBuild.GRCh38: "hg38", GenomeBuild.GRCh37: "hg19"}[self]


GRCH37_CONTIGS: tuple[tuple[str, int], ...] = (
    ("1", 249250621),
    ("2", 243199373),
    ("3", 198022430),
    ("4", 191154276),
    ("5", 180915260),
    ("6", 171115067),
    ("7", 159138663),
    ("8", 146364022),
    ("9", 141213431),
    ("10", 135534747),
    ("11", 135006516),
    ("12", 133851895),
    ("13", 115169878),
    ("14", 107349540),
    ("15", 102531392),
    ("16", 90354753),
    ("17", 81195210),
    ("18", 78077248),
    ("19", 59128983),
    ("20", 63025520),
    ("21", 48129895),
    ("22", 51304566),
    ("X", 155270560),
    ("Y", 59373566),
    ("MT", 16569),
)
"""Primary-assembly contig names and lengths for GRCh37 (Ensembl naming)."""

GRCH38_CONTIGS: tuple[tuple[str, int], ...] = (
    ("1", 248956422),
    ("2", 242193529),
    ("3", 198295559),
    ("4", 190214555),
    ("5", 181538259),
    ("6", 170805979),
    ("7", 159345973),
    ("8", 145138636),
    ("9", 138394717),
    ("10", 133797422),
    ("11", 135086622),
    ("12", 133275309),
    ("13", 114364328),
    ("14", 107043718),
    ("15", 101991189),
    ("16", 90338345),
    ("17", 83257441),
    ("18", 80373285),
    ("19", 58617616),
    ("20", 64444167),
    ("21", 46709983),
    ("22", 50818468),
    ("X", 156040895),
    ("Y", 57227415),
    ("MT", 16569),
)
"""Primary-assembly contig names and lengths for GRCh38 (Ensembl naming)."""

_CONTIGS_BY_BUILD: dict[GenomeBuild, tuple[tuple[str, int], ...]] = {
    GenomeBuild.GRCh37: GRCH37_CONTIGS,
    GenomeBuild.GRCh38: GRCH38_CONTIGS,
}


def contigs_for(build: GenomeBuild) -> tuple[tuple[str, int], ...]:
    """
    Return the ordered primary-assembly contigs for a build.

    Args:
        build: The genome build to look up.

    Returns:
        A tuple of (contig name, length) pairs in canonical karyotype order.
    """
    return _CONTIGS_BY_BUILD[build]


def contig_order(build: GenomeBuild) -> dict[str, int]:
    """
    Return a mapping from contig name to its sort index for a build.

    Args:
        build: The genome build to look up.

    Returns:
        A mapping from contig name to a zero-based rank in karyotype order.
    """
    return {name: index for index, (name, _length) in enumerate(contigs_for(build))}
