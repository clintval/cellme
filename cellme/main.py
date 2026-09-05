import logging
import sys
from collections.abc import Callable

import defopt

from cellme.builds import GenomeBuild
from cellme.cbioportal import CellLineError
from cellme.tools.truth_track import truth_track

_tools: list[Callable[..., None]] = [
    truth_track,
]


def setup_logging(level: str = "INFO") -> None:
    """Set up basic logging to print to the console."""
    logging.basicConfig(
        level=level,
        format="%(asctime)s %(name)s:%(funcName)s:%(lineno)s [%(levelname)s]: %(message)s",
    )


def run() -> None:
    """Set up logging, then hand over to defopt for running the command line tool."""
    setup_logging()
    logger = logging.getLogger("cellme")
    logger.info("Executing: " + " ".join(sys.argv))
    (command,) = _tools
    try:
        # Parse --build through GenomeBuild so its hidden GRCh aliases are accepted
        # while the help text still advertises only the canonical hg38 / hg19 keys.
        defopt.run(command, argv=sys.argv[1:], version=True, parsers={GenomeBuild: GenomeBuild})
    except CellLineError as error:
        print(f"Error: {error}", file=sys.stderr)
        raise SystemExit(1) from error
    logger.info("Finished executing successfully.")
