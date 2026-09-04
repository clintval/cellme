import logging
import sys
from collections.abc import Callable

import defopt

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
    defopt.run(command, argv=sys.argv[1:], version=True)
    logger.info("Finished executing successfully.")
