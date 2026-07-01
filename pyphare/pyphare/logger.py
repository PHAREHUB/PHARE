import os
import logging

LOG_LEVEL = os.environ.get("PYPHARE_LOG_LEVEL", None)
log_levels = {"INFO": 20, "WARNING": 30, "ERROR": 40, "DEBUG": 10}  # maps enum values

if LOG_LEVEL:
    logging.basicConfig(
        format="%(asctime)s,%(msecs)03d %(levelname)-8s [%(filename)s:%(lineno)d] %(message)s",
        datefmt="%Y-%m-%dT%H:%M:%S",
        level=LOG_LEVEL,
    )


def getLogger(name, level=LOG_LEVEL):
    logger = logging.getLogger(__name__)
    if LOG_LEVEL:
        level = log_levels[level] if isinstance(level, str) else level
        logger.setLevel(level)
    return logger
