import os
import logging

logging.basicConfig(level=logging.DEBUG)

LOG_LEVEL = os.environ.get("PYPHARE_LOG_LEVEL", "INFO")
log_levels = {"INFO": 20, "ERROR": 40, "DEBUG": 10}  # maps enum values


def getLogger(name, level=LOG_LEVEL):
    logger = logging.getLogger(__name__)
    level = log_levels[level] if isinstance(level, str) else level
    logger.setLevel(level)
    return logger
