#
#
#

import os
from pathlib import Path

from pyphare.logger import getLogger

logger = getLogger(__name__)


def have_phlop():
    from importlib.util import find_spec

    try:
        return find_spec("phlop.dict") is not None
    except (ImportError, ModuleNotFoundError):
        return False


def valdict(**kwargs):
    if not have_phlop():
        return dict

    from phlop.dict import ValDict  # pylint: disable=import-error

    return ValDict(**kwargs)


_globals = valdict(stats_man=None)


def monitoring_yaml_file(cpplib):
    path = Path(".phare") / "stats" / f"rank.{cpplib.mpi_rank()}.yaml"
    path.parent.mkdir(exist_ok=True, parents=True)
    return path


def setup_monitoring(cpplib, interval=10):
    if not have_phlop():
        return

    from phlop.app import stats_man as sm  # pylint: disable=import-error

    _globals.stats_man = sm.AttachableRuntimeStatsManager(
        valdict(yaml=monitoring_yaml_file(cpplib), interval=interval),
        dict(rank=cpplib.mpi_rank()),
    ).start()


def monitoring_shutdown(cpplib):
    if not have_phlop():
        return

    if _globals.stats_man:
        _globals.stats_man.kill().join()


def timing_setup(cpplib):
    if cpplib.mpi_rank() == 0:
        try:
            Path(".phare/timings").mkdir(parents=True, exist_ok=True)
        except FileNotFoundError:
            logger.error(f"Couldn't find timing dir from {os.getcwd() }")
