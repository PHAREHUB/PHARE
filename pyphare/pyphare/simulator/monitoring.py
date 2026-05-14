#
# Resource monitoring requires phlop
#  python3 -m pip install phlop
#
#

import os
from pathlib import Path
from dataclasses import dataclass

from pyphare import cpp
from pyphare.logger import getLogger

logger = getLogger(__name__)


@dataclass
class MonitoringOptions:
    interval: int = 0
    rank_modulo: int = 1  # all ranks

    def __post_init__(self):
        if self.rank_modulo == 0:
            raise ValueError("MonitoringOptions rank_modulo cannot be zero")

    def active(self):
        return self.interval > 0 and cpp.mpi_rank() % self.rank_modulo == 0

    @staticmethod
    def FROM(input):
        if type(input) is MonitoringOptions:
            return input
        if type(input) is int:
            return MonitoringOptions(interval=input)
        if type(input) is bool:
            return MonitoringOptions(interval=100)  # seconds
        raise ValueError("MonitoringOptions not constructible from ", input)


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


def monitoring_yaml_file():
    path = Path(".phare") / "stats" / f"rank.{cpp.mpi_rank()}.yaml"
    path.parent.mkdir(exist_ok=True, parents=True)
    return path


def setup_monitoring(input):
    if not have_phlop():
        return

    from phlop.app import stats_man as sm  # pylint: disable=import-error

    options = MonitoringOptions.FROM(input)
    if not options.active():
        return

    _globals.stats_man = sm.AttachableRuntimeStatsManager(
        valdict(yaml=monitoring_yaml_file(), interval=options.interval),
        dict(rank=cpp.mpi_rank()),
    ).start()


def monitoring_shutdown():
    if not have_phlop():
        return

    cpp.mpi_barrier()  # force similar end time
    if _globals.stats_man:
        _globals.stats_man.kill().join()


def timing_setup():
    if cpp.mpi_rank() == 0:
        try:
            Path(".phare/timings").mkdir(parents=True, exist_ok=True)
        except FileNotFoundError:
            logger.error(f"Couldn't find timing dir from {os.getcwd() }")
