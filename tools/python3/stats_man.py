import os
import signal
import sys
import threading
import time
from datetime import datetime
from multiprocessing import Process
from pathlib import Path

import psutil
import yaml

SCRIPT_DIRECTORY = os.path.dirname(os.path.realpath(__file__))
ROOT_DIRECTORY = Path(SCRIPT_DIRECTORY).parent.parent

config = dict(
    dir=ROOT_DIRECTORY / ".runtime_log",
    interval=2,  # seconds
)

_globals = dict()

def check_pid(pid):
    """Check For the existence of a unix pid."""
    try:
        os.kill(pid, 0)
    except OSError:
        return False
    else:
        return True


def cleanup():
    output = config["dir"]
    output.mkdir(parents=True, exist_ok=True)

    with open(output / f"stats.{_globals['pid']}.yaml", "w") as file:
        file.write(yaml.dump(_globals["captures"]))


def signal_handler(sig, frame):
    _globals["sub_run"] = False
    cleanup()
    _globals.clear()
    sys.exit(0)


def capture_now(pid):
    now = datetime.utcnow().timestamp()
    process = psutil.Process(pid=pid)
    mem = int(process.memory_info().rss / 1024**2)  # bytes -> MB # drop decimals
    fds = len(process.open_files())
    return dict(
        mem=mem,
        fds=fds,
        now=now,
    )


class RuntimeStatsManager:
    def __init__(self, pid):
        self.pid = pid
        self.p = Process(target=RuntimeStatsManager._run, args=(pid, os.getpid()))
        self.t = threading.Thread(target=RuntimeStatsManager._thread, args=(self,))
        self.p.start()
        self.t.start()

    def __del__(self):
        pid = self.p.pid
        os.kill(pid, signal.SIGINT)

        # hangs if exitcode is not checked?
        while self.p.exitcode is None and check_pid(pid):
            time.sleep(1)

    @staticmethod
    def _thread(self):
        while threading.main_thread().is_alive():
            time.sleep(1)
        _globals.clear()

    @staticmethod
    def _run(pid, host_pid):
        _globals["pid"] = pid
        _globals["host_pid"] = host_pid
        _globals["captures"] = []
        _globals["sub_run"] = True
        signal.signal(signal.SIGINT, signal_handler)
        while _globals["sub_run"]:
            time.sleep(config["interval"])
            _globals["captures"] += [capture_now(pid)]


def attach_to_this_process():
    _globals["this_proc"] = RuntimeStatsManager(os.getpid())
