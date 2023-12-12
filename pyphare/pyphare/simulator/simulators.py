"""



"""
import os

import sys
import json
import time
import inspect
import importlib
import concurrent

from enum import IntEnum
from pathlib import Path
from datetime import datetime
from multiprocessing import Process
from multiprocessing import shared_memory, Lock

import numpy as np
import pyphare.pharein as ph
from pyphare.pharesee.run import Run
from pyphare.pharesee.hierarchy import hierarchy_compare
from pyphare.simulator.simulator import Simulator


_globals = dict(servers=dict(), busy=0)
phare_runs_dir = Path(os.getcwd()) / "phare_runs"  # can't contain period "."

lock = Lock()

shared_size = (200, 200)  # kinda arbitrary but over allocated


class SharedSimulationStateEnum(IntEnum):
    CURRENT_TIME = 0
    IS_BUSY = 7
    CAN_ADVANCE = 8
    IS_FINISHED = 9


def create_shared_block(n_sims):
    a = np.zeros(shared_size, dtype=np.int64)
    shm = shared_memory.SharedMemory(create=True, size=a.nbytes)
    np_array = np.ndarray(a.shape, dtype=np.int64, buffer=shm.buf)
    return shm, np_array


def atomic_set(sim_id, shr_name, pos, val):
    existing_shm = shared_memory.SharedMemory(name=shr_name)
    np_array = np.ndarray(shared_size, dtype=np.int64, buffer=existing_shm.buf)
    lock.acquire()

    np_array[sim_id][pos] = val

    lock.release()
    existing_shm.close()


def poll(sim_id, shr_name):
    while True:
        wait_time = state_machine(sim_id, shr_name)
        time.sleep(wait_time)


def state_machine(sim_id, shr_name):
    result = 2
    finished = is_finished()

    existing_shm = shared_memory.SharedMemory(name=shr_name)
    np_array = np.ndarray(shared_size, dtype=np.int64, buffer=existing_shm.buf)

    lock.acquire()  # atomic operations below
    should_advance = np_array[sim_id][SharedSimulationStateEnum.CAN_ADVANCE] == 1
    if should_advance:
        np_array[sim_id][SharedSimulationStateEnum.IS_BUSY] = 1
        np_array[sim_id][SharedSimulationStateEnum.CAN_ADVANCE] = 0
    lock.release()  # atomic operations above

    if not finished and should_advance:
        advance_sim()

    lock.acquire()  # atomic operations below
    if should_advance:
        np_array[sim_id][SharedSimulationStateEnum.IS_BUSY] = 0
    np_array[sim_id][SharedSimulationStateEnum.IS_FINISHED] = is_finished()
    np_array[sim_id][SharedSimulationStateEnum.CURRENT_TIME] = 1e9 * current_time()
    lock.release()  # atomic operations above

    existing_shm.close()

    if should_advance and finished:
        print("FINISHED", sim_id)
        exit(0)

    return result


def set_env(dic):
    for k, v in dic.items():
        os.environ[k] = v


def prepend_python_path(val):
    sys.path = [val] + sys.path


def build_diag_dir(sim_id):
    return str(phare_runs_dir / f"diags_{os.environ['PHARE_LOG_TIME']}-ID={sim_id}")


def init_sim(sim):
    ph.global_vars.sim = sim

    sim.diag_options["options"]["dir"] = build_diag_dir(os.environ["SIM_ID"])
    Path(sim.diag_options["options"]["dir"]).mkdir(parents=True, exist_ok=True)
    _globals["simulator"] = Simulator(sim).initialize()
    return str(_globals["simulator"].currentTime())


def advance_sim():
    _globals["busy"] = 1
    _globals["simulator"].advance()
    _globals["busy"] = 0


def is_finished():
    return _globals["simulator"].finished()


def current_time():
    return _globals["simulator"].currentTime()


def init_simulation(sim_id, sim, shr_name, dic):
    set_env(dic)

    init_sim(sim)

    poll(sim_id, shr_name)


def start_server_process(sim_id, sim, shr_name, dic):
    servers = _globals["servers"]
    assert sim_id not in servers

    if "build_dir" in dic:
        prepend_python_path(dic["build_dir"])

    _globals["servers"][sim_id] = Process(
        target=init_simulation, args=(sim_id, sim, shr_name, dic)
    )
    _globals["servers"][sim_id].start()

    try_count = 5
    for i in range(1, try_count + 1):
        time.sleep(0.5)
        try:
            assert servers[sim_id].exitcode is None  # or it crashed/exited early
            return
        except Exception as e:
            if i == try_count:
                raise e


def stop_servers():
    for k, v in _globals["servers"].items():
        v.terminate()


def build_dir_path(path):
    p = Path(os.path.realpath(path))
    if not p.exists():
        p = Path(os.getcwd()) / path
    assert p.exists()
    return str(p)


class Simulators:
    def __init__(self, starting_sim_id=10):
        self.simulations = []
        self.simulation_configs = []
        self.states = dict(init=False)
        self.starting_sim_id = starting_sim_id
        self.log_time = datetime.now().strftime("%m_%d_%Y_%H_%M_%S")
        os.environ["PHARE_LOG_TIME"] = self.log_time

        # loaded during init
        self.thread_pool = None
        self.shr = None
        self.shared_np_array = None

    def register(self, simulation, build_dir=None, diag_dir=None):
        self.simulations += [simulation]
        self.simulation_configs += [dict(build_dir=build_dir, diag_dir=diag_dir)]
        ph.global_vars.sim = None

    def init(self):
        shr, np_array = create_shared_block(len(self.simulations))
        self.shr = shr
        self.shared_np_array = np_array
        self._state_machine_set_per_simulation(SharedSimulationStateEnum.IS_BUSY, 1)

        self.thread_pool = concurrent.futures.ThreadPoolExecutor(
            max_workers=len(self.simulations)
        )

        for i, simulation in enumerate(self.simulations):
            sim_id = self.starting_sim_id + i
            simulation_config = self.simulation_configs[i]

            init_dict = dict(SIM_ID=str(sim_id), PHARE_LOG_TIME=self.log_time)
            if simulation_config["build_dir"]:
                init_dict["build_dir"] = build_dir_path(simulation_config["build_dir"])

            start_server_process(
                sim_id,
                simulation,
                shr.name,
                init_dict,
            )
        self.states["init"] = True

    def _state_machine_list_for_value(self, offset):
        existing_shm = shared_memory.SharedMemory(name=self.shr.name)
        np_array = np.ndarray(shared_size, dtype=np.int64, buffer=existing_shm.buf)
        lock.acquire()

        results = np.ndarray(len(self.simulations))
        for i, simulation in enumerate(self.simulations):
            sim_id = self.starting_sim_id + i
            results[i] = np_array[sim_id][offset]

        lock.release()
        existing_shm.close()
        return results

    def _state_machine_set_per_simulation(self, offset, val):
        existing_shm = shared_memory.SharedMemory(name=self.shr.name)
        np_array = np.ndarray(shared_size, dtype=np.int64, buffer=existing_shm.buf)
        lock.acquire()

        for i, simulation in enumerate(self.simulations):
            sim_id = self.starting_sim_id + i
            np_array[sim_id][offset] = val

        lock.release()
        existing_shm.close()

    def wait_for_simulations(self):
        # we have to wait for all simulations to finish the current timestep
        while True:
            if (
                np.sum(
                    self._state_machine_list_for_value(
                        SharedSimulationStateEnum.IS_BUSY
                    )
                )
                == 0
            ):
                return
            time.sleep(1)

    def advance(self, compare=False):
        if not self.states["init"]:
            self.init()

        atomic_set(0, self.shr.name, 0, int(compare))
        self._state_machine_set_per_simulation(SharedSimulationStateEnum.CAN_ADVANCE, 1)

        if compare:
            self._state_machine_set_per_simulation(SharedSimulationStateEnum.IS_BUSY, 1)
            # we have to wait for all simulations to finish the current timestep
            self.wait_for_simulations()
            self.compare()

    def compare(self):
        for i, simulation in enumerate(self.simulations):
            assert (
                _globals["servers"][self.starting_sim_id + i].exitcode is None
            ), "or it crashed early"

        sim_times = self._state_machine_list_for_value(
            SharedSimulationStateEnum.CURRENT_TIME
        )

        times = {
            i: (0.0 + sim_times[i]) / 1e9  # :( it's an int so after decimal is dropped
            for i, simulation in enumerate(self.simulations)
        }

        if len(times) < 2:
            return

        for i, simulation in enumerate(self.simulations):
            for j in range(i + 1, len(times)):
                if times[i] == times[j]:
                    run0 = Run(
                        build_diag_dir(self.starting_sim_id + i),
                        single_hier_for_all_quantities=True,
                    )
                    run1 = Run(
                        build_diag_dir(self.starting_sim_id + j),
                        single_hier_for_all_quantities=True,
                    )

                    print(
                        f"comparing {self.starting_sim_id + i} and {self.starting_sim_id + j} at time {times[i]}"
                    )
                    assert hierarchy_compare(
                        run0.GetAllAvailableQties(times[i], []),
                        run1.GetAllAvailableQties(times[i], []),
                    )
                    print("OK!")

    def run(self, compare=False):
        if not self.states["init"]:
            self.init()

        while self.simulations:
            self.advance(compare=compare)

            is_busy = self._state_machine_list_for_value(
                SharedSimulationStateEnum.IS_BUSY
            )
            is_finished = self._state_machine_list_for_value(
                SharedSimulationStateEnum.IS_FINISHED
            )

            # trim finished simulations
            self.simulations = [
                sim
                for i, sim in enumerate(self.simulations)
                if is_busy[i] or not is_finished[i]
            ]

            print("running simulations", len(self.simulations))
            time.sleep(1)

        self.kill()

    def __del__(self):
        self.kill()

    def kill(self):
        time.sleep(2)

        if self.shr:
            self.shr.close()
            self.shr.unlink()
            self.shr = None

        stop_servers()
