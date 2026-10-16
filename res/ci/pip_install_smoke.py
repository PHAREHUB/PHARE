#!/usr/bin/env python3

import os
import sys
from pathlib import Path

assert "PYTHONPATH" not in os.environ

import pyphare
import pyphare.pharein as ph
from pyphare import cpp
from pyphare.cpp import validate
from pyphare.simulator.simulator import Simulator, startMPI

ph.NO_GUI()

REPOSITORY = Path(__file__).resolve().parents[2]


def outside_repository(path):
    return not Path(path).resolve().is_relative_to(REPOSITORY)


def check_install(editable):
    import pybindlibs

    for path in pybindlibs.__path__:
        assert outside_repository(path), path
    if not editable:
        assert outside_repository(pyphare.__file__), pyphare.__file__
    assert "PHARE_CONFIG_ERROR" not in cpp.build_config()
    validate.check_build_config_is_runtime_compatible()


def config():
    sim = ph.Simulation(
        interp_order=1,
        refined_particle_nbr=2,
        cells=20,
        dl=0.2,
        time_step_nbr=5,
        final_time=0.005,
    )

    ph.MaxwellianFluidModel(
        bx=lambda x: 1.0,
        protons={
            "charge": 1,
            "density": lambda x: 1.0,
            "vbulkx": lambda x: 0.0,
            "vbulky": lambda x: 0.0,
            "vbulkz": lambda x: 0.0,
            "vthx": lambda x: 1.0,
            "vthy": lambda x: 1.0,
            "vthz": lambda x: 1.0,
            "nbr_part_per_cell": 50,
        },
    )
    ph.ElectronModel(closure="isothermal", Te=0.12)

    return sim


def main():
    check_install("--editable" in sys.argv[1:])
    startMPI()
    rank = cpp.mpi_rank()
    Simulator(config()).run()
    if rank == 0:
        print("pip install smoke: OK")


if __name__ == "__main__":
    main()
