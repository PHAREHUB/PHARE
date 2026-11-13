#!/usr/bin/env python3

import os
import re
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


def check_bundled_libraries():
    import pybindlibs

    lib_dir = Path(pybindlibs.__file__).resolve().parent / "lib"
    if not lib_dir.is_dir() or not sys.platform.startswith("linux"):
        return
    bundled = {path.name for path in lib_dir.iterdir()}
    loaded = set()
    for line in Path("/proc/self/maps").read_text().splitlines():
        fields = line.split(maxsplit=5)
        if len(fields) < 6:
            continue
        path = Path(fields[5].strip())
        if re.match(r"lib(phare|SAMRAI)_\w+\.so", path.name) and path.name in bundled:
            loaded.add(path.resolve())
    assert loaded, f"no library of {lib_dir} is loaded"
    for path in loaded:
        assert path.parent == lib_dir, f"{path} loaded instead of the one in {lib_dir}"


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
    check_bundled_libraries()
    if rank == 0:
        print("pip install smoke: OK")


if __name__ == "__main__":
    main()
