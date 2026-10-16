
=====================
Run PHARE from python
=====================

PHARE can run as a python script.

Installing PHARE
----------------

PHARE requires a minimmum version of python 3.11 to run properly.
Make sure `python3` shows the version is at least 3.11.
From the PHARE source directory, `pip install .` builds PHARE with CMake and installs the `pyphare`
python package together with the compiled simulator modules, in the active (virtual) environment:

.. code-block:: bash

     python3 -m venv phare_venv
     source phare_venv/bin/activate
     pip install .

or with `uv <https://docs.astral.sh/uv/>`_, which also creates the environment:

.. code-block:: bash

     uv sync
     source .venv/bin/activate

Python package dependencies are listed in `pyproject.toml`.
By default SAMRAI is built from `subprojects/` and every permutation of `res/sim/all.txt` is compiled, which takes a long time.
To use a pre-installed SAMRAI, or to build only some permutations:

.. code-block:: bash

     pip install . -C cmake.define.SAMRAI_ROOT=/path/to/SAMRAI/install \
                   -C cmake.define.PHARE_PERMUTATIONS=/path/to/permutations.txt

With uv, pass the definitions through the `SKBUILD_CMAKE_DEFINE` environment variable (semicolon-separated `KEY=VALUE`):

.. code-block:: bash

     SKBUILD_CMAKE_DEFINE="SAMRAI_ROOT=/path/to/SAMRAI/install;PHARE_PERMUTATIONS=/path/to/permutations.txt" uv sync

`SAMRAI_ROOT` and `PHARE_PERMUTATIONS` are not remembered between installs: pass them each time.
uv does not detect changes to these options, nor to a permutation file outside the repository.
After changing them, run `uv sync --reinstall-package pyphare` (pip reinstalls on every `pip install .`).

The number of simulator modules compiled concurrently is limited according to the available memory.
Override it with `-C cmake.define.PHARE_MODULE_COMPILE_JOBS=<n>` (`0` disables the limit) or
`-C cmake.define.PHARE_MODULE_COMPILE_MEMORY_MB=<mb>` (memory assumed per compile job, default 2048).

The build directory `_skbuild/` is kept between installs, so rebuilds are incremental.
Do not run `ninja` or `cmake --build` by hand inside `_skbuild/`: after a CMake file changes,
the automatic reconfigure does not have your `-C` options and falls back to the defaults.
Reinstall with pip or uv instead.

Building PHARE with your own CMake build
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Developers who keep their own CMake build should install only the python dependencies, since a plain `uv sync`
builds PHARE with every permutation:

.. code-block:: bash

     uv sync --no-install-project
     source .venv/bin/activate

or, with pip, `pip3 install -r requirements.txt` in a virtual environment.

Do not mix your own CMake build with an environment where `uv sync` or `uv run` installed PHARE.
That install is editable, and its import finder takes precedence over `PYTHONPATH`:
the installed `pybindlibs` is used instead of your build (a non-editable `pip install .` does not do this).
`uv sync` and `uv run` also manage `.venv` as the project environment, and may recreate it
(for example with a uv-managed python), dropping packages installed there by hand.
Use a separate environment for your own build, or keep uv from touching it:
`uv sync --no-install-project`, `uv run --no-sync <cmd>`, or the environment variable `UV_NO_SYNC=1`
for scripts that call `uv run`.
A `.venv` shared between several checkouts (for example symlinked into git worktrees) is synced by whichever checkout runs uv.

Installing on a cluster / overriding defaults
---------------------------------------------

PHARE is always compiled locally against your compilers, MPI and parallel HDF5: there are no prebuilt wheels.
Load the same modules before installing and before running.
Choose compilers with the `CC`, `CXX` and `FC` environment variables, CMake options with
`-C cmake.define.<VAR>=...` (or the `CMAKE_ARGS` environment variable) and the build type with `-C cmake.build-type=...`.
For example:

.. code-block:: bash

     module load cray-python gcc-native cray-hdf5-parallel
     FC=mpif90 pip install . \
         -C cmake.build-type=Release \
         -C cmake.define.CMAKE_CXX_FLAGS="-DNDEBUG -g0 -O3 -march=native -mtune=native"

- The persistent CMake cache keeps compilers, flags and modules: delete `_skbuild/` after changing them
  (only `SAMRAI_ROOT` and `PHARE_PERMUTATIONS` are forgotten automatically).
- `-march=native` binaries only run on the CPU type they were built on: build on the compute-node architecture.
- At run time, launch with the `mpirun`/`srun` of the MPI PHARE was built with. The launcher of a different MPI
  starts N independent single-rank runs, and a different MPI on `LD_LIBRARY_PATH` takes precedence over the recorded RUNPATH.


Running PHARE
-------------
After `pip install .` or `uv sync`, PHARE is accessible to python: there is nothing to set.

With your own CMake build only, let python know where to find PHARE.
Assuming PHARE source directory is in  `/path/to/PHARE`, and the build directory is
`/path/to/build/`:


.. code-block:: bash

     export PYTHONPATH=/path/to/PHARE/pyphare:/path/to/build:$PYTHONPATH

Set both paths together, or not at all.

Write a :doc:`simulation_inputs` and run the following command:


.. code-block:: bash

      python3 /path/to/my_script.py


Executable PHARE initialization script
--------------------------------------

Running PHARE from python basically consists in declaring desired blocks
detailed in a :doc:`simulation_inputs` file, and adding a `main` function to run
the simulation.

If running the simulation from python, you need to define a `main` function that
basically runs the simulation. Here is a small example:


.. code-block:: python

    from pyphare.simulator.simulator import Simulator
    from pyphare.pharein import Simulation
    from pyphare.pharein import MaxwellianFluidModel


    # define your simulation parameters
    # and initial condition...

    sim = Simulation(
    # ...
    )

    MaxwellianFluidModel(
    # ...
    )


    def main():
        simulator = Simulator(sim)
        simulator.run()



The Simulator
-------------

The simulator is used to initialize and run the simulation.
It is given the `Simulation` object as an argument, as in the 
simple example above. The simplest way to start the simulation is
then to call the `run()` method.


.. autoclass:: pyphare.simulator.simulator.Simulator
   :members:

