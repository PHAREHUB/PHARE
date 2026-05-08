
=====================
Run PHARE from Python
=====================

PHARE runs as a Python script. There is no standalone binary — all simulations are
configured and launched from Python.


Python dependencies
-------------------

PHARE requires a minimum version of Python 3.8.
Make sure ``python3`` shows the version is at least 3.8.
Python package dependencies are listed in the ``requirements.txt`` file.

Install dependencies for the user:

.. code-block:: bash

    pip3 install --user -r requirements.txt

Install dependencies in a virtual environment:

.. code-block:: bash

    python3 -m venv phare_venv
    source phare_venv/bin/activate
    pip3 install -r requirements.txt


Running PHARE
-------------

First, make sure PHARE is accessible to Python.
Assuming the PHARE source directory is at ``/path/to/PHARE`` and the build directory is
``/path/to/build``, set the Python path as follows:

.. code-block:: bash

    export PYTHONPATH=/path/to/PHARE/pyphare:/path/to/build:$PYTHONPATH

Write a :doc:`simulation_inputs` script, then run it:

.. code-block:: bash

    python3 /path/to/my_script.py

**MPI execution**

To run with multiple MPI processes, use ``mpirun`` (or ``mpiexec``):

.. code-block:: bash

    mpirun -np 4 python3 -Ou /path/to/my_script.py

The ``-Ou`` flags are recommended for MPI runs:

- ``-O`` — basic Python optimizations (removes ``assert`` statements and ``__debug__`` guards).
- ``-u`` — forces unbuffered stdout and stderr, so output from all ranks is written
  immediately rather than being held in per-rank buffers. This produces cleaner,
  interleaved output and avoids missing log lines if the process is killed early.

.. note::

    Only MPI rank 0 prints progress lines during ``simulator.run()``. All C++ log output
    is redirected to per-rank ``.log/`` files by default (see `Log files`_ below).


Executable PHARE initialization script
---------------------------------------

A PHARE script has two parts:

1. **Configuration** — declare physics objects (``Simulation``, ``MaxwellianFluidModel``,
   diagnostics, etc.) at module scope. See :doc:`simulation_inputs` for the full API.
2. **Execution** — a ``main()`` function that creates a ``Simulator`` and runs it.

.. code-block:: python

    import pyphare.pharein as ph
    from pyphare.pharein import MaxwellianFluidModel
    from pyphare.simulator.simulator import Simulator

    sim = ph.Simulation(
        time_step_nbr=1000,
        final_time=10.0,
        cells=(100,),
        dl=(0.2,),
        # ... other parameters
    )

    MaxwellianFluidModel(
        # ... species parameters
    )


    def main():
        simulator = Simulator(sim)
        simulator.initialize()
        simulator.run()


    if __name__ == "__main__":
        main()

.. note::

    Always guard the ``main()`` call with ``if __name__ == "__main__":`` so that the
    module-level configuration objects are created when the file is imported, but the
    simulation only runs when the script is executed directly.


The Simulator lifecycle
-----------------------

The ``Simulator`` class manages the full simulation lifecycle. Understanding its phases
allows fine-grained control, for example to inject analysis between time steps.

**Phase 1 — Construction**: ``Simulator(sim)``

Binds the simulator to a ``Simulation`` configuration object. No C++ or SAMRAI
resources are allocated yet.

.. code-block:: python

    simulator = Simulator(sim)

**Phase 2 — Initialization**: ``simulator.initialize()``

SAMRAI creates the patch hierarchy, populates initial conditions from Python
initializers, and performs the first diagnostic dump if one is scheduled at ``t=0``.
Calling ``initialize()`` is optional before ``run()`` or ``advance()`` — both methods
call it automatically if it has not been called yet.

.. code-block:: python

    simulator.initialize()

**Phase 3 — Advancing**: ``simulator.advance(dt)`` or ``simulator.run()``

- ``advance(dt)`` advances the simulation by one time step of size ``dt``. If ``dt``
  is omitted it uses the simulation's configured time step. Diagnostics are dumped
  automatically after each step when ``auto_dump=True`` (the default).
- ``run()`` loops ``advance()`` from the current time to the configured end time,
  printing per-step timing to rank 0, and calls ``reset()`` on completion.

.. code-block:: python

    # Manual step-by-step loop
    while simulator.currentTime() < sim.final_time:
        simulator.advance()

    # Or simply:
    simulator.run()

**Phase 4 — Reset**: ``simulator.reset()``

Releases all C++ and SAMRAI resources. Called automatically by ``run()`` and at
process exit. Call it explicitly if you need to free resources before the script ends.

.. code-block:: python

    simulator.reset()

**Complete example**

.. code-block:: python

    import pyphare.pharein as ph
    from pyphare.simulator.simulator import Simulator

    sim = ph.Simulation(
        time_step_nbr=500,
        final_time=5.0,
        cells=(128,),
        dl=(0.1,),
    )

    # ... define models and diagnostics ...


    def main():
        simulator = Simulator(sim)
        simulator.initialize()

        # custom step-by-step loop
        while simulator.currentTime() < sim.final_time:
            simulator.advance()

        simulator.reset()


    if __name__ == "__main__":
        main()


Post-advance callbacks
----------------------

Pass a callable as ``post_advance`` to execute Python code after every time step.
The callback receives the current simulation time as its only argument. It is invoked
after diagnostics are dumped for that step.

.. code-block:: python

    def my_callback(time):
        print(f"Advanced to t={time:.4f}")

    simulator = Simulator(sim, post_advance=my_callback)
    simulator.run()

This is useful for lightweight online analysis or progress logging without modifying
the main run loop.

.. warning::

    The callback runs on every rank. If you perform I/O inside it, guard the write
    with an MPI rank check to avoid all ranks writing the same data:

    .. code-block:: python

        from pyphare import cpp

        def my_callback(time):
            if cpp.mpi_rank() == 0:
                print(f"t = {time:.4f}")

        simulator = Simulator(sim, post_advance=my_callback)
        simulator.run()


Dry run mode
------------

A dry run validates the simulation configuration — parsing inputs and setting up C++
data structures — without actually advancing the simulation. This is useful for
catching configuration errors quickly.

Enable via the ``Simulation`` constructor:

.. code-block:: python

    sim = ph.Simulation(
        # ...
        dry_run=True,
    )

Or via the environment variable, without changing the script:

.. code-block:: bash

    PHARE_DRY_RUN=1 python3 my_script.py

In dry run mode ``initialize()``, ``advance()``, and ``run()`` all return immediately
without performing any computation.


Log files
---------

By default (``log_to_file=True``), C++ output from each MPI rank is redirected to a
separate file in the ``.log/`` directory created in the working directory:

- ``RANK_FILES`` mode (default) — one file per rank, e.g., ``.log/rank_0.log``.
- ``DATETIME_FILES`` mode — one file per rank with a timestamp suffix.

The mode is controlled by the ``PHARE_LOG`` environment variable:

.. code-block:: bash

    # Default: per-rank files
    PHARE_LOG=RANK_FILES mpirun -np 4 python3 -Ou my_script.py

    # Timestamped files
    PHARE_LOG=DATETIME_FILES mpirun -np 4 python3 -Ou my_script.py

    # Print everything to stdout (useful for debugging single-rank runs)
    PHARE_LOG=CLI python3 my_script.py

    # Suppress all C++ output
    PHARE_LOG=NULL python3 my_script.py

To disable log files entirely, pass ``log_to_file=False`` to ``Simulator``:

.. code-block:: python

    simulator = Simulator(sim, log_to_file=False)


API reference
-------------

.. autoclass:: pyphare.simulator.simulator.Simulator
   :members:
