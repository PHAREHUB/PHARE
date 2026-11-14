=======================
Running a simulation
=======================

**Goal:** run your input script, on one core or on several, and know what the
run prints and writes.

PHARE runs as a Python script: the script describes the simulation, then asks a
``Simulator`` to run it. Python and its packages must be installed first (see
:doc:`/getting_started/install`).

.. contents:: On this page
   :local:
   :depth: 1


Minimal example
---------------

First tell Python where PHARE is: its Python package ``pyphare`` and your build
folder (see :doc:`/getting_started/build`). Assuming the PHARE source is in
``/path/to/PHARE`` and the build folder is ``/path/to/build``:

.. code-block:: bash

   export PYTHONPATH=/path/to/PHARE/pyphare:/path/to/build:$PYTHONPATH

Then write a script that declares a :doc:`simulation`, a model
(:doc:`hybrid_inputs`), and :doc:`diagnostics`, and ends with these two lines
(the complete script is in :doc:`/getting_started/first_simulation`):

.. code-block:: python

   from pyphare.simulator.simulator import Simulator

   # ... ph.Simulation(...), ph.MaxwellianFluidModel(...), ph.ElectronModel(...)

   Simulator(sim).run()

Run it from an empty folder, because PHARE writes its files where you start it:

.. code-block:: bash

   python3 /path/to/my_script.py

If you prefer, put the run in a ``main`` function:

.. code-block:: python

   def main():
       Simulator(sim).run()


   if __name__ == "__main__":
       main()


What happens when you run
-------------------------

``run()`` sets everything up on its own if you have not done it, then advances
in time until ``final_time``. Diagnostics are written at the times you asked
for.

.. mermaid::

   flowchart TD
     A["Your script<br/>sim = ph.Simulation(...)"] --> B["Simulator(sim).run()"]
     B --> C["setup: start MPI, redirect C++ logs,<br/>build the simulation from your inputs"]
     C --> D["initialize: create the particles and fields,<br/>write diagnostics for the initial time"]
     D --> E{"current time<br/>< final_time ?"}
     E -- yes --> F["advance one time step"]
     F --> G["write diagnostics (and restarts) due at this time<br/>then call post_advance, if any"]
     G --> E
     E -- no --> H["print timing summary, end"]

While it runs, rank 0 prints one line per time step (current time, wall-clock
time of the step, total elapsed time) and, at the end, the mean and total
advance time.

You will see a few messages at startup. They are harmless:

* ``Allowing patches smaller than the given smallest patch size to prevent overlaps``
  comes from the grid library (SAMRAI). It is informational.
* ``Warning: PHARE was not built with the configurator active``: an optional
  build feature is off. It does not change the results.


The Simulator object
--------------------

``Simulator(sim, **options)`` takes the ``Simulation`` object, and optionally
the keyword options below. They have good defaults.

.. list-table::
   :header-rows: 1

   * - Option
     - Meaning
     - Default
   * - ``print_one_line``
     - If ``True``, the per-step line overwrites the previous one instead of
       scrolling.
     - ``False``
   * - ``post_advance``
     - A Python function ``f(time)``, called after a step in which diagnostics
       were written. It lets you analyse the data just written before the next
       step. ``time`` is the current simulation time.
     - ``None``
   * - ``log_to_file``
     - If ``True``, text printed by the C++ code is sent to files in the
       ``.log/`` folder, one per MPI rank (see below).
     - ``True``

``run()`` has two options of its own:

.. list-table::
   :header-rows: 1

   * - Option
     - Meaning
     - Default
   * - ``plot_times``
     - If ``True``, saves the wall-clock time of each step as
       ``timestep_times.png`` (rank 0) when the run ends.
     - ``False``
   * - ``monitoring``
     - Resource monitoring (CPU and memory per MPI rank), written to
       ``.phare/stats/rank.<N>.yaml``. Give a number of seconds between
       samples, e.g. ``monitoring=10``. If left unset, the environment variable
       ``PHARE_SIM_MON`` decides. Needs the optional ``phlop`` package
       (``python3 -m pip install phlop``); without it nothing is recorded and
       no error is raised.
     - ``None`` (off)

Example using all of them:

.. code-block:: python

   def after_dump(time):
       print(f"data written at t = {time}")

   simulator = Simulator(sim, print_one_line=True, post_advance=after_dump)
   simulator.run(plot_times=True, monitoring=10)

.. note::

   Pass a number to ``monitoring``. ``monitoring=True`` is accepted, but in
   this version it samples every second (not every 100 s as the code comment
   says).

``Simulator`` also has the methods ``setup()``, ``initialize()`` and
``advance()``. ``run()`` calls them for you; call them yourself only if you
want to step the simulation by hand.


Logs: the ``.log/`` folder and ``PHARE_LOG``
---------------------------------------------

By default (``log_to_file=True``) the C++ output goes to ``.log/<rank>.out``, one
file per MPI rank, in the folder where you started the run. The folder is
created for you. The environment variable ``PHARE_LOG`` changes this:

* ``RANK_FILES`` (default): ``.log/<rank>.out``
* ``DATETIME_FILES``: one file per rank named with the start date and time,
  for example ``.log/2026-11-13-22:54:46_0_of_1.out``
* ``CLI``: no files, C++ output appears in the terminal
* ``NULL``: C++ output is discarded

.. code-block:: bash

   PHARE_LOG=CLI python3 my_script.py

The messages printed by Python (the per-step lines, the start and end times) are
not affected: they always go to the terminal. See :doc:`environment_variables`.


Running with MPI
----------------

To use several cores, start the same script through ``mpirun``:

.. code-block:: bash

   mpirun -n 4 python3 /path/to/my_script.py

Nothing changes in the script. The domain is split among the ranks (processes),
and each rank writes its own log in ``.log/``. Only rank 0 prints the per-step
lines. Diagnostics from all ranks end up in the same output files.

.. note::

   A restart must be run with the same number of ranks as the run that wrote
   it (see :doc:`restarts`).


Checking a script without running it (dry run)
----------------------------------------------

A dry run builds and checks the simulation (models, inputs, grid) but never
creates the particles and never advances: ``initialize()``, ``advance()`` and
``run()`` return immediately. It is a fast way to catch typos in a script.
Turn it on in the script,

.. code-block:: python

   sim = ph.Simulation(
       # ...
       dry_run=True,
   )

or from the shell, without editing the script:

.. code-block:: bash

   PHARE_DRY_RUN=1 python3 my_script.py

A dry run may still create the empty folders ``.log/`` and the diagnostics
folder.


Running PHARE as an executable
------------------------------

PHARE can also be built as a stand-alone program, ``phare-exe``, which embeds
Python. It is not built by default: configure the build with
``-DPHARE_EXE=ON``. You then give it your input script as its only argument,
for example ``phare-exe my_script.py``; the script is read as a Python module
and must therefore be importable from where you run it. Most users should use
the ``python3`` route described above.

.. admonition:: In preparation
   :class: in-preparation

   This page describes the executable from the source code only; it has not
   been tested here.


Common mistakes
---------------

.. list-table::
   :header-rows: 1

   * - Message
     - Cause
     - Fix
   * - ``PHARE_LOG invalid type, valid keys are RANK_FILES/DATETIME_FILES/CLI/NULL``
       (followed by ``Error in Simulator.setup()``)
     - ``PHARE_LOG`` has another value
     - use one of the four keys, in capitals
   * - ``ModuleNotFoundError: No module named 'pyphare'``
     - ``PYTHONPATH`` does not contain the PHARE source
     - see "Minimal example" above
   * - no ``.phare/stats`` folder after ``monitoring=10``
     - the ``phlop`` package is not installed
     - ``python3 -m pip install phlop``
   * - a run that ends at once and writes nothing
     - ``PHARE_DRY_RUN=1`` is set in your shell
     - ``unset PHARE_DRY_RUN``


.. seealso:: How this works inside: :doc:`/dev/one_time_step`


Reference
---------

.. autoclass:: pyphare.simulator.simulator.Simulator
   :members:
