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

.. plot::
   :include-source: False
   :caption: Setup and initialize first, then one advance per time step. Squares: diagnostics
             written (t = 0 by initialize); post_advance is called after each step that writes.

   from phare_draw import *
   fig, (a,) = strip(1, ["Simulator(sim).run()"], height=2.6, phone=True)
   arrow(a, (0, 0), (10.8, 0), color=C["grey"], lw=1.2)
   text(a, 11.2, 0, "t", color=C["grey"])
   for k in range(11):
       a.plot([k, k], [-0.1, 0.1], color=C["grey"], lw=1.0)
   text(a, 0, -0.5, "0", color=C["grey"])
   text(a, 10, -0.5, "end", color=C["grey"])
   brace(a, 2, 3, 0.25, "1 step", up=True)
   for t in (0, 4, 8):
       mark(a, "dual", t, 1.3, color=C["E"], size=10)
   text(a, 4, 1.9, "post_advance", color=C["ion"])
   arrow(a, (8, 1.9), (8, 1.55), color=C["ion"], lw=1.2)
   text(a, 0.2, 2.6, "setup + initialize", ha="left")
   arrow(a, (0, 2.35), (0, 1.55), color=C["ink"], lw=1.2)
   finish(fig, [a])

.. autoclass:: pyphare.simulator.simulator.Simulator
   :members:
