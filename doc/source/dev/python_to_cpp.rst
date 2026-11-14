.. _dev-python-to-cpp:

==========================
From Python script to C++
==========================

.. container:: breadcrumb

   You are here: :doc:`Developer guide </dev/index>` › From Python script to C++

**In one sentence:** the Python input script fills one C++ dictionary of
parameters, the :term:`PHAREDict`; then Python imports the compiled module
that matches the run's main options, and that module builds the C++
simulator from the dictionary.

The picture
-----------

.. flow::

   Input script | ph.Simulation(...), models, diagnostics
   Simulator.setup()
   1 · populateDict() | pyphare.pharein
   typed add_* calls
   PHAREDict | C++ singleton, filled by Python
   2 · cpp_lib(sim) | imports pybindlibs.cpp_<options>
   3 · make_hierarchy() | cpp_etc: SAMRAI hierarchy, reads PHAREDict | /dev/resources_manager
   4 · make_simulator(hier) | Simulator<options>, reads PHAREDict
   then
   initialize(), advance() in a loop | | /dev/one_time_step

Steps 1 to 4 happen in that order inside ``Simulator.setup()``
(:src:`pyphare/pyphare/simulator/simulator.py#L105-L121`).


How it works
------------

#. **The script declares the run.** ``ph.Simulation(...)``, the model, the
   electrons and the diagnostics are plain Python objects. Nothing reaches
   C++ yet.
#. **The parameters go to C++.** ``populateDict()``
   (:src:`pyphare/pyphare/pharein/__init__.py#L75-L91`) walks these objects
   and writes every value into the ``PHAREDict`` under a path such as
   ``simulation/AMR/tag_buffer``. It writes through the small pybind module
   ``dictator`` (:src:`src/initializer/dictator.cpp#L41-L61`), which has one
   function per C++ type (``add_int``, ``add_double``, ``add_string``, ...).
   Initial profiles such as the density are Python functions: they are
   stored as callable objects (``addInitFunction1D/2D/3D``) and called later
   from C++.
#. **Python picks the compiled module.** The main options (dimension,
   interpolation order, refined particle number, and the MHD scheme for MHD
   runs) give a module name such as ``pybindlibs.cpp_1_1_2``
   (``simulator_id`` and ``cpp_lib``,
   :src:`pyphare/pyphare/cpp/__init__.py#L8-L81`). Python imports it.
#. **C++ builds the mesh hierarchy.** ``cpp_etc.make_hierarchy()`` calls
   ``PHARE::make_hierarchy`` (:src:`src/simulator/simulator_runtime.hpp#L51-L63`),
   which reads the dimension and the grid from the dictionary and creates the
   SAMRAI hierarchy.
#. **C++ builds the simulator.** The module's ``make_simulator`` creates
   ``Simulator<options>`` and calls its ``setup`` with the dictionary
   (``makeSimulator``, :src:`src/simulator/simulator.hpp#L629-L635`).
   Every C++ object then reads its own part of the dictionary in its
   constructor (solvers, messengers, diagnostics, ...).


The tricks
----------

One compiled module per combination of options
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** Dimension, interpolation order, refined particle number and the
MHD scheme are template parameters: they fix the types of almost every C++
object, so the compiler can produce fast loops. They cannot be chosen at run
time.

**Fix.** Each line of :src:`res/sim/all.txt` is one :term:`build
permutation`. At build time, :src:`src/python3/CMakeLists.txt#L9-L58` turns
each line into its own module, compiled from the same source file with the
line as a macro (``PHARE_SIM_STR``). In that module,
``resolve_simulator_options()`` turns the macro into the ``SimOpts`` of the
module (:src:`src/python3/cpp_simulator.hpp#L34-L38`). A line has 3 fields
for a hybrid build (``dim,interp,refined_particle_nbr``), 5 for an MHD build
(``dim,reconstruction,limiter,riemann,hall``), 7 when it carries both.

**What you see when it is missing.** If the run needs a combination that was
not built, ``cpp_lib`` raises ``ModuleNotFoundError`` with the exact line to
add to ``res/sim/all.txt`` and the list of modules that exist
(:src:`pyphare/pyphare/cpp/__init__.py#L55-L69`). Add the line, rebuild.

The dictionary is a singleton
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem.** Python fills the parameters before any C++ simulator exists,
and many C++ objects, created at different times, need to read them.

**Fix.** ``PHAREDictHandler::INSTANCE()`` holds one dictionary for the whole
process (:src:`src/initializer/data_provider.hpp#L77-L96`). Python writes
into it; C++ constructors read from it. ``stop()`` deletes it. The Python
``Simulator`` calls it (``clearDict``,
:src:`pyphare/pyphare/pharein/__init__.py#L66-L72`) when it is reset, so a
second simulation in the same Python process does not see stale values.

The stand-alone executable takes another road
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``phare`` executable (:src:`src/phare/phare.cpp`) runs an input script
without the Python ``Simulator``. It starts an embedded Python interpreter
that imports the script and fills the dictionary
(``PythonDataProvider::read``,
:src:`src/initializer/python_data_provider.hpp#L31-L35`), then builds the
hierarchy and finds the simulator with ``makeSimulatorAtRuntime``
(:src:`src/simulator/simulator_runtime.hpp#L19-L32`). That function loops
over a fixed, short list of combinations compiled into the executable,
``phare_exe_default_simulators()``
(:src:`src/simulator/simulator_def.hpp#L169-L184`): 1D, 2D and 3D, orders 1
to 3, the default refined particle number, hybrid only. Use the Python
``Simulator`` for any other combination: for a run outside this list,
``makeSimulatorAtRuntime`` finds no simulator and returns an empty pointer.


Go deeper
---------

* :doc:`/dev/one_time_step`: what ``advance()`` does once the simulator exists.
* :doc:`/dev/resources_manager`: how the hierarchy holds PHARE data.

Related
-------

* :doc:`/dev/big_map`: where these layers sit in the code.
* :doc:`/user_guide/simulation`: the ``Simulation`` block, from the user side.
* :doc:`/dev/working_on_phare`: building PHARE and choosing permutations.
