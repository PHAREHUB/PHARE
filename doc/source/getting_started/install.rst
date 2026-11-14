=============
Install PHARE
=============

**Goal:** get PHARE's source code and everything it needs, so that you can
:doc:`build <build>` it and run :doc:`your first simulation <first_simulation>`.

PHARE is used in two parts: a compiled C++ core, which does the computation, and
a Python package, ``pyphare``, with which you write input scripts and read the
results. You install the tools for both here, then compile the core on the
:doc:`build` page.

.. mermaid::

   flowchart TD
     A["1. system tools<br/>compiler, CMake, MPI, HDF5, Python"] --> B["2. get the source<br/>git clone"]
     B --> C["3. Python packages<br/>uv sync or pip"]
     C --> D["4. build the C++ core"]
     click D href "build.html"


1. System tools
---------------

Install these with your system's package manager, or load them as modules on a
cluster:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Tool
     - Why PHARE needs it
   * - ``git``
     - to get the source, and the libraries PHARE downloads while building
   * - CMake 3.20.1 or newer
     - prepares the build
   * - ``make`` or ``ninja``
     - runs the compilation
   * - a C++20 compiler (and C and Fortran compilers)
     - compiles PHARE and the SAMRAI library it uses for the mesh
   * - MPI
     - runs a simulation on several processes
   * - parallel HDF5
     - writes the output files (a serial HDF5 builds, but PHARE then runs on
       one process only)
   * - ``m4``
     - used while compiling SAMRAI (usually present on Linux)
   * - Python 3.11 or newer, with its development files (headers)
     - input scripts, analysis, and the link between Python and C++

On a laptop, "development files" means a package often called
``python3-devel`` or ``python3-dev``.


2. Get the source
-----------------

.. code-block:: bash

   git clone https://github.com/PHAREHUB/PHARE

This creates a ``PHARE`` folder. The libraries PHARE depends on (SAMRAI for
the adaptive mesh, HighFive for HDF5 files, cppdict for the input dictionary,
and pybind11 for the link to Python if it is not already installed) are
downloaded automatically into ``PHARE/subprojects`` the first time you
build, so the machine needs internet access at that step.


3. Python packages
------------------

The Python packages PHARE needs are listed in ``pyproject.toml`` and
``requirements.txt`` (``numpy``, ``scipy``, ``h5py``, ``matplotlib`` and a few
others). Choose one of these three ways, from inside the ``PHARE`` folder.

With `uv <https://docs.astral.sh/uv/>`_ (creates a virtual environment in
``.venv``):

.. code-block:: bash

   uv sync
   source .venv/bin/activate

With ``pip``, in a virtual environment:

.. code-block:: bash

   python3 -m venv phare_venv
   source phare_venv/bin/activate
   python3 -m pip install -r requirements.txt

With ``pip``, for your user only:

.. code-block:: bash

   python3 -m pip install --user -r requirements.txt

.. important::

   Use the **same Python** to build PHARE and to run it. The compiled core is
   made for one Python version (for example 3.13) and does not load in
   another one. If you use a virtual environment, activate it before building.


Next
----

:doc:`build` compiles PHARE and checks that Python can find it.
