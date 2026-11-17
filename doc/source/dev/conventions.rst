.. _dev-conventions:

===========
Conventions
===========

.. container:: breadcrumb

   You are here: :doc:`/dev/index` › :doc:`/dev/working_on_phare` › Conventions.

**In one sentence:** the rules the code base follows, so that a change looks like the code around it
and does not break at run time on another machine.

.. admonition:: In preparation
   :class: in-preparation

   C++ conventions (layout of ``src/``, headers, classes, naming, comments) are being written in
   :pr:`1463` (C++ coding conventions), following :issue:`1317` (no style conventions). This page will
   be updated when they are agreed. Until then, ``doc/conventions.md`` in the repository is the
   reference file.

Formatting and lint
-------------------

* C++: ``.clang-format`` at the repository root (version pinned in ``pyproject.toml``).
* Python: lint rules in ``ruff.toml``. While ``ruff.toml`` exists, ruff ignores the ``[tool.ruff]``
  table of ``pyproject.toml``.

Python dependencies
-------------------

Third-party Python packages are declared in ``pyproject.toml``: what a run needs in
``[project].dependencies``, what tests, docs and tools need in ``[dependency-groups]``. The CI image
and Read the Docs do not install from ``pyproject.toml`` yet, so ``requirements.txt`` (run) and
``doc/source/requirements.txt`` (docs) repeat these lists: a new package goes in both places.
Fewer dependencies is better; weigh the cost of each new one.

Python imports: some packages only inside functions
---------------------------------------------------

**Why.** ``h5py`` and ``mpi4py`` load the system HDF5 and MPI libraries as soon as they are
imported. These may not be the libraries PHARE was built with, and the mismatch breaks the run.
``scipy.optimize`` needs system libraries that may be missing. A file that imports them at the top
fails as soon as it is imported, even by the test scanner, before any of its functions is used.

**Rule.** Never import ``h5py``, ``mpi4py`` or ``scipy.optimize`` at the top of a file; import them
in the function that uses them. ``numpy`` and similar widely used packages are fine at the top.

.. code-block:: python

   def read_times(path):
       import h5py  # here, not at the top of the file

       with h5py.File(path, "r") as f:
           return list(f["t"].keys())
