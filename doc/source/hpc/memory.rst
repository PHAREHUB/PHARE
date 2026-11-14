======
Memory
======

**Goal:** know how much memory your run uses on each rank, and which input
choices make it grow.

Measure it
----------

Run with the resource monitor (needs ``python3 -m pip install phlop==0.0.45``, see
:doc:`measuring_performance`; without ``phlop`` nothing is recorded and no
error is shown):

.. code-block:: python

   Simulator(sim).run(monitoring=10)   # one sample every 10 seconds

Then read the per-rank files and plot the memory of each rank:

.. code-block:: python

   from pathlib import Path
   import yaml
   import matplotlib.pyplot as plt

   fig, ax = plt.subplots(figsize=(6, 3.5))
   for path in sorted(Path(".phare/stats").glob("rank.*.yaml")):
       stats = yaml.safe_load(path.read_text())
       col = stats["headers"].index("mem_used_mb")
       mem = [float(s["v"].split(",")[col]) for s in stats["snapshots"]]
       ax.plot(mem, label=f"rank {stats['rank']}")
       print(f"rank {stats['rank']}: peak {max(mem):.0f} MB")

   ax.set_xlabel("sample (one every 'interval' seconds)", fontsize=11)
   ax.set_ylabel("memory used (MB)", fontsize=11)
   ax.legend(fontsize=11)
   fig.tight_layout()
   fig.savefig("memory_per_rank.png")

``mem_used_mb`` is the resident memory of the process (the RAM it really
occupies), in MB. Compare the peak with the memory available per core on
your machine: ranks on the same node share the node's memory.

Where memory goes
-----------------

.. mermaid::

   flowchart TD
     R["memory of one rank"] --> P["particles<br/>(usually the largest)"]
     R --> F["fields<br/>(with ghost cells)"]
     R --> S["mesh bookkeeping<br/>(SAMRAI)"]
     P --> P1["own particles +<br/>copies near patch<br/>and level borders"]

In a hybrid run, particles usually dominate. A :term:`rank` stores the particles of its
patches, plus copies of neighbour particles near patch borders
(:term:`patch ghost particles <patch ghost particle>`) and near the borders
of refined levels (:term:`level ghost particles <level ghost particle>`).
Fields are usually smaller, but every :term:`patch` carries extra
:term:`ghost cells <ghost cell>` around it, so many small patches cost more
than a few large ones. The last part, the lists of patches and the
communication plans, is kept by :term:`SAMRAI`, the AMR library PHARE is
built on.

What makes it grow
------------------

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Input
     - Effect on memory
   * - ``nbr_part_per_cell``
     - Particle memory grows in proportion.
   * - Refined area (``max_nbr_levels``, tagging threshold)
     - Each refined cell carries its own particles, on top of the coarse
       level, which keeps its particles underneath: refined regions add
       memory as they grow.
   * - Number of populations
     - Each population has its own particles.
   * - ``smallest_patch_size`` / ``largest_patch_size``
     - Smaller patches: more ghost cells and border copies per useful cell.
   * - Number of ranks
     - More ranks: less memory per rank, but not proportionally (border
       copies and per-rank overheads remain).

Common mistakes
---------------

* **No file in** ``.phare/stats/`` — ``phlop`` is not installed: monitoring
  silently does nothing — ``python3 -m pip install phlop``.
* **Passing** ``monitoring=True`` — samples every second, which makes long
  YAML files — give the interval in seconds, e.g. ``monitoring=60``.
* **Comparing the peak of one rank with the memory of a node** — ranks on
  the same node share its memory — multiply by the number of ranks per
  node.

.. admonition:: In preparation
   :class: in-preparation

   Measured memory per particle and per cell, and memory per rank for
   reference runs.

.. seealso:: How data are laid out in memory: :doc:`/dev/hpc_internals`.
