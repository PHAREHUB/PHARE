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

.. plot::
   :include-source: False
   :caption: On a refined level, each rank also stores level ghost particles: fine particles split from coarse ones, in a thin layer just outside the refined region.

   from phare_draw import *
   import random
   random.seed(3)
   fig, (b,) = strip(1, ["Refined region over level 0: level ghost particles just outside it"], height=3.0)
   grid(b, (0, 0, 10, 5))
   level(b, (3, 1, 7, 4))
   coarse = [(random.uniform(0.05, 9.95), random.uniform(0.05, 4.95)) for _ in range(45)]
   fine = [(random.uniform(3.05, 6.95), random.uniform(1.05, 3.95)) for _ in range(45)]
   ring = [p for p in [(random.uniform(2.5, 7.5), random.uniform(0.5, 4.5)) for _ in range(200)]
           if not (2.75 < p[0] < 7.25 and 0.75 < p[1] < 4.25)][:26]
   particles(b, coarse, color=C["grey"], size=18)
   particles(b, fine, color=C["ion"], size=26)
   particles(b, ring, color=C["ghost_edge"], size=30, marker="s")
   particles(b, [(0.0, -0.7)], color=C["grey"], size=18); text(b, 0.25, -0.7, "level 0", color=C["grey"], ha="left")
   b.add_patch(__import__("matplotlib").patches.Rectangle((2.05, -0.95), 0.5, 0.5, fc="none", ec=C["L1"], lw=2.5, zorder=3))
   particles(b, [(2.3, -0.7)], color=C["ion"], size=26); text(b, 2.75, -0.7, "refined level", color=C["ion"], ha="left")
   particles(b, [(5.6, -0.7)], color=C["ghost_edge"], size=30, marker="s"); text(b, 5.85, -0.7, "level ghost particles", color=C["ghost_edge"], ha="left")
   finish(fig, [b])

In a hybrid run, particles usually take most of the memory: each rank stores the particles of its
patches and, on refined levels, the :term:`level ghost particles <level ghost particle>` drawn above.
Particles that leave a patch during a step (:term:`patch ghost particles <patch ghost particle>`) are
held only until they are handed to the neighbouring patch. Fields are usually smaller, but every
:term:`patch` carries :term:`ghost cells <ghost cell>`, so many small patches cost more than a few
large ones. :term:`SAMRAI` also keeps the patch lists and communication plans.

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
     - Smaller patches: more ghost cells per useful cell.
   * - Number of ranks
     - More ranks: less memory per rank, but not proportionally (ghost
       cells and per-rank overheads remain).

Common mistakes
---------------

* **No file in** ``.phare/stats/`` — ``phlop`` is not installed: monitoring
  silently does nothing — ``python3 -m pip install phlop==0.0.45``.
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
