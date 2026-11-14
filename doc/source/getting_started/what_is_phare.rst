=============
What is PHARE
=============

**Goal:** know in five minutes what PHARE computes, with which models, and
where to go next.

PHARE simulates collisionless plasmas, such as the solar wind or a planetary
magnetosphere, on a mesh that refines itself where the solution needs it.
It is written in C++ for speed and driven from Python: you describe a run in
a Python script, run it, and read the results back in Python.


PHARE in one picture
--------------------

.. plot::

   import numpy as np
   import matplotlib.pyplot as plt
   from matplotlib.patches import Rectangle

   rng = np.random.default_rng(3)
   plt.rcParams.update({"font.size": 12})
   fig, axes = plt.subplots(3, 1, figsize=(6.5, 9.5))

   # 1. hybrid: ion particles + electron fluid on a mesh
   ax = axes[0]
   x = np.linspace(0, 8, 200)
   ax.fill_between(x, 0, 4, color="tab:blue", alpha=0.12)
   for xi in range(9):
       ax.plot([xi, xi], [0, 4], color="0.75", lw=0.8)
   for yi in range(5):
       ax.plot([0, 8], [yi, yi], color="0.75", lw=0.8)
   px, py = rng.uniform(0, 8, 160), rng.uniform(0, 4, 160)
   ax.scatter(px, py, s=10, color="tab:red", zorder=3)
   ax.set_title("Hybrid model: ions are particles (red dots),\n"
                "electrons are a fluid on the mesh (blue)")
   ax.set_xlim(0, 8); ax.set_ylim(0, 4); ax.set_aspect("equal"); ax.axis("off")

   # 2. MHD: one fluid, cell values
   ax = axes[1]
   xx, yy = np.meshgrid(np.arange(8) + 0.5, np.arange(4) + 0.5)
   rho = 1 + 0.5 * np.tanh((xx - 4) / 1.0) * np.cos(yy / 2)
   ax.imshow(rho, extent=(0, 8, 0, 4), origin="lower", cmap="Blues",
             vmin=0, vmax=2)
   for xi in range(9):
       ax.plot([xi, xi], [0, 4], color="0.6", lw=0.8)
   for yi in range(5):
       ax.plot([0, 8], [yi, yi], color="0.6", lw=0.8)
   ax.set_title("MHD model: the whole plasma is one fluid,\n"
                "described by values on the mesh")
   ax.set_xlim(0, 8); ax.set_ylim(0, 4); ax.set_aspect("equal"); ax.axis("off")

   # 3. AMR: a finer patch where the solution needs it
   ax = axes[2]
   for xi in range(9):
       ax.plot([xi, xi], [0, 4], color="0.6", lw=0.8)
   for yi in range(5):
       ax.plot([0, 8], [yi, yi], color="0.6", lw=0.8)
   x0, y0, w, h = 3, 1, 3, 2
   ax.add_patch(Rectangle((x0, y0), w, h, fill=True, color="tab:orange",
                          alpha=0.15, lw=0))
   for k in range(2 * w + 1):
       ax.plot([x0 + k / 2, x0 + k / 2], [y0, y0 + h], color="tab:orange", lw=0.8)
   for k in range(2 * h + 1):
       ax.plot([x0, x0 + w], [y0 + k / 2, y0 + k / 2], color="tab:orange", lw=0.8)
   ax.add_patch(Rectangle((x0, y0), w, h, fill=False, color="tab:orange", lw=2.5))
   ax.text(4.5, 3.35, "level 1: cells 2 times smaller", ha="center",
           color="tab:orange")
   ax.set_title("Mesh refinement: level 0 covers the domain,\n"
                "finer levels only where needed")
   ax.set_xlim(0, 8); ax.set_ylim(0, 4); ax.set_aspect("equal"); ax.axis("off")

   fig.tight_layout()


The two models
--------------

**Hybrid** (the default). Ions are followed as many
:term:`macro-particles <macro-particle>`, so their full velocity
distribution is kept: beams, temperature anisotropy, non-Maxwellian tails.
Electrons are a massless fluid; the electric field comes from their momentum
equation (:term:`Ohm's law`), assuming :term:`quasineutrality`. This is the
model to use when ion kinetic physics matters, at scales from a fraction of
the ion inertial length upward. See :doc:`/theory/hybridpic`.

**MHD**. The plasma is a single conducting fluid, solved with a
:term:`finite volume` method. It is much cheaper and is the model to use for
large scales, where ion kinetic effects do not matter. See
:doc:`/theory/mhd`.

You choose the model in the input script, and a run uses one model on all
its levels: the two cannot be combined in one run yet (for example MHD on
coarse levels and hybrid on fine levels). See
:doc:`/user_guide/choosing_a_model`.


Mesh refinement
---------------

The domain is first covered by a coarse mesh, *level 0*. Where the solution
needs more resolution, PHARE adds finer :term:`levels <level>`, each 2 times
finer than the one below and made of rectangular :term:`patches <patch>`.
Finer regions can be fixed by you, or follow the solution as it evolves
(:term:`tagging`). This is :term:`adaptive mesh refinement <AMR>`, done with
the :term:`SAMRAI` library. See :doc:`/user_guide/amr`.

PHARE runs in parallel on many cores with :term:`MPI`. See
:doc:`/hpc/running_at_scale`.


What a run looks like
---------------------

#. **Describe** the run in a Python script: domain and mesh, time step,
   model and initial conditions, what to write to disk.
#. **Run** it, on a laptop or on a cluster.
#. **Read** the output files back with PHARE's Python tools and plot them.

:doc:`first_simulation` does all three on a small example.


Where to go next
----------------

* :doc:`install`, then :doc:`build`: get PHARE on your machine.
* :doc:`first_simulation`: your first run, step by step.
* :doc:`/user_guide/simulation`: every option of a run.
* :doc:`/dev/index`: how PHARE works inside, for developers.
