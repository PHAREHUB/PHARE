
Tutorial: Adaptive Mesh Refinement with Tagging
================================================

This tutorial demonstrates PHARE's adaptive mesh refinement (AMR) tagging
system using a propagating tangential discontinuity as the test case. Where
the previous tutorials fix refined regions with static boxes, tagging lets
the solver decide at runtime which cells need refinement — the grid follows
the physics automatically.

.. contents:: Contents
   :local:
   :depth: 2


Static vs adaptive refinement
------------------------------

PHARE supports two strategies for driving mesh refinement.

**Static boxes** (``refinement="boxes"``) require you to specify the refined
regions explicitly before the simulation starts:

.. code-block:: python

    ph.Simulation(
        refinement_boxes={
            "L0": {"B0": [(50,), (150,)]},
        },
        ...
    )

This is simple and free of runtime overhead, but it works only when you know
in advance where the interesting physics will be. A propagating discontinuity,
a drifting current sheet, or a spontaneously forming shock will quickly leave
the static boxes behind.

**Adaptive tagging** (``refinement="tagging"``) evaluates a criterion — by
default the normalised gradient of the magnetic field magnitude — every
coarse time step, marks cells above a threshold, and creates or removes
refined patches dynamically:

.. code-block:: python

    ph.Simulation(
        refinement="tagging",
        max_nbr_levels=3,
        tagging_threshold=0.5,
        ...
    )

The refined region tracks whatever structure exceeds the gradient threshold,
with no intervention from the user.

.. list-table:: When to choose each strategy
   :header-rows: 1
   :widths: 20 40 40

   * - Strategy
     - Good for
     - Not ideal for
   * - ``"boxes"``
     - Known equilibria, convergence tests, reproducibility studies
     - Moving structures, turbulence, spontaneous instabilities
   * - ``"tagging"``
     - Propagating discontinuities, spontaneous reconnection, shocks
     - Runs where the refined topology must not change (e.g. exact restarts)


Setting up a tagged simulation
-------------------------------

The test case is a one-dimensional tangential discontinuity (TD): a thin
layer where :math:`B_y` reverses sign, carried across the periodic domain
at bulk velocity :math:`v_x = 2`. This reproduces the scenario used in
``tests/functional/tdtagged/td1dtagged.py``.

Create a file ``td_tagged_1d.py``:

.. code-block:: python

    #!/usr/bin/env python3

    import numpy as np
    import pyphare.pharein as ph
    from pyphare.simulator.simulator import Simulator

    ph.NO_GUI()


    diag_dir = "phare_outputs/td_tagged"

    final_time = 20.0
    time_step  = 0.04
    timestamps = np.arange(0, final_time + time_step, final_time / 50)


    # ----------------------------------------------------------------
    # Smooth step function — used to build the By tanh profile.
    # ----------------------------------------------------------------
    def S(x, x0, l):
        return 0.5 * (1.0 + np.tanh((x - x0) / l))


    def config():
        sim = ph.Simulation(
            # --------------------------------------------------------
            # 1-D domain: 200 cells, cell size 1 di.
            # --------------------------------------------------------
            cells=200,
            dl=1.0,
            time_step=time_step,
            final_time=final_time,
            boundary_types="periodic",
            interp_order=1,
            # --------------------------------------------------------
            # Patch sizes: 20 cells keeps patches small enough that the
            # tagging algorithm can concentrate refinement around the
            # two thin discontinuity layers without wasting cells in
            # the quiet background plasma.
            # --------------------------------------------------------
            smallest_patch_size=20,
            largest_patch_size=20,
            # --------------------------------------------------------
            # AMR: adaptive tagging with up to 3 levels.
            # The default tagging_threshold=0.1 is used here; see
            # the "Performance tips" section for guidance on tuning it.
            # --------------------------------------------------------
            refinement="tagging",
            max_nbr_levels=3,
            nesting_buffer=1,
            clustering="tile",
            # --------------------------------------------------------
            # Hyper-resistivity damps grid-scale noise near the
            # discontinuity without broadening the physical layer.
            # --------------------------------------------------------
            hyper_resistivity=0.01,
            diag_options={
                "format": "phareh5",
                "options": {"dir": diag_dir, "mode": "overwrite"},
            },
        )

        # Retrieve the domain length after the Simulation object is built.
        L = sim.simulation_domain()[0]  # 200 di

        # ----------------------------------------------------------------
        # Density: uniform background.
        # ----------------------------------------------------------------
        def density(x):
            return 1.0

        # ----------------------------------------------------------------
        # Magnetic field: two tangential discontinuities in By placed at
        # L/4 and 3L/4.  The S() steps produce +1 between the layers and
        # -1 outside, so By reverses sign twice across the domain.
        # Bz provides an out-of-plane guide field.
        # ----------------------------------------------------------------
        def bx(x):
            return 0.0

        def by(x):
            v1 = -1.0
            v2 =  1.0
            return v1 + (v2 - v1) * (S(x, L * 0.25, 1) - S(x, L * 0.75, 1))

        def bz(x):
            return 0.5

        def b2(x):
            return bx(x) ** 2 + by(x) ** 2 + bz(x) ** 2

        # ----------------------------------------------------------------
        # Temperature: derived from total pressure balance with K = 1.
        # ----------------------------------------------------------------
        def T(x):
            K = 1.0
            return (K - 0.5 * b2(x)) / density(x)

        # ----------------------------------------------------------------
        # Bulk velocity: the TD is carried at vx = 2 in the +x direction.
        # ----------------------------------------------------------------
        def vx(x):
            return 2.0

        def vy(x):
            return 0.0

        def vz(x):
            return 0.0

        # Thermal speed: isotropic Maxwellian with local temperature T(x).
        def vthx(x):
            return T(x)

        def vthy(x):
            return T(x)

        def vthz(x):
            return T(x)

        ph.MaxwellianFluidModel(
            bx=bx,
            by=by,
            bz=bz,
            protons={
                "charge": 1,
                "density": density,
                "vbulkx": vx,
                "vbulky": vy,
                "vbulkz": vz,
                "vthx": vthx,
                "vthy": vthy,
                "vthz": vthz,
            },
        )

        ph.ElectronModel(closure="isothermal", Te=0.12)

        # ----------------------------------------------------------------
        # Standard diagnostics.
        # ----------------------------------------------------------------
        for quantity in ["E", "B"]:
            ph.ElectromagDiagnostics(quantity=quantity, write_timestamps=timestamps)

        for quantity in ["charge_density", "bulkVelocity"]:
            ph.FluidDiagnostics(quantity=quantity, write_timestamps=timestamps)

        # ----------------------------------------------------------------
        # MetaDiagnostics: write the refinement tag field at every output
        # time.  The output file records which cells were tagged above the
        # threshold at each snapshot, making it easy to verify that
        # refinement tracks the discontinuities.
        # ----------------------------------------------------------------
        ph.MetaDiagnostics(quantity="tags", write_timestamps=timestamps)

        return sim


    if __name__ == "__main__":
        Simulator(config()).run()


AMR parameter reference
~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   * - Parameter
     - Default
     - Effect
   * - ``refinement``
     - ``"boxes"``
     - Set to ``"tagging"`` to enable runtime criterion-based refinement.
   * - ``max_nbr_levels``
     - 1
     - Total number of AMR levels including the coarse root level.
       Level 1 means no refinement; level 3 means two nested fine levels.
       Each level halves the cell size (refinement ratio = 2).
   * - ``tagging_threshold``
     - 0.1
     - Normalised gradient of :math:`|\mathbf{B}|` above which a cell is
       tagged for refinement.  Lower values tag more cells and produce
       wider refined regions; higher values restrict refinement to the
       sharpest gradients.
   * - ``nesting_buffer``
     - 0
     - Minimum gap in coarse cells between the edge of any refined patch
       and the edge of the enclosing level.  A value of 1 ensures refined
       patches never touch the coarse-level boundary, which avoids
       interpolation artefacts near the coarse boundary.
   * - ``tag_buffer``
     - 1
     - Number of cells by which the set of tagged cells is grown before
       clustering into patch boxes.  A value of 1 adds a one-cell margin
       around each tagged region, which reduces the risk of the refined
       patch boundary catching up with a moving structure between two
       regridding steps.
   * - ``clustering``
     - ``"tile"``
     - Algorithm that groups tagged cells into rectangular patch boxes.
       See "Performance tips" below.
   * - ``smallest_patch_size``
     - interp-dependent
     - Minimum number of cells in a patch along each direction.  Smaller
       values allow finer-grained load balancing but increase the number
       of patches and the associated MPI communication overhead.
   * - ``largest_patch_size``
     - ``None``
     - Maximum number of cells in a patch.  Setting this prevents a
       single patch from dominating a rank's workload at the cost of
       producing more patches overall.


Monitoring refinement
----------------------

The ``MetaDiagnostics`` block writes the cell tagging field alongside the
usual physics diagnostics.  Load it with ``pharesee`` to inspect which cells
drove refinement at any snapshot:

.. code-block:: python

    from pyphare.pharesee.hierarchy import hierarchy_from

    diag_dir = "phare_outputs/td_tagged"
    tags = hierarchy_from(h5_filename=f"{diag_dir}/meta_tags.h5", times="10.0000000000")

    # Inspect the patch layout at t = 10.
    for ilvl, level in tags.levels().items():
        print(f"Level {ilvl}: {len(level.patches)} patches")
        for patch in level.patches:
            print(f"  origin={patch.origin}  upper={patch.box.upper}")

Each patch object exposes the tag array stored under the ``"tags"`` dataset,
which is 1 where the B-gradient exceeded ``tagging_threshold`` and 0
elsewhere.


Analyzing results
------------------

Load the magnetic field output and compare the tagged run against a uniform
grid reference to see how adaptive refinement improves the resolution of the
discontinuity profile.

.. code-block:: python

    import numpy as np
    import matplotlib.pyplot as plt
    from pyphare.pharesee.run import Run
    from pyphare.pharesee.hierarchy import hierarchy_from

    diag_dir = "phare_outputs/td_tagged"
    run = Run(diag_dir)

    plot_time = 10.0   # physical time at which to compare

    # ----------------------------------------------------------------
    # Merge all patches from all levels onto a single interpolated grid.
    # interp="linear" fills coarse-fine interfaces smoothly.
    # ----------------------------------------------------------------
    B = run.GetB(plot_time, merged=True, interp="linear")
    x_by = B["By"][1][0]
    by   = B["By"][0](x_by)

    # ----------------------------------------------------------------
    # Also load the raw (un-merged) hierarchy to draw patch boundaries.
    # ----------------------------------------------------------------
    BH = run.GetB(plot_time)

    fig, ax = plt.subplots(figsize=(10, 4))

    ax.plot(x_by, by, color="royalblue", label=r"$B_y$ (tagged, merged)")

    # Overlay coloured spans showing the extent of each AMR level.
    level_colors = ["#cccccc", "#aaddff", "#ffddaa"]
    for ilvl, level in BH.levels().items():
        for patch in level.patches:
            x0 = patch.origin[0]
            x1 = (patch.box.upper[0] + 1) * patch.layout.dl[0]
            ax.axvspan(
                x0, x1,
                color=level_colors[ilvl % len(level_colors)],
                ec="k",
                alpha=0.3,
                ymin=ilvl / (len(BH.levels()) + 1),
                ymax=(ilvl + 1) / (len(BH.levels()) + 1),
                label=f"Level {ilvl} patches" if patch is level.patches[0] else None,
            )

    ax.set_xlabel("x")
    ax.set_ylabel(r"$B_y$")
    ax.set_title(f"Tangential discontinuity with adaptive refinement,  t = {plot_time}")
    ax.legend()
    plt.tight_layout()
    plt.show()

At :math:`t = 10` the two discontinuities (initially at :math:`x = 50` and
:math:`x = 150`) have propagated by :math:`v_x \cdot t = 20` cells to
:math:`x \approx 70` and :math:`x \approx 170`. The shaded spans confirm
that the level-2 (and level-3) patches have followed them rather than
remaining at their initial positions.


Tracking the discontinuity position
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The adaptive grid updates every coarse time step, so the refined region
remains centred on the discontinuity throughout the run.  A simple way to
verify this is to track the :math:`B_y = 0` crossing as a function of time:

.. code-block:: python

    import numpy as np
    from pyphare.pharesee.run import Run
    from pyphare.pharesee.hierarchy.fromh5 import get_times_from_h5

    diag_dir = "phare_outputs/td_tagged"
    run = Run(diag_dir)

    times = get_times_from_h5(f"{diag_dir}/EM_B.h5")

    crossings = []
    for t in times:
        B  = run.GetB(t, merged=True, interp="linear")
        x  = B["By"][1][0]
        by = B["By"][0](x)
        # Find the zero crossing nearest to the expected position.
        expected = (50.0 + 2.0 * t) % 200.0
        idx = np.argmin(np.abs(x - expected))
        # Linearly interpolate to sub-cell precision.
        if idx > 0 and idx < len(by) - 1:
            xi   = np.interp(0.0, [by[idx - 1], by[idx + 1]], [x[idx - 1], x[idx + 1]])
        else:
            xi = x[idx]
        crossings.append(xi)

    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    ax.plot(times, crossings, "o-")
    ax.plot(times, (50.0 + 2.0 * np.array(times)) % 200.0, "k--", label="expected")
    ax.set_xlabel("time")
    ax.set_ylabel("discontinuity position")
    ax.legend()
    plt.show()

The measured positions should follow the dashed line closely.  Deviations
larger than a fraction of :math:`\Delta x` at the finest level indicate
either that the tagging threshold is too high (the discontinuity escapes the
refined zone before regridding) or that ``tag_buffer`` is too small.


Performance tips
-----------------

Patch size trade-offs
~~~~~~~~~~~~~~~~~~~~~~

Every patch boundary costs MPI communication (ghost-cell exchanges) and
SAMRAI bookkeeping.  Very small patches (e.g. ``smallest_patch_size=5``)
produce many patches and can saturate the communication layer even though
total cell count is modest.  Very large patches waste cells: a 200-cell
patch covering a 5-cell discontinuity resolves the quiet background
unnecessarily.

As a rule of thumb for 1-D problems:

- Set ``smallest_patch_size`` to roughly 4–5 times the physical width of
  the structure you want to resolve.
- Set ``largest_patch_size`` to a few times the smallest size; a ratio of
  2–4 is typical.
- For 2-D and 3-D, the same guidance applies per direction.

Clustering algorithm
~~~~~~~~~~~~~~~~~~~~~

``clustering="tile"`` (default) partitions the tagged-cell domain into a
regular tiling of rectangular boxes.  It is fast, scales well on many
ranks, and produces uniform patch sizes that simplify load balancing.
Its drawback is that tile patches can include untagged cells, so the refined
level can be somewhat wider than strictly necessary.

``clustering="berger"`` uses the Berger–Rigoutsos algorithm to fit patches
more tightly around the set of tagged cells.  This reduces wasted cells but
is slower and can produce irregular patch geometries that complicate load
balancing.

Prefer ``"tile"`` (the default) unless you have a specific reason to
minimise the number of refined cells, for example in a very memory-constrained
run.

Tagging threshold tuning
~~~~~~~~~~~~~~~~~~~~~~~~~

``tagging_threshold`` is a normalised gradient: the local gradient of
:math:`|\mathbf{B}|` divided by the maximum gradient on the level at that
time step.  Its range is always (0, 1].

- A threshold near 0 tags almost every cell that has any gradient at all,
  effectively covering the whole domain with refined patches — expensive and
  generally not useful.
- A threshold near 1 tags only the single sharpest feature, which can
  cause a thin discontinuity to occasionally escape the refined zone between
  two regridding cycles.
- Values in the range 0.1–0.5 are typical starting points.  Start at 0.5
  and lower gradually while monitoring the output: as long as the discontinuity
  profile is well resolved on the finest level, there is no need to refine
  further.

.. note::

   The ``tagging_threshold`` default is 0.1 when ``refinement="tagging"``
   and ``max_nbr_levels`` is not explicitly set to 1.  This conservative
   value ensures that most gradient features are captured; adjust upward for
   cheaper runs once you have established that the physics is resolved.


See also
--------

- :doc:`../theory/amr` — description of the AMR algorithm, level hierarchy,
  and coarse-fine interpolation strategy used in PHARE
- :doc:`../usage/simulation_inputs` — complete reference for all
  ``ph.Simulation`` parameters
- :doc:`harris_2d` — 2-D reconnection tutorial that also uses
  ``refinement="tagging"`` to track thin current sheets automatically
