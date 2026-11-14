===================
Grid and Yee layout
===================

**Goal:** know where each quantity lives in a cell, so that you read, plot and compare
PHARE's fields at the right positions.

PHARE uses a :term:`staggered grid`, the :term:`Yee layout`: the components of a vector do not all
sit at the same place in a cell. This is not a detail of output: it is what makes the
discrete curl and divergence behave like the real ones (for example, the divergence of B
stays zero when B is advanced with the curl of E).

Minimal example
---------------

The table of positions is available in Python, without running a simulation:

.. code-block:: python

   from pyphare.core.gridlayout import yee_centering

   for qty in ["Bx", "By", "Bz", "Ex", "Ey", "Ez", "rho", "mhdRho"]:
       print(qty, [yee_centering[d][qty] for d in "xyz"])

It prints, for each quantity, whether it is *primal* or *dual* along x, y and z::

   Bx ['primal', 'dual', 'dual']
   By ['dual', 'primal', 'dual']
   Bz ['dual', 'dual', 'primal']
   Ex ['dual', 'primal', 'primal']
   Ey ['primal', 'dual', 'primal']
   Ez ['primal', 'primal', 'dual']
   rho ['primal', 'primal', 'primal']
   mhdRho ['dual', 'dual', 'dual']

What it means
-------------

Along one direction, with cell size :math:`\Delta x`:

* :term:`primal` values sit on the cell boundaries, at :math:`x_i = i\,\Delta x`;
* :term:`dual` values sit at the cell centres, at :math:`x_i = (i + \tfrac12)\,\Delta x`.

A quantity is primal or dual in each direction separately. In two dimensions this gives
four places in a cell:

.. plot::
   :include-source: False
   :caption: One cell in 2D. Nodes (corners, primal–primal): hybrid ion moments and
      :math:`E_z`, :math:`J_z`. Middle of the vertical sides (primal in x, dual in y):
      :math:`B_x`, :math:`E_y`, :math:`J_y`. Middle of the horizontal sides (dual in x,
      primal in y): :math:`B_y`, :math:`E_x`, :math:`J_x`. Centre (dual–dual): :math:`B_z`
      and, in MHD, the fluid quantities.

   import matplotlib.pyplot as plt

   fig, ax = plt.subplots(figsize=(5, 5))
   ax.plot([0, 1, 1, 0, 0], [0, 0, 1, 1, 0], color="0.3", lw=1.5)
   kw = dict(fontsize=15, ha="center", va="center")
   for x, y in [(0, 0), (1, 0), (0, 1), (1, 1)]:
       ax.plot(x, y, "o", ms=12, color="#2a7ab9", zorder=3)
   ax.text(0.05, -0.15, "n, V, Ez", **kw)
   for x in (0, 1):
       ax.plot(x, 0.5, "s", ms=10, color="#e07a1f", zorder=3)
   ax.text(0.0, 0.66, "Bx, Ey", **kw)
   for y in (0, 1):
       ax.plot(0.5, y, "D", ms=9, color="#3a9d5d", zorder=3)
   ax.text(0.5, 1.1, "By, Ex", **kw)
   ax.plot(0.5, 0.5, "*", ms=16, color="#8e44ad", zorder=3)
   ax.text(0.5, 0.33, "Bz\nMHD: ρ, ρV,\nP, Etot", **kw)
   ax.annotate("", xy=(1, -0.3), xytext=(0, -0.3),
               arrowprops=dict(arrowstyle="<->", color="0.4"))
   ax.text(0.5, -0.38, "Δx", **kw)
   ax.set_xlim(-0.35, 1.35)
   ax.set_ylim(-0.5, 1.25)
   ax.set_aspect("equal")
   ax.axis("off")
   fig.tight_layout()

In 3D the same rule gives the usual Yee cell: B on the face centres, E and J on the edge
midpoints, hybrid moments on the corners, MHD fluid quantities at the cell centre.

So, over a patch of :math:`N` cells in one direction, a primal quantity has :math:`N+1`
values and a dual one has :math:`N`.

Hybrid and MHD differ in one place
----------------------------------

.. list-table::
   :header-rows: 1

   * - Quantity
     - Hybrid
     - MHD
   * - :math:`\mathbf{B}`
     - faces
     - faces
   * - :math:`\mathbf{E}`, :math:`\mathbf{J}`
     - edges
     - edges
   * - density, velocity, pressure
     - **nodes** (ion moments, deposited by the particles)
     - **cell centres** (finite-volume averages)

Ghost cells
-----------

Every :term:`patch` carries a few extra layers of values around its domain, the :term:`ghost cells <ghost cell>`,
so that its stencils can read past its border. Their number is fixed by the numerical
scheme you choose:

.. list-table::
   :header-rows: 1

   * - Model and choice
     - Ghost layers on each side
   * - hybrid, ``interp_order=1``
     - 2
   * - hybrid, ``interp_order=2`` or ``3``
     - 4
   * - MHD, ``reconstruction="Constant"``
     - 2
   * - MHD, any other reconstruction
     - 4

You do not set this number; it matters if you read raw patch arrays with their ghosts.

What you get when you read the data
-----------------------------------

The output files keep each quantity where the solver stores it. When you read them with
``pyphare``, the getters of ``Run`` for B, E, J and the MHD quantities (``GetB``, ``GetE``,
``GetJ``, ``GetMHDV``, ...) by default move every component to the nodes
(``all_primal=True``), so all components share the same positions and array sizes. Pass
``all_primal=False`` to get the values at their own Yee positions, for example to check
:math:`\nabla\cdot\mathbf{B}` exactly. With ``merged=True`` the values are always the raw
ones. See :doc:`/data_analysis/reading_data`.

Common mistakes
---------------

* **Combining raw components as if they were at the same place.** With
  ``all_primal=False``, :math:`B_x` and :math:`B_y` are half a cell apart, and along a
  direction where one is primal and the other dual their arrays differ by one value. Use
  the default getters, or average to a common place yourself.
* **Measuring a divergence or a curl on node values.** Values moved to the nodes are
  averages; the discrete :math:`\nabla\cdot\mathbf{B}` is only zero on the raw face values.
* **Using one coordinate array for all raw quantities.** Each quantity has its own
  coordinates; take them from that quantity's patch data (its ``x``, ``y``, ``z``).

.. seealso:: How this works inside: :doc:`/dev/gridlayout`.
