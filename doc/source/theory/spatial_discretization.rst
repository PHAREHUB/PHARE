
======================
Spatial discretization
======================

PHARE discretizes fields and current density on a structured Cartesian grid.
The physical equations (see :doc:`hybridpic`) involve curl and gradient operators
that must be approximated numerically on that grid. The choice of where each
field component lives on the grid — its *centering* — determines the accuracy
and conservation properties of those operators. PHARE follows the Yee lattice
convention.


The Yee lattice
---------------

On a uniform grid with cell size :math:`\Delta x` (and :math:`\Delta y`,
:math:`\Delta z` in higher dimensions), not all field components are stored at
the same location within a cell. Instead, each component is offset to a
specific position, either at a cell *node* (primal centering) or at a cell
*face/edge midpoint* (dual centering).

In 1D, a primal point :math:`i` sits at position :math:`x_i = i\,\Delta x`,
while the corresponding dual point sits halfway between two primal nodes at
:math:`x_{i+1/2} = (i + \tfrac{1}{2})\,\Delta x`. In 2D and 3D each direction
is independently primal or dual.

This staggering is not arbitrary. It ensures that the curl finite-difference
operators (Faraday and Ampere) are second-order accurate without requiring any
averaging, and that :math:`\nabla\cdot\mathbf{B} = 0` is maintained to machine
precision by construction.


Field centering table
---------------------

The centering of each quantity in PHARE follows the standard Yee convention.
"P" denotes primal and "D" denotes dual centering along the corresponding axis.

**1D**

+----------+-----+
| Quantity | x   |
+==========+=====+
| Bx       | P   |
+----------+-----+
| By       | D   |
+----------+-----+
| Bz       | D   |
+----------+-----+
| Ex       | D   |
+----------+-----+
| Ey       | P   |
+----------+-----+
| Ez       | P   |
+----------+-----+
| Jx       | D   |
+----------+-----+
| Jy       | P   |
+----------+-----+
| Jz       | P   |
+----------+-----+
| rho      | P   |
+----------+-----+

**2D**

+----------+--------+--------+
| Quantity | x      | y      |
+==========+========+========+
| Bx       | P      | D      |
+----------+--------+--------+
| By       | D      | P      |
+----------+--------+--------+
| Bz       | D      | D      |
+----------+--------+--------+
| Ex       | D      | P      |
+----------+--------+--------+
| Ey       | P      | D      |
+----------+--------+--------+
| Ez       | P      | P      |
+----------+--------+--------+
| Jx       | D      | P      |
+----------+--------+--------+
| Jy       | P      | D      |
+----------+--------+--------+
| Jz       | P      | P      |
+----------+--------+--------+
| rho      | P      | P      |
+----------+--------+--------+

**3D**

+----------+--------+--------+--------+
| Quantity | x      | y      | z      |
+==========+========+========+========+
| Bx       | P      | D      | D      |
+----------+--------+--------+--------+
| By       | D      | P      | D      |
+----------+--------+--------+--------+
| Bz       | D      | D      | P      |
+----------+--------+--------+--------+
| Ex       | D      | P      | P      |
+----------+--------+--------+--------+
| Ey       | P      | D      | P      |
+----------+--------+--------+--------+
| Ez       | P      | P      | D      |
+----------+--------+--------+--------+
| Jx       | D      | P      | P      |
+----------+--------+--------+--------+
| Jy       | P      | D      | P      |
+----------+--------+--------+--------+
| Jz       | P      | P      | D      |
+----------+--------+--------+--------+
| rho      | P      | P      | P      |
+----------+--------+--------+--------+

The pattern is: each component of :math:`\mathbf{B}` is primal along its own
direction and dual along the others. Each component of :math:`\mathbf{E}` and
:math:`\mathbf{J}` is dual along its own direction and primal along the others.
The scalar density :math:`\rho` is primal in every direction.


Finite difference operators
---------------------------

Faraday's law
~~~~~~~~~~~~~

The Maxwell-Faraday equation advances the magnetic field:

.. math::

   \frac{\partial \mathbf{B}}{\partial t} = -\nabla \times \mathbf{E}

In 1D the x-component of the curl of any vector field vanishes (no y or z
derivatives exist), so only :math:`B_y` and :math:`B_z` are advanced:

.. math::

   \frac{\partial B_y}{\partial t} = -\frac{\partial E_z}{\partial x}, \qquad
   \frac{\partial B_z}{\partial t} = +\frac{\partial E_y}{\partial x}

Both right-hand sides involve the derivative of a primal field (:math:`E_z` and
:math:`E_y` are primal in x) evaluated at a dual location (where :math:`B_y`
and :math:`B_z` live). The centered finite difference is:

.. math::

   \left.\frac{\partial E_z}{\partial x}\right|_{i+1/2}
   = \frac{E_z^{i+1} - E_z^{i}}{\Delta x}

where the superscript denotes the primal index. This is a simple two-point
stencil, second-order accurate, requiring no averaging because the centering
of the output matches the centering demanded by the update equation.

Ampere's law
~~~~~~~~~~~~

The static Maxwell-Ampere equation provides the current density from the
magnetic field (displacement current is neglected in the hybrid model):

.. math::

   \mu_0 \mathbf{j} = \nabla \times \mathbf{B}

In 1D this gives only :math:`J_y` and :math:`J_z` (the x-component of the
curl vanishes):

.. math::

   \mu_0 J_y = +\frac{\partial B_z}{\partial x}, \qquad
   \mu_0 J_z = -\frac{\partial B_y}{\partial x}

Here the derivative of a dual field (:math:`B_z`, :math:`B_y` are dual in x)
is evaluated at a primal location (where :math:`J_y`, :math:`J_z` live):

.. math::

   \left.\frac{\partial B_z}{\partial x}\right|_{i}
   = \frac{B_z^{i+1/2} - B_z^{i-1/2}}{\Delta x}

Again, a two-point stencil, second-order accurate.


Gradient operator
-----------------

The generalized Ohm's law (see :doc:`hybridpic`) contains the electron
pressure gradient term :math:`\nabla P_e / (en)`. The electron pressure
:math:`P_e = n T_e` is primal in all directions (it lives at the same nodes
as :math:`\rho`). The gradient component along x is required at a dual
location (where :math:`E_x` lives):

.. math::

   \left.\frac{\partial P_e}{\partial x}\right|_{i+1/2}
   = \frac{P_e^{i+1} - P_e^{i}}{\Delta x}

This is the same primal-to-dual derivative stencil used in the Faraday
operator.


Ghost cells
-----------

When a computational domain is split into patches for adaptive mesh refinement,
each patch does not hold information about its neighbors during a stencil
evaluation. Ghost cells (also called halo or guard cells) are extra layers of
cells appended at each side of a patch boundary. They are filled from
neighboring patches or from coarser-level interpolation before any differential
operator is applied, so that stencils at the physical boundary of a patch can
use the same formulas as interior points.

The number of ghost layers required on each side equals the half-width of the
widest stencil used by the particle interpolation, which in turn equals the
interpolation order set by ``interp_order`` (see
:doc:`../usage/simulation_inputs`):

- order 1 → 1 ghost cell per side
- order 2 → 2 ghost cells per side
- order 3 → 3 ghost cells per side

A schematic of a 1D patch with two ghost layers on each side::

    |  ghost  |  ghost  | cell_0 | cell_1 | ... | cell_N | cell_N+1 |  ghost  |  ghost  |
    ^         ^         ^                                  ^          ^         ^
  g_{-2}   g_{-1}    x_0                               x_{N+1}    g_{N+2}   g_{N+3}

Physical mesh cells run from ``cell_0`` to ``cell_N+1``; the ghost layers on
each side are filled by the messenger layer of PHARE (see
:doc:`../theory/amr`) before each solve step.
