=============
MHD numerics
=============

How the MHD solver advances one step: the finite-volume scheme and how to choose its
options, then each part in detail. The grid, the time step, boundaries, mesh refinement and
the divergence of B are shared with the hybrid solver: see :doc:`common`.

.. toctree::
   :maxdepth: 1

   mhd_finite_volume
   mhd_reconstruction
   mhd_riemann_solvers
   mhd_time_integration
