================
Hybrid numerics
================

How the hybrid solver (kinetic ions as macro-particles, fluid electrons) advances one step:
the time scheme, Ohm's law on the grid (written so that it conserves energy), and the particles
(shapes, pusher, interpolation, deposit). The grid,
the time step and mesh refinement are shared with MHD (domain boundaries other than
periodic exist for MHD only): see :doc:`common`.

.. toctree::
   :maxdepth: 1

   hybrid_time_integration
   hybrid_ohm
   particles
