
==============
Load balancing
==============

On a parallel run, work (particles and cells) is split across MPI ranks by
patch. As the simulation evolves - particles move, AMR levels change shape -
that split can become uneven. ``LoadBalancer`` controls how and when PHARE
corrects for that. It is optional: without one, load balancing still runs
with the defaults below.

.. autoclass:: pyphare.pharein.LoadBalancer

.. code-block:: python

    from pyphare.pharein import LoadBalancer

    LoadBalancer(mode="nppc", tol=0.05, every=200)

Use ``mode="nppc"`` (particles per rank) for a Hybrid simulation - it's the
quantity that actually drives compute cost there. ``mode="homogeneous"``
(cells per rank) is closer to what an MHD-only run needs, since there are no
particles.
