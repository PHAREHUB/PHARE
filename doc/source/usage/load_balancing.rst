
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

Two different scopes
---------------------

Not all parameters act on the same part of the hierarchy:

* ``mode`` shapes how load is measured everywhere - at any level, whenever
  that level's patches are redistributed across ranks.
* ``tol`` does too, but is also reused as the threshold for a second,
  separate mechanism: proactively rebalancing the coarsest level (L0)
  between regrids.
* ``auto``, ``every`` and ``on_init`` control only that second mechanism.
  They have no effect on finer levels - those are only ever redistributed
  across ranks as part of their own AMR regrid, on their own cadence.

In practice, this makes ``auto``/``every``/``on_init`` mostly relevant to
Hybrid simulations: L0 covers the whole domain and typically holds most of
the particles, so it is the level whose imbalance matters most between
regrids.
