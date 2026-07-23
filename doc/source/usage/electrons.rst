
==============
Electron model
==============

In a Hybrid simulation, ions are kinetic but electrons are a fluid: their
momentum equation (via a chosen pressure closure) is used to compute the
electric field, assuming quasineutrality. ``ElectronModel`` sets that
closure, and is required in every Hybrid simulation.

.. autoclass:: pyphare.pharein.ElectronModel

.. code-block:: python

    from pyphare.pharein import ElectronModel

    ElectronModel(closure="isothermal", Te=0.2)

Not needed for an MHD simulation (see :doc:`models`).
