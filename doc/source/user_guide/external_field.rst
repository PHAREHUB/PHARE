=======================
External magnetic field
=======================

Either model can be given an external (background) magnetic field B0. It is
declared with one of the blocks below, at most once, right after the
``Simulation`` and before the model; declaring it after the model, or a second
time, is an error. Without such a block there is no B0: it is neither
allocated nor dumped.

B0 is computed on every level (and recomputed after regrids and restarts),
dumped with ``ElectromagDiagnostics(quantity="B0", ...)`` and written to
restart files. No solver uses it yet: declaring it does not change how the
simulation evolves, and it does not add to the ``bx``, ``by``, ``bz`` given to
the model.

.. autoclass:: pyphare.pharein.DipoleExternalField

.. autoclass:: pyphare.pharein.UserDefinedExternalField

.. code-block:: python

    import pyphare.pharein as ph

    sim = ph.Simulation(...)

    def az(x, y, t):
        return x * y * (1.0 + t)

    def daz_dt(x, y, t):
        return x * y

    ph.UserDefinedExternalField(a0z=az, da0z_dt=daz_dt)

    ph.MaxwellianFluidModel(...)
    ph.ElectromagDiagnostics(quantity="B0", write_timestamps=timestamps)
