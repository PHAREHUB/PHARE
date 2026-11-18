=======================
External magnetic field
=======================

Either model can be given an external (background) magnetic field B0. It is
declared with one of the blocks below, right after the
``Simulation`` and before the model; declaring it after the model, or a second
time, is an error. Without such a block there is no B0: it is neither
allocated nor dumped.

No solver uses it yet.

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
