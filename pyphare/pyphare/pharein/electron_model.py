from . import global_vars


class IsothermalClosure(object):
    closure_name = "isothermal"

    def __init__(self, **kwargs):
        self.Te = kwargs.get("Te", IsothermalClosure._defaultTe())

    @staticmethod
    def _defaultTe():
        return 0.1

    def dict_path(self):
        return {"name/": IsothermalClosure.closure_name, "Te": self.Te}

    @staticmethod
    def name():
        return IsothermalClosure.closure_name


class ElectronModel(object):
    """
    ElectronModel sets the closure used to compute the fluid electron
    pressure in a Hybrid simulation (kinetic ions, fluid electrons). This
    pressure enters the generalized Ohm's law used to advance the electric
    field. Required in every Hybrid simulation.

    **Usage example:**

    .. code-block:: python

        from pyphare.pharein import ElectronModel

        ElectronModel(closure="isothermal", Te=0.2)

    **Parameters**:

        * **closure** (``str``), currently only "isothermal" is implemented.
        * **Te** (``float``), default=0.1, electron temperature. With the isothermal
          closure, this temperature is constant in both space and time, and the
          electron pressure is simply `Te` times the electron density.
    """

    def __init__(self, **kwargs):
        if kwargs["closure"] == "isothermal":
            self.closure = IsothermalClosure(**kwargs)
        else:
            self.closure = None

        global_vars.sim.set_electrons(self)

    def dict_path(self):
        return [
            ("electrons/pressure_closure/" + k, v)
            for k, v in self.closure.dict_path().items()
        ]
