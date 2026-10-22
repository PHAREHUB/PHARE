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


class PolytropicClosure:
    closure_name = "polytropic"

    def __init__(self, **kwargs):
        self.dim = global_vars.sim.ndim
        self.Pe = kwargs.get("Pe", self._defaultPe())
        self.Gamma = kwargs.get("Gamma", PolytropicClosure._defaultGamma())

    def _defaultPe(self):
        if self.dim == 1:
            return lambda x: 0.1
        if self.dim == 2:
            return lambda x, y: 0.1
        if self.dim == 3:
            return lambda x, y, z: 0.1
        return 0.1

    @staticmethod
    def _defaultGamma():
        return 5.0 / 3.0

    def dict_path(self):
        return {
            "name/": PolytropicClosure.closure_name,
            "Pe": self.Pe,
            "Gamma": self.Gamma,
        }

    @staticmethod
    def name():
        return PolytropicClosure.closure_name


class ElectronModel(object):
    """
    ElectronModel sets the closure used to compute the fluid electron
    pressure in a Hybrid simulation (kinetic ions, fluid electrons). This
    pressure enters the generalized Ohm's law used to advance the electric
    field. Required in every Hybrid simulation.

    **Usage example:**

    .. code-block:: python

        import pyphare.pharein as ph

        ph.ElectronModel(closure="isothermal", Te=0.2)

    **Parameters**:

        * **closure** (``str``), currently only "isothermal" is implemented.
        * **Te** (``float``), default=0.1, electron temperature. With the isothermal
          closure, this temperature is constant in both space and time, and the
          electron pressure is simply `Te` times the electron density.
    """

    def __init__(self, **kwargs):
        if kwargs["closure"] == "isothermal":
            self.closure = IsothermalClosure(**kwargs)
        elif kwargs["closure"] == "polytropic":
            self.closure = PolytropicClosure(**kwargs)
        else:
            self.closure = None

        global_vars.sim.set_electrons(self)

    def dict_path(self):
        return [
            ("electrons/pressure_closure/" + k, v)
            for k, v in self.closure.dict_path().items()
        ]
