from . import global_vars


class IsothermalClosure:
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
        if "Pe" not in kwargs:
            raise ValueError(
                "ElectronModel: the polytropic closure needs the initial electron pressure 'Pe'"
            )
        for name in ("gamma", "Gamma"):
            if name in kwargs:
                raise ValueError(
                    f"ElectronModel: '{name}' is not an option of the closure, "
                    "the polytropic index is the 'gamma' given to Simulation"
                )
        self.Pe = kwargs["Pe"]
        self.gamma = global_vars.sim.gamma

    def dict_path(self):
        return {
            "name/": PolytropicClosure.closure_name,
            "Pe": self.Pe,
            "gamma": self.gamma,
        }

    @staticmethod
    def name():
        return PolytropicClosure.closure_name


class ElectronModel:
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

        * **closure** (``str``), "isothermal" or "polytropic".
        * **Te** (``float``), default=0.1, electron temperature, isothermal closure
          only. This temperature is constant in both space and time, and the
          electron pressure is simply `Te` times the electron density.
        * **Pe** (``callable``), required with the polytropic closure: the initial
          electron pressure, a function of the coordinates (same arguments as the
          density).

        The polytropic index of the polytropic closure is the ``gamma`` given to
        ``Simulation`` (default 5/3), which is also the adiabatic index of the MHD fluid.
    """

    def __init__(self, **kwargs):
        if kwargs["closure"] == "isothermal":
            self.closure = IsothermalClosure(**kwargs)
        elif kwargs["closure"] == "polytropic":
            self.closure = PolytropicClosure(**kwargs)
        else:
            raise ValueError(
                f"ElectronModel: unknown closure {kwargs['closure']!r}, "
                "expected 'isothermal' or 'polytropic'"
            )

        global_vars.sim.set_electrons(self)

    def dict_path(self):
        return [
            ("electrons/pressure_closure/" + k, v)
            for k, v in self.closure.dict_path().items()
        ]
