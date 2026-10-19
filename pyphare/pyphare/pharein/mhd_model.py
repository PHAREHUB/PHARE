import numbers
from . import global_vars


class MHDModel(object):
    """
    MHDModel sets up the MHD simulation type: the initial density, velocity,
    magnetic field and pressure profiles of the single MHD fluid.

    Used for an MHD simulation, i.e. one declared with
    ``Simulation(model_options=["MHDModel"], ...)``. For a Hybrid simulation
    (the default), use :class:`~pyphare.pharein.MaxwellianFluidModel` instead
    - PHARE runs either Hybrid or MHD, never both together in the same run.

    **Usage example:**

    .. code-block:: python

        ph.MHDModel(
            density=density,
            vx=vx, vy=vy, vz=vz,
            bx=bx, by=by, bz=bz,
            p=pressure,
        )

    **Parameters**:

        * **density** (*function or number*), default=1, mass density profile.
        * **vx**, **vy**, **vz** (*function or number*), default=(1, 0, 0), velocity profile, per direction.
        * **bx**, **by**, **bz** (*function or number*), default=(1, 0, 0), total magnetic field profile, per direction.
        * **b1x**, **b1y**, **b1z** (*function or number*), default=(0, 0, 0), magnetic field profile minus the external field
          of ``Simulation(external_field=...)``, per direction. Exclusive with **bx**, **by**, **bz**.
        * **p** (*function or number*), default=1, thermal pressure profile.

    Each parameter is a function of the spatial coordinates (one argument per
    dimension: `f(x)`, `f(x, y)` or `f(x, y, z)`), or a plain number for a
    spatially uniform value.
    """

    def defaulter(self, input, value):
        if input is not None:
            if isinstance(input, numbers.Real):
                value = input
            elif not callable(input):
                raise TypeError(
                    f"model profile must be a number or a function, got {type(input).__name__}"
                )
            else:
                import inspect

                params = list(inspect.signature(input).parameters.values())
                assert len(params)
                param_per_dim = len(params) == self.dim
                has_vargs = params[0].kind == inspect.Parameter.VAR_POSITIONAL
                assert param_per_dim or has_vargs
                return input
        if self.dim == 1:
            return lambda x: value + x * 0
        if self.dim == 2:
            return lambda x, y: value
        if self.dim == 3:
            return lambda x, y, z: value

    def __init__(
        self,
        density=None,
        vx=None,
        vy=None,
        vz=None,
        bx=None,
        by=None,
        bz=None,
        p=None,
        b1x=None,
        b1y=None,
        b1z=None,
    ):
        if global_vars.sim is None:
            raise RuntimeError("A simulation must be declared before a model")

        if global_vars.sim.model is not None:
            raise RuntimeError("A model is already created")

        self.dim = global_vars.sim.ndim

        density = self.defaulter(density, 1.0)
        vx = self.defaulter(vx, 1.0)
        vy = self.defaulter(vy, 0.0)
        vz = self.defaulter(vz, 0.0)
        p = self.defaulter(p, 1.0)

        b1_given = any(b is not None for b in (b1x, b1y, b1z))
        if b1_given and any(b is not None for b in (bx, by, bz)):
            raise ValueError(
                "MHDModel: provide either (bx, by, bz) or (b1x, b1y, b1z), not both"
            )

        if b1_given:
            bx = self.defaulter(b1x, 0.0)
            by = self.defaulter(b1y, 0.0)
            bz = self.defaulter(b1z, 0.0)
        else:
            bx = self.defaulter(bx, 1.0)
            by = self.defaulter(by, 0.0)
            bz = self.defaulter(bz, 0.0)

        self.model_dict = {}

        self.model_dict.update(
            {
                "density": density,
                "vx": vx,
                "vy": vy,
                "vz": vz,
                "bx": bx,
                "by": by,
                "bz": bz,
                "p": p,
                "b_is_total": not b1_given,
            }
        )

        global_vars.sim.set_model(self)
