from pyphare import cpp

from . import global_vars
from .vector_potential import magnetic_non_periodic, resolve_magnetic_init


class MHDModel(object):
    """
    MHDModel sets the initial MHD state.

    **Parameters**:

        * **density** (*function*): mass density
        * **vx**, **vy**, **vz** (*function*): velocity components
        * **bx** (*function*): magnetic field in x direction
        * **by** (*function*): magnetic field in y direction
        * **bz** (*function*): magnetic field in z direction
        * **p** (*function*): pressure
        * **ax** (*function*): vector potential in x direction
        * **ay** (*function*): vector potential in y direction
        * **az** (*function*): vector potential in z direction

    Giving any of ``ax, ay, az`` sets B to the discrete curl of A on the Yee grid,
    so that the initial div B is at round-off. Missing ``a*`` are 0. In this mode,
    ``bx, by`` (and ``bz`` in 3D) cannot be given; in 2D, ``bz`` may be given
    instead of ``ax, ay``. A is not supported in 1D. A need not be periodic, only
    its curl does (e.g. a Harris sheet Az).
    """

    def defaulter(self, input, value):
        if input is not None:
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
        ax=None,
        ay=None,
        az=None,
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
        magnetic = resolve_magnetic_init(
            self.dim, self.defaulter, bx, by, bz, ax, ay, az
        )
        p = self.defaulter(p, 1.0)

        self.model_dict = {}

        self.model_dict.update(
            {
                "density": density,
                "vx": vx,
                "vy": vy,
                "vz": vz,
                "p": p,
            }
        )

        if magnetic["mode"] == "b":
            self.model_dict.update({k: magnetic[k] for k in ("bx", "by", "bz")})
        else:
            self.model_dict["vector_potential"] = {
                k: magnetic[k] for k in ("ax", "ay", "az")
            }
            self.model_dict["direct_b"] = magnetic["direct"]

        should_validate = not any(
            [global_vars.sim.dry_run, global_vars.sim.is_from_restart()]
        )
        self.validated = False
        if should_validate:
            self.validate(global_vars.sim)
            self.validated = True

        global_vars.sim.set_model(self)

    def validate(self, sim, atol=1e-15):
        """periodicity of the curl of A (A mode only)"""
        if "vector_potential" not in self.model_dict:
            return

        not_periodic = magnetic_non_periodic(sim, self.model_dict, atol)
        if not_periodic:
            cpp.print_rank0(
                "Warning: Simulation is periodic but some functions are not : ",
                not_periodic,
            )
            if sim.strict:
                raise RuntimeError("Simulation is not periodic")
