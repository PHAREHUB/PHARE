import numpy as np

from . import global_vars
from .boundary import SuperMagnetofastInflowBoundary


class MHDModel(object):
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
        self, density=None, vx=None, vy=None, vz=None, bx=None, by=None, bz=None, p=None
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
        bx = self.defaulter(bx, 1.0)
        by = self.defaulter(by, 0.0)
        bz = self.defaulter(bz, 0.0)
        p = self.defaulter(p, 1.0)

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
            }
        )

        self.validate_inflow_normal_b(global_vars.sim)

        global_vars.sim.set_model(self)

    def validate_inflow_normal_b(self, sim, rtol=1e-6):
        domain = sim.simulation_domain()
        b_functions = [self.model_dict[name] for name in ("bx", "by", "bz")]

        for location, bc in sim.domain_boundaries.items():
            if not isinstance(bc, SuperMagnetofastInflowBoundary):
                continue

            normal = "xyz".index(location[0])
            face = 0.0 if location[1:] == "lower" else domain[normal]
            axes = [
                (
                    np.array([face])
                    if idir == normal
                    else (np.arange(sim.cells[idir]) + 0.5) * sim.dl[idir]
                )
                for idir in range(sim.ndim)
            ]
            coords = [c.ravel() for c in np.meshgrid(*axes, indexing="ij")]
            bn = np.broadcast_to(
                np.asarray(b_functions[normal](*coords), dtype=float), coords[0].shape
            )

            expected = bc.B[normal]
            scale = np.linalg.norm(bc.B) or 1.0
            deviation = np.abs(bn - expected)
            if not np.all(deviation <= rtol * scale):
                raise ValueError(
                    f"the initial {'xyz'[normal]} magnetic field on the inflow boundary "
                    f"'{location}' must equal the boundary value B[{normal}]={expected}, since "
                    f"the normal magnetic field stays at its initial value on an inflow face; "
                    f"max deviation {np.max(deviation)}"
                )
