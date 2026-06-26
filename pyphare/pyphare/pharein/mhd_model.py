from . import global_vars


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
        self,
        density=None,
        vx=None,
        vy=None,
        vz=None,
        bx=None,
        by=None,
        bz=None,
        p=None,
        b0x=None,
        b0y=None,
        b0z=None,
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
        # b0x, b0y, b0z prescribe the static background field B0 (default zero, which reduces the
        # split formulation to classical MHD).
        b0x = self.defaulter(b0x, 0.0)
        b0y = self.defaulter(b0y, 0.0)
        b0z = self.defaulter(b0z, 0.0)
        # The grid stores the TOTAL field B = B0 + B1 under "bx/by/bz" (the C++ initializes B1 with
        # it, then subtracts B0). The user prescribes EITHER the total field directly (bx/by/bz) OR
        # the perturbation (b1x/b1y/b1z), in which case the total is B0 + B1.
        b1_given = any(b is not None for b in (b1x, b1y, b1z))
        if b1_given:
            if any(b is not None for b in (bx, by, bz)):
                raise ValueError("MHDModel: provide either (bx,by,bz) or (b1x,b1y,b1z), not both")
            b1x = self.defaulter(b1x, 0.0)
            b1y = self.defaulter(b1y, 0.0)
            b1z = self.defaulter(b1z, 0.0)
            bx = lambda *xyz: b0x(*xyz) + b1x(*xyz)
            by = lambda *xyz: b0y(*xyz) + b1y(*xyz)
            bz = lambda *xyz: b0z(*xyz) + b1z(*xyz)
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
                "b0x": b0x,
                "b0y": b0y,
                "b0z": b0z,
            }
        )

        global_vars.sim.set_model(self)
