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
        a0z=None,
        a1z=None,
    ):
        if global_vars.sim is None:
            raise RuntimeError("A simulation must be declared before a model")

        if global_vars.sim.model is not None:
            raise RuntimeError("A model is already created")

        self.dim = global_vars.sim.ndim

        # --- vector-potential init (2D only) -------------------------------------
        # B = curl(A_z z_hat): Bx = dA_z/dy, By = -dA_z/dx, Bz = 0. Computed on the C++ side with
        # the discrete curl so that div B = 0 to machine precision. a0z drives B0, a1z drives B1;
        # either, both, or neither may be given (independent modes), with the component-wise init
        # as the default fallback.
        b0_from_potential = a0z is not None
        b1_from_potential = a1z is not None
        if (b0_from_potential or b1_from_potential) and self.dim != 2:
            raise ValueError(
                "MHDModel vector-potential init (a0z/a1z) is only supported in 2D"
            )
        if b0_from_potential and any(b is not None for b in (b0x, b0y)):
            raise ValueError(
                "MHDModel: a0z (B0 from vector potential) is exclusive with b0x/b0y"
            )

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
        a0z = self.defaulter(a0z, 0.0)
        a1z = self.defaulter(a1z, 0.0)
        # The grid stores the TOTAL field B = B0 + B1 under "bx/by/bz" (the C++ initializes B1 with
        # it, then subtracts B0). The user prescribes EITHER the total field directly (bx/by/bz) OR
        # the perturbation (b1x/b1y/b1z), in which case the total is B0 + B1. When B0 comes from a
        # potential the Python total cannot fold it in, so it stores B1 directly and the C++ skips
        # the subtraction.
        b1_given = any(b is not None for b in (b1x, b1y, b1z))
        b_total_given = any(b is not None for b in (bx, by, bz))
        if b1_from_potential and (b1_given or b_total_given):
            raise ValueError(
                "MHDModel: a1z (B1 from vector potential) is exclusive with bx/by/bz and b1x/b1y/b1z"
            )
        if b1_from_potential:
            # B1 is built on the C++ side from a1z; the "magnetic" dict is unused for B1 but its
            # keys must exist, so fill them with zeros.
            bx = self.defaulter(None, 0.0)
            by = self.defaulter(None, 0.0)
            bz = self.defaulter(None, 0.0)
        elif b1_given:
            if b_total_given:
                raise ValueError("MHDModel: provide either (bx,by,bz) or (b1x,b1y,b1z), not both")
            b1x = self.defaulter(b1x, 0.0)
            b1y = self.defaulter(b1y, 0.0)
            b1z = self.defaulter(b1z, 0.0)
            bx = lambda *xyz: b0x(*xyz) + b1x(*xyz)
            by = lambda *xyz: b0y(*xyz) + b1y(*xyz)
            bz = lambda *xyz: b0z(*xyz) + b1z(*xyz)
        else:
            # When B0 comes from a potential the stored field is B1 (no subtraction), so an
            # unspecified field means no perturbation (default 0); otherwise it is the classical
            # total field whose default is a uniform Bx = 1.
            bx = self.defaulter(bx, 0.0 if b0_from_potential else 1.0)
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
                "a0z": a0z,
                "a1z": a1z,
                "b0_init_mode": "potential" if b0_from_potential else "components",
                "b1_init_mode": "potential" if b1_from_potential else "components",
            }
        )

        global_vars.sim.set_model(self)
