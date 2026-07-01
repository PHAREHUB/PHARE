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
        ax=None,
        ay=None,
        az=None,
        a0x=None,
        a0y=None,
        a0z=None,
        a1x=None,
        a1y=None,
        a1z=None,
    ):
        if global_vars.sim is None:
            raise RuntimeError("A simulation must be declared before a model")

        if global_vars.sim.model is not None:
            raise RuntimeError("A model is already created")

        self.dim = global_vars.sim.ndim

        # --- vector-potential init -----------------------------------------------
        # B = curl(A) with a full 3D vector potential A = (Ax, Ay, Az), computed on the C++ side
        # with the discrete curl so that div B = 0 to machine precision. This mirrors the direct-B
        # interface:
        #   a0x/a0y/a0z -> B0 = curl(a0)   (background; exclusive with b0x/b0y/b0z)
        #   a1x/a1y/a1z -> B1 = curl(a1)   (perturbation; exclusive with bx/by/bz and b1x/b1y/b1z)
        #   ax/ay/az    -> plain alias for the perturbation with B0 = 0 (folded into a1; exclusive
        #                  with a0*/a1*)
        # Either, both, or neither of a0/a1 may be given (independent modes), with the component-
        # wise init as the default fallback. The legacy 2D scalar keys a0z/a1z are just the
        # z-component (x/y default to zero), so pre-existing 2D scripts are unchanged.
        a_given = any(a is not None for a in (ax, ay, az))
        a0_given = any(a is not None for a in (a0x, a0y, a0z))
        a1_given = any(a is not None for a in (a1x, a1y, a1z))
        if a_given and (a0_given or a1_given):
            raise ValueError(
                "MHDModel: plain (ax,ay,az) is exclusive with a0*/a1* (it is the perturbation "
                "potential with B0 = 0)"
            )
        if a_given:
            # Plain a is the perturbation potential with no background: fold it into a1, keep B0 = 0.
            a1x, a1y, a1z = ax, ay, az
            a1_given = True

        b0_from_potential = a0_given
        b1_from_potential = a1_given

        density = self.defaulter(density, 1.0)
        vx = self.defaulter(vx, 1.0)
        vy = self.defaulter(vy, 0.0)
        vz = self.defaulter(vz, 0.0)
        p = self.defaulter(p, 1.0)
        # b0x, b0y, b0z prescribe the static background field B0 (default zero, which reduces the
        # split formulation to classical MHD). Capture whether a component B0 was given before it
        # is defaulted, so a b0-only run keeps B1 = 0 (see the field-handling below).
        b0_given = any(b is not None for b in (b0x, b0y, b0z))
        if b0_from_potential and b0_given:
            raise ValueError(
                "MHDModel: a0 (B0 from vector potential) is exclusive with b0x/b0y/b0z"
            )
        b0x = self.defaulter(b0x, 0.0)
        b0y = self.defaulter(b0y, 0.0)
        b0z = self.defaulter(b0z, 0.0)
        a0x = self.defaulter(a0x, 0.0)
        a0y = self.defaulter(a0y, 0.0)
        a0z = self.defaulter(a0z, 0.0)
        a1x = self.defaulter(a1x, 0.0)
        a1y = self.defaulter(a1y, 0.0)
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
                "MHDModel: a1 (B1 from vector potential) is exclusive with bx/by/bz and b1x/b1y/b1z"
            )
        if b1_from_potential:
            # B1 is built on the C++ side from a1; the "magnetic" dict is unused for B1 but its
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
        elif b0_from_potential:
            # B0 from a potential: the stored field is B1 directly (no subtraction); an unspecified
            # field means no perturbation.
            bx = self.defaulter(bx, 0.0)
            by = self.defaulter(by, 0.0)
            bz = self.defaulter(bz, 0.0)
        elif b0_given and not b_total_given:
            # Only a component background B0 was given (no total bx/by/bz): the total equals B0, so
            # B1 = total - B0 = 0. (Falling through to the classical default below would set bx=1
            # and spuriously make B1 = (1,0,0) - B0.)
            bx, by, bz = b0x, b0y, b0z
        else:
            # Classical default total field: uniform Bx = 1.
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
                "a0x": a0x,
                "a0y": a0y,
                "a0z": a0z,
                "a1x": a1x,
                "a1y": a1y,
                "a1z": a1z,
                "b0_init_mode": "potential" if b0_from_potential else "components",
                "b1_init_mode": "potential" if b1_from_potential else "components",
            }
        )

        global_vars.sim.set_model(self)
