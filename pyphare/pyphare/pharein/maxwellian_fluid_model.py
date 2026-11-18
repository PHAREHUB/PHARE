import numbers

import numpy as np

from pyphare import cpp
from pyphare.core import phare_utilities
from pyphare.core.box import Box
from pyphare.core.gridlayout import HybridGridLayoutFor
from pyphare.pharein import global_vars


class MaxwellianFluidModel:
    """
    MaxwellianFluidModel sets up the Hybrid simulation type: the magnetic field
    profile and one or more ion populations, each assumed to locally follow a
    Maxwellian velocity distribution described by its density, bulk velocity and
    thermal velocity moments. Every population is loaded as macro-particles
    from these moments; the resulting particles are pushed kinetically while
    the field is advanced from the fluid electron/Ohm's law closure set by
    :class:`~pyphare.pharein.ElectronModel`.

    Used for a Hybrid simulation (the default `model_options`). For an MHD
    simulation, use :class:`~pyphare.pharein.MHDModel` instead.

    **Usage example:**

    .. code-block:: python

        ph.MaxwellianFluidModel(
            bx=bx,
            by=by,
            bz=bz,
            protons={"charge": 1,
                     "density": density,
                     "vbulkx": vx,
                     "vbulky": vy,
                     "vbulkz": vz,
                     "vthx": vthx,
                     "vthy": vthy,
                     "vthz": vthz,
                     "nbr_part_per_cell": 100},

            beam={"charge": 1,
                     "density": beam_density,
                     "vbulkx": vx_beam,
                     "vbulky": vy_beam,
                     "vbulkz": vz_beam,
                     "vthx": vthx_beam,
                     "vthy": vthy_beam,
                     "vthz": vthz_beam,
                     "nbr_part_per_cell": 500}
        )

    **Magnetic field parameters**:

        * **bx** (*function*): magnetic field in x direction
        * **by** (*function*): magnetic field in y direction
        * **bz** (*function*): magnetic field in z direction

    Each of `bx`, `by`, `bz` is a function of the spatial coordinates (one
    argument per dimension: `f(x)`, `f(x, y)` or `f(x, y, z)`), or a plain
    number for a spatially uniform component (default: `bx=1`, `by=bz=0`).

    **Ion populations**:

    Every other keyword argument names one ion population; its value is a
    dict of that population's parameters, all accepted by `add_population`:

        * **charge** (``float``), default=1, electric charge of one particle of this population.
        * **mass** (``float``), default=1, mass of one particle of this population.
        * **density** (*function or number*), default=1, particle density profile.
        * **vbulkx**, **vbulky**, **vbulkz** (*function or number*), default=0, bulk velocity profile, per direction.
        * **vthx**, **vthy**, **vthz** (*function or number*), default=1, thermal velocity profile, per direction.
        * **nbr_part_per_cell** (``int``), default=100, number of macro-particles loaded per cell for this population.
        * **density_cut_off** (``float``), default=1e-16, cells where `density` evaluates below this value are not populated with particles.
        * **init** (``dict``), optional ``{"seed": <int>}`` to fix the random seed used to draw this population's particle positions/velocities (default: unseeded, i.e. different particles on every run). Each cell draws from its own stream, keyed by the seed, the population name and the cell, so a seeded run gives the same particles whatever the patch layout or number of MPI ranks (for a given compiler and standard library). Two populations with the same seed still get different particles, and renaming a population changes its particles.

    All profile functions (`density`, `vbulk*`, `vth*`, and `bx`/`by`/`bz` above)
    take the spatial coordinates as arguments and must return a value (or an
    array of values, when called with array coordinates) of the same shape.

    **Periodicity check**: when the domain is periodic (the default), every
    profile function above is checked for periodicity at construction time.
    A profile that isn't actually periodic prints a warning; if
    ``Simulation(strict=True)`` was set, it raises instead. This check is
    skipped for `dry_run` (see :class:`~pyphare.pharein.Simulation`)
    simulations and when restarting.
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

    def __init__(self, bx=None, by=None, bz=None, **kwargs):
        if global_vars.sim is None:
            raise RuntimeError("A simulation must be declared before a model")

        if global_vars.sim.model is not None:
            raise RuntimeError("A model is already created")

        self.dim = global_vars.sim.ndim
        bx = self.defaulter(bx, 1.0)
        by = self.defaulter(by, 0.0)
        bz = self.defaulter(bz, 0.0)

        self.model_dict = {"model": "model", "model_name": "custom"}

        self.model_dict.update({"bx": bx, "by": by, "bz": bz})

        self.populations = list(kwargs.keys())
        for population in self.populations:
            self.add_population(population, **kwargs[population])

        should_validate = not any(
            [global_vars.sim.dry_run, global_vars.sim.is_from_restart()]
        )
        self.validated = False
        if should_validate:
            # on every rank so a strict failure stops them all, rather than rank 0 alone
            self.validate(global_vars.sim)
            self.validated = True

        global_vars.sim.set_model(self)

    # ------------------------------------------------------------------------------

    def nbr_populations(self):
        """
        returns the number of species currently registered in the model
        """
        return len(self.populations)

    # ------------------------------------------------------------------------------

    def add_population(
        self,
        name,
        charge=1.0,
        mass=1.0,
        nbr_part_per_cell=100,
        density=None,
        vbulkx=None,
        vbulky=None,
        vbulkz=None,
        vthx=None,
        vthy=None,
        vthz=None,
        init=None,
        density_cut_off=1e-16,
    ):
        """
        Add an ion population, named `name`, to this model.

        This is what the population dicts passed to
        :class:`MaxwellianFluidModel`'s constructor (e.g. `protons={...}`) are
        forwarded to; see the class docstring for the full parameter list
        (`charge`, `mass`, `density`, `vbulkx/y/z`, `vthx/y/z`,
        `nbr_part_per_cell`, `density_cut_off`, `init`).
        """

        init = init if init is not None else {}

        init_keys = ["seed"]
        wrong_keys = phare_utilities.not_in_keywords_list(init_keys, **init)
        if len(wrong_keys) > 0:
            raise ValueError(
                "Model Error: invalid init arguments - " + " ".join(wrong_keys)
            )
        init["seed"] = init.get("seed", None)

        density = self.defaulter(density, 1.0)

        vbulkx = self.defaulter(vbulkx, 0.0)
        vbulky = self.defaulter(vbulky, 0.0)
        vbulkz = self.defaulter(vbulkz, 0.0)

        vthx = self.defaulter(vthx, 1.0)
        vthy = self.defaulter(vthy, 1.0)
        vthz = self.defaulter(vthz, 1.0)

        new_population = {
            name: {
                "charge": charge,
                "mass": mass,
                "density": density,
                "vx": vbulkx,
                "vy": vbulky,
                "vz": vbulkz,
                "vthx": vthx,
                "vthy": vthy,
                "vthz": vthz,
                "nbrParticlesPerCell": nbr_part_per_cell,
                "init": init,
                "density_cut_off": density_cut_off,
            }
        }

        keys = self.model_dict.keys()
        if name in keys:
            raise ValueError("population already registered")

        self.model_dict.update(new_population)

    # ------------------------------------------------------------------------------

    def to_dict(self):
        self.model_dict["nbr_ion_populations"] = self.nbr_populations()
        return self.model_dict

    # ------------------------------------------------------------------------------

    def validate(self, sim, rtol=1e-12, transverse_samples=64):
        """
        Warn (raise under sim.strict) if a user function is not periodic along a periodic
        direction: f at the nodes around the lower border, within the ghost width, against f
        one domain length further, on the whole transverse plane (every primal and dual node in
        2D, a grid of at most transverse_samples^2 of them in 3D). A jump larger than rtol times
        the largest |value| in its group, over the domain (B, a population's density, a population's bulk
        and thermal velocities) is reported with the function, the direction and where it is.
        A periodic node shared by two patches is evaluated on both sides, so such a jump makes
        their copies differ from t = 0.
        """
        if sim.ndim not in (1, 2, 3):
            raise ValueError("Unknown dimension")
        domain_box = Box([0] * sim.ndim, sim.cells)
        layout = HybridGridLayoutFor(
            domain_box, domain_box.lower, sim.dl, sim.interp_order
        )
        nbrGhosts = layout.nbrGhostsPrimal(sim.interp_order)
        directions = ["X", "Y", "Z"]
        domain = sim.simulation_domain()

        def transverse(idir):
            nodes = np.arange(0, sim.cells[idir] + 1) * sim.dl[idir]
            points = np.sort(np.concatenate((nodes, nodes[:-1] + 0.5 * sim.dl[idir])))
            if sim.ndim == 3 and points.size > transverse_samples:
                points = points[:: -(-points.size // transverse_samples)]
            return points

        # (function, name, group): a jump is relative to the largest |f| of its group, so
        # that a component close to 0 everywhere is not judged on its own round-off
        functions = [
            (self.model_dict[b], b.capitalize(), "B") for b in ("bx", "by", "bz")
        ]
        for pop in self.populations:
            functions += [(self.model_dict[pop]["density"], f"{pop} density", pop)]
            functions += [
                (self.model_dict[pop][fn], f"{pop} {fn}", f"{pop} velocity")
                for fn in ("vx", "vy", "vz", "vthx", "vthy", "vthz")
            ]

        # the size of each group over the domain, on at most 256 nodes per direction
        everywhere = np.meshgrid(
            *[transverse(d)[:: -(-transverse(d).size // 256)] for d in range(sim.ndim)],
            indexing="ij",
        )
        scales = {}
        for f, _, group in functions:
            scales[group] = max(scales.get(group, 0.0), np.abs(f(*everywhere)).max())

        not_periodic = []
        for idir in range(sim.ndim):
            if not sim.periodicities[idir]:
                continue
            jumps = []
            for f, name, group in functions:
                qty = name if name in ("Bx", "By", "Bz") else "rho"
                border = np.arange(-nbrGhosts, nbrGhosts) * sim.dl[idir]
                if layout.qtyIsDual(qty, directions[idir]):
                    border = border + 0.5 * sim.dl[idir]
                axes = [border if d == idir else transverse(d) for d in range(sim.ndim)]
                left = np.meshgrid(*axes, indexing="ij")
                right = [
                    c + domain[idir] if d == idir else c for d, c in enumerate(left)
                ]
                fL = np.broadcast_to(f(*left), left[0].shape)
                fR = np.broadcast_to(f(*right), left[0].shape)
                jump = np.abs(fL - fR)
                at = np.unravel_index(jump.argmax(), jump.shape)
                jumps += [(name, group, jump[at], tuple(float(c[at]) for c in left))]
                scales[group] = max(scales[group], np.abs(fL).max(), np.abs(fR).max())
            not_periodic += [
                f"{name} along {directions[idir]}: jump {jump:.3e} (largest |{group}|"
                f" {scales[group]:.3e}) between {where} and its image"
                for name, group, jump, where in jumps
                if jump > rtol * scales[group]
            ]

        if not_periodic:
            cpp.print_rank0(
                "Warning: Simulation is periodic but some functions are not:\n  "
                + "\n  ".join(not_periodic)
            )
            if sim.strict:
                raise RuntimeError("Simulation is not periodic")
