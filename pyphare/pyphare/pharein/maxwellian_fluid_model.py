import numpy as np

from pyphare import cpp
from pyphare.core import phare_utilities
from pyphare.core.box import Box
from pyphare.core.gridlayout import HybridGridLayoutFor
from pyphare.pharein import global_vars


class MaxwellianFluidModel(object):
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
        * **init** (``dict``), optional ``{"seed": <int>}`` to fix the random seed used to draw this population's particle positions/velocities (default: unseeded, i.e. different particles on every run).

    All profile functions (`density`, `vbulk*`, `vth*`, and `bx`/`by`/`bz` above)
    take the spatial coordinates as arguments and must return a value (or an
    array of values, when called with array coordinates) of the same shape.

    **Periodicity check**: when the domain is periodic (the default), every
    profile function above is checked for periodicity at construction time.
    A profile that isn't actually periodic prints a warning; if
    ``Simulation(strict=True)`` was set, it raises instead. This check is
    skipped for `dry_run` simulations and when restarting.
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
        init={},
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

        init_keys = ["seed"]
        wrong_keys = phare_utilities.not_in_keywords_list(init_keys, **init)
        if len(wrong_keys) > 0:
            raise ValueError(
                "Model Error: invalid init arguments - " + " ".join(wrong_keys)
            )
        init["seed"] = init["seed"] if "seed" in init else None

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

    def validate(self, sim, atol=1e-15):
        phare_utilities.debug_print(f"validating dim={sim.ndim}")
        if sim.ndim == 1:
            self.validate1d(sim, atol)
        elif sim.ndim == 2:
            self.validate2d(sim, atol)
        elif sim.ndim == 3:
            self.validate3d(sim, atol)
        else:
            raise ValueError("Unknown dimension")

    def validate1d(self, sim, atol):
        domain_box = Box([0] * sim.ndim, sim.cells)
        layout = HybridGridLayoutFor(
            domain_box, domain_box.lower, sim.dl, sim.interp_order
        )
        nbrDualGhosts = layout.nbrGhostsPrimal(sim.interp_order)
        nbrPrimalGhosts = layout.nbrGhostsPrimal(sim.interp_order)
        directions = ["X"]
        domain = sim.simulation_domain()
        bx = self.model_dict["bx"]
        by = self.model_dict["by"]
        bz = self.model_dict["bz"]
        is_periodic = True

        dual_left = (np.arange(-nbrDualGhosts, nbrDualGhosts) + 0.5) * sim.dl[0]
        dual_right = dual_left + domain[0]
        primal_left = np.arange(-nbrPrimalGhosts, nbrPrimalGhosts) * sim.dl[0]
        primal_right = primal_left + domain[0]

        direction = directions[0]

        for b_i, b_name in zip((bx, by, bz), ("Bx", "By", "Bz")):
            if layout.qtyIsDual(b_name, direction):
                xL, xR = dual_left, dual_right
            else:
                xL, xR = primal_left, primal_right

            is_periodic &= np.allclose(b_i(xL), b_i(xL), atol=atol, rtol=0)

        for pop in self.populations:
            functions = ("vx", "vy", "vz", "vthx", "vthy", "vthz")
            xL, xR = primal_left, primal_right

            for fn in functions:
                f = self.model_dict[pop][fn]
                is_periodic &= np.allclose(f(xL), f(xR), atol=atol, rtol=0)

        if not is_periodic:
            cpp.print_rank0(
                "Warning: Simulation is periodic but some functions are not"
            )
            if sim.strict:
                raise RuntimeError("Simulation is not periodic")

    def validate2d(self, sim, atol):
        domain_box = Box([0] * sim.ndim, sim.cells)
        layout = HybridGridLayoutFor(
            domain_box, domain_box.lower, sim.dl, sim.interp_order
        )
        nbrDualGhosts = layout.nbrGhostsPrimal(sim.interp_order)
        nbrPrimalGhosts = layout.nbrGhostsPrimal(sim.interp_order)
        directions = ["X", "Y"]
        domain = sim.simulation_domain()
        bx = self.model_dict["bx"]
        by = self.model_dict["by"]
        bz = self.model_dict["bz"]
        is_periodic = True
        not_periodic = []

        def getCoord(L, R, idir):
            if idir == 0:
                return (L, np.zeros_like(L)), (R, np.zeros_like(R))
            else:
                return (np.zeros_like(L), L), (np.zeros_like(R), R)

        phare_utilities.debug_print("2d periodic validation")
        for idir in np.arange(sim.ndim):
            phare_utilities.debug_print("validating direction ...", idir)
            if sim.boundary_types[idir] == "periodic":
                phare_utilities.debug_print(f"direction {idir} is periodic?")
                dual_left = (np.arange(-nbrDualGhosts, nbrDualGhosts) + 0.5) * sim.dl[
                    idir
                ]
                dual_right = dual_left + domain[idir]
                primal_left = (
                    np.arange(-nbrPrimalGhosts, nbrPrimalGhosts) * sim.dl[idir]
                )
                primal_right = primal_left + domain[idir]

                direction = directions[idir]

                for b_i, b_name in zip((bx, by, bz), ("Bx", "By", "Bz")):
                    if layout.qtyIsDual(b_name, direction):
                        L, R = dual_left, dual_right
                    else:
                        L, R = primal_left, primal_right
                    coordsL, coordsR = getCoord(L, R, idir)
                    check = np.allclose(b_i(*coordsL), b_i(*coordsR), atol=atol, rtol=0)

                    if not check:
                        not_periodic += [(b_name, idir)]
                    is_periodic &= check

                for pop in self.populations:
                    functions = ("vx", "vy", "vz", "vthx", "vthy", "vthz")
                    L, R = primal_left, primal_right
                    coordsL, coordsR = getCoord(L, R, idir)

                    for fn in functions:
                        f = self.model_dict[pop][fn]
                        fL = f(*coordsL)
                        fR = f(*coordsR)
                        check = np.allclose(fL, fR, atol=atol, rtol=0)
                        phare_utilities.debug_print(
                            f"checked {fn} : fL = {fL} and fR = {fR} and check = {check}"
                        )
                        if not check:
                            not_periodic += [(fn, idir)]
                        is_periodic &= check

        if not is_periodic:
            phare_utilities.debug_print(
                "Warning: Simulation is periodic but some functions are not : ",
                not_periodic,
            )
            if sim.strict:
                raise RuntimeError("Simulation is not periodic")

    def validate3d(self, sim, atol):
        domain_box = Box([0] * sim.ndim, sim.cells)
        layout = HybridGridLayoutFor(
            domain_box, domain_box.lower, sim.dl, sim.interp_order
        )
        nbrDualGhosts = layout.nbrGhostsPrimal(sim.interp_order)
        nbrPrimalGhosts = layout.nbrGhostsPrimal(sim.interp_order)
        directions = ["X", "Y", "Z"]
        domain = sim.simulation_domain()
        bx = self.model_dict["bx"]
        by = self.model_dict["by"]
        bz = self.model_dict["bz"]
        is_periodic = True
        not_periodic = []

        def getCoord(L, R, idir):
            if idir == 0:
                return (L, np.zeros_like(L), np.zeros_like(L)), (
                    R,
                    np.zeros_like(R),
                    np.zeros_like(R),
                )
            elif idir == 1:
                return (np.zeros_like(L), L, np.zeros_like(L)), (
                    np.zeros_like(R),
                    R,
                    np.zeros_like(R),
                )
            else:
                return (np.zeros_like(L), np.zeros_like(L), L), (
                    np.zeros_like(R),
                    np.zeros_like(R),
                    R,
                )

        phare_utilities.debug_print("3d periodic validation")
        for idir in np.arange(sim.ndim):
            phare_utilities.debug_print("validating direction ...", idir)
            if sim.boundary_types[idir] == "periodic":
                phare_utilities.debug_print(f"direction {idir} is periodic?")
                dual_left = (np.arange(-nbrDualGhosts, nbrDualGhosts) + 0.5) * sim.dl[
                    idir
                ]
                dual_right = dual_left + domain[idir]
                primal_left = (
                    np.arange(-nbrPrimalGhosts, nbrPrimalGhosts) * sim.dl[idir]
                )
                primal_right = primal_left + domain[idir]

                direction = directions[idir]

                for b_i, b_name in zip((bx, by, bz), ("Bx", "By", "Bz")):
                    if layout.qtyIsDual(b_name, direction):
                        L, R = dual_left, dual_right
                    else:
                        L, R = primal_left, primal_right

                    coordsL, coordsR = getCoord(L, R, idir)
                    check = np.allclose(b_i(*coordsL), b_i(*coordsR), atol=atol, rtol=0)
                    if not check:
                        not_periodic += [(b_name, idir)]
                    is_periodic &= check

                for pop in self.populations:
                    functions = ("vx", "vy", "vz", "vthx", "vthy", "vthz")
                    L, R = primal_left, primal_right
                    coordsL, coordsR = getCoord(L, R, idir)

                    for fn in functions:
                        f = self.model_dict[pop][fn]
                        fL = f(*coordsL)
                        fR = f(*coordsR)
                        check = np.allclose(fL, fR, atol=atol, rtol=0)
                        phare_utilities.debug_print(
                            f"checked {fn} : fL = {fL} and fR = {fR} and check = {check}"
                        )
                        if not check:
                            not_periodic += [(fn, idir)]
                        is_periodic &= check

        if not is_periodic:
            cpp.print_rank0(
                "Warning: Simulation is periodic but some functions are not : ",
                not_periodic,
            )
            if sim.strict:
                raise RuntimeError("Simulation is not periodic")
