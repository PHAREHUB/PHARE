import numpy as np

from ..core import phare_utilities
from . import global_vars


def all_timestamps(sim):
    init_time = sim.start_time()
    # consume the authoritative integer step count rather than re-deriving it from
    # final_time/dt (fp-fragile). last stamp = init + nbr*dt = sim.final_time exactly.
    return sim.time_step * np.arange(sim.time_step_nbr + 1) + init_time


# ------------------------------------------------------------------------------


def diagnostics_checker(func):
    def wrapper(diagnostics_object, name, **kwargs):
        mandatory_keywords = ["write_timestamps", "quantity"]

        # check if some mandatory keywords are not missing
        missing_mandatory_kwds = phare_utilities.check_mandatory_keywords(
            mandatory_keywords, **kwargs
        )
        if len(missing_mandatory_kwds) > 0:
            raise RuntimeError(
                "Error: missing mandatory parameters : "
                + ", ".join(missing_mandatory_kwds)
            )

        one_of_required = ["elapsed_timestamps", "write_timestamps"]
        if not any([k in kwargs for k in one_of_required]):
            raise RuntimeError(
                "Error: missing parameters - one required: "
                + ", ".join(one_of_required)
            )

        accepted_keywords = ["path", "population_name", "flush_every"]
        accepted_keywords += mandatory_keywords

        # check that all passed keywords are in the accepted keyword list
        wrong_kwds = phare_utilities.not_in_keywords_list(accepted_keywords, **kwargs)
        if len(wrong_kwds) > 0:
            raise RuntimeError("Error: invalid arguments - " + " ".join(wrong_kwds))

        try:
            # just take mandatory arguments from the dict
            # since if we arrived here we are sure they are there

            kwargs["path"] = kwargs.get("path", "./")

            return func(diagnostics_object, name, **kwargs)

        except ValueError as msg:
            print(msg)

    return wrapper


# ------------------------------------------------------------------------------
def validate_timestamps(clazz, key, **kwargs):
    sim = global_vars.sim

    init_time = sim.start_time()

    timestamps = phare_utilities.np_array_ify(kwargs.get(key, []))

    if np.any(timestamps < init_time):
        timestamps = timestamps[timestamps >= init_time]
        print(f"Warning: some timestamps below ({init_time}) are filtered")

    if np.any(timestamps > sim.final_time):
        raise RuntimeError(
            f"Error: timestamp({sim.time_step_nbr}) cannot be greater than simulation.final_time({sim.final_time}))"
        )
    if not np.all(np.diff(timestamps) >= 0):
        raise RuntimeError(f"Error: {clazz}.{key} not in ascending order)")
    if not np.all(
        np.abs(timestamps / sim.time_step - np.rint(timestamps / sim.time_step) < 1e-9)
    ):
        raise RuntimeError(
            f"Error: {clazz}.{key} is inconsistent with simulation.time_step"
        )

    return timestamps


def validate_elapsed_timestamps(clazz, key, **kwargs):
    import datetime

    timestamps = phare_utilities.np_array_ify(kwargs.get(key, []))

    timestamps = [
        int(ts.total_seconds()) if isinstance(ts, datetime.timedelta) else ts
        for ts in phare_utilities.np_array_ify(timestamps)
    ]

    if not np.all(np.diff(timestamps) >= 0):
        raise RuntimeError(f"Error: {clazz}.{key} not in ascending order)")

    return timestamps


# ------------------------------------------------------------------------------


def try_cpp_dep_vers():
    try:
        from pyphare.cpp import cpp_etc_lib

        return cpp_etc_lib().phare_deps()
    except ImportError:
        return {}


def try_cpp_build_config():
    try:
        from pyphare import cpp

        return cpp.build_config()
    except ImportError:
        return {}


class Diagnostics(object):
    """
    Base parameters shared by every diagnostics block below
    (:class:`ElectromagDiagnostics`, :class:`FluidDiagnostics`,
    :class:`ParticleDiagnostics`, :class:`MHDDiagnostics`,
    :class:`MetaDiagnostics`, :class:`InfoDiagnostics`). A diagnostics block
    writes one physical `quantity` to disk at a chosen set of times.

    **Common parameters**:

        * **quantity** (``str``), mandatory, which quantity to write - valid values depend on the diagnostics type, see each class below.
        * **write_timestamps** (``list`` of ``float``), simulation times at which to write this quantity. Required unless `elapsed_timestamps` is given instead.
        * **elapsed_timestamps** (``list`` of ``float`` seconds or ``datetime.timedelta``), write this quantity every time this much wall-clock time has elapsed since the simulation started, instead of at fixed simulation times. Can be given instead of, or in addition to, `write_timestamps`.
        * **path** (``str``), default="./", directory the diagnostics files for this block are written to.
        * **flush_every** (``int``), default=1, write to disk every `flush_every` dump (1 flushes every time, safest but slowest; larger values buffer more dumps in memory before writing, trading safety against a crash for I/O throughput).

    Declaring two diagnostics blocks of the same type with the same
    `quantity` (and, where relevant, the same `population_name`) does not
    create two independent diagnostics: it extends the existing one's list
    of timestamps with the new ones instead.
    """

    h5_flush_never = 0
    cpp_dep_vers = try_cpp_dep_vers()

    @diagnostics_checker
    def __init__(self, name, **kwargs):
        if global_vars.sim is None:
            raise RuntimeError("A simulation must be created before adding diagnostics")

        self.name = name
        self.path = kwargs["path"]

        self.write_timestamps = validate_timestamps(
            self.__class__.__name__, "write_timestamps", **kwargs
        )

        self.elapsed_timestamps = validate_elapsed_timestamps(
            self.__class__.__name__, "elapsed_timestamps", **kwargs
        )
        # for now every diagnostics needing computation (like momentum tensor)
        # is computing at the same time as written.
        # later this parameter can evolve to allow for different timestamps
        # that will depend on the type of diagnostics.
        self.compute_timestamps = self.write_timestamps

        self.attributes = kwargs.get("attributes", {})

        build_config = try_cpp_build_config()
        self.attributes["git_hash"] = build_config.get(
            "GIT_HASH", "git hash not available"
        )

        for dep, dep_ver in Diagnostics.cpp_dep_vers.items():
            self.attributes[f"{dep}_version"] = dep_ver

        for key in self.attributes:
            self.attributes[key] = self.attributes[key]

        self.quantity = None  # set in next line, stops pylint complaining
        self._setSubTypeAttributes(**kwargs)
        self.flush_every = kwargs.get(
            "flush_every", 1
        )  # flushes every dump, safe, but costly

        if self.flush_every < 0:
            raise RuntimeError(
                f"{(self.__class__.__name__,)}.flush_every cannot be negative"
            )

        self.__extent = None

        # if a diag already is registered we just concatenate the timestamps
        addIt = True
        registered_diags = global_vars.sim.diagnostics
        for diagname, diag in registered_diags.items():
            if self.quantity == diag.quantity:
                print(
                    f"{diag.name} already registered {self.quantity}, merging timestamps"
                )
                my_times = self.write_timestamps
                existing_times = diag.write_timestamps
                new_times = np.concatenate((my_times, existing_times))
                new_times.sort()
                mask = np.ones(len(new_times), dtype=bool)
                mask[1:] = (
                    np.diff(new_times) > 1e-10
                )  # assumed smaller than any realistic dt
                global_vars.sim.diagnostics[diagname].write_timestamps = new_times[mask]
                addIt = False
                break  # there can be only one

        if addIt:
            global_vars.sim.add_diagnostics(self)

    def extent(self):
        return self.__extent

    def _setSubTypeAttributes(self, **kwargs):  # stop pyline complaining
        raise RuntimeError("Never to be called, defined in subclass")

    def to_dict(self, type, **kwargs):
        return {
            "name": self.name,
            "type": type,
            "quantity": self.quantity,
            "write_timestamps": self.write_timestamps,
            "elapsed_timestamps": self.elapsed_timestamps,
            "compute_timestamps": self.compute_timestamps,
            "path": self.path,
            **kwargs,
        }


# ------------------------------------------------------------------------------
class MHDDiagnostics(Diagnostics):
    """
    Writes MHD fluid quantities to disk. Only meaningful in an MHD
    simulation (``Simulation(model_options=["MHDModel"])``).

    **Usage example:**

    .. code-block:: python

        from pyphare.pharein import MHDDiagnostics

        MHDDiagnostics(quantity="rho", write_timestamps=timestamps)

    **quantity**, one of:

        * "rho" - mass density
        * "V" - velocity
        * "P" - thermal pressure
        * "rhoV" - momentum density
        * "Etot" - total energy density
    """

    mhd_quantities = ["rho", "V", "P", "rhoV", "Etot"]
    type = "mhd"

    def __init__(self, **kwargs):
        super(MHDDiagnostics, self).__init__(
            MHDDiagnostics.type
            + str(global_vars.sim.count_diagnostics(MHDDiagnostics.type)),
            **kwargs,
        )

    def _setSubTypeAttributes(self, **kwargs):
        if kwargs["quantity"] not in MHDDiagnostics.mhd_quantities:
            error_msg = "Error: '{}' not a valid mhd diagnostics : " + ", ".join(
                MHDDiagnostics.mhd_quantities
            )
            raise ValueError(error_msg.format(kwargs["quantity"]))
        else:
            self.quantity = "/mhd/" + kwargs["quantity"]

        self.attributes["heat_capacity_ratio"] = global_vars.sim.gamma

    def to_dict(self):
        return {
            "name": self.name,
            "type": MHDDiagnostics.type,
            "quantity": self.quantity,
            "write_timestamps": self.write_timestamps,
            "compute_timestamps": self.compute_timestamps,
            "path": self.path,
        }


# ------------------------------------------------------------------------------
class ElectromagDiagnostics(Diagnostics):
    """
    Writes the electric or magnetic field to disk.

    **Usage example:**

    .. code-block:: python

        from pyphare.pharein import ElectromagDiagnostics

        ElectromagDiagnostics(quantity="E", write_timestamps=timestamps)
        ElectromagDiagnostics(quantity="B", write_timestamps=timestamps)

    **quantity**, one of: "E" (electric field), "B" (magnetic field).
    """

    em_quantities = ["E", "B"]
    type = "electromag"

    def __init__(self, **kwargs):
        super(ElectromagDiagnostics, self).__init__(
            ElectromagDiagnostics.type
            + str(global_vars.sim.count_diagnostics(ElectromagDiagnostics.type)),
            **kwargs,
        )

    def _setSubTypeAttributes(self, **kwargs):
        if kwargs["quantity"] not in ElectromagDiagnostics.em_quantities:
            error_msg = "Error: '{}' not a valid electromag diagnostics : " + ", ".join(
                ElectromagDiagnostics.em_quantities
            )
            raise ValueError(error_msg.format(kwargs["quantity"]))
        else:
            self.quantity = "/EM_" + kwargs["quantity"]

    def to_dict(self):
        return super().to_dict(type(self).type)


# ------------------------------------------------------------------------------


def population_in_model(population):
    return population in [p for p in global_vars.sim.model.populations]


class FluidDiagnostics_(Diagnostics):
    """:meta private:"""

    fluid_quantities = [
        "density",
        "charge_density",
        "mass_density",
        "flux",
        "bulkVelocity",
        "momentum_tensor",
    ]
    type = "fluid"

    def __init__(self, **kwargs):
        super(FluidDiagnostics_, self).__init__(
            FluidDiagnostics_.type
            + str(global_vars.sim.count_diagnostics(FluidDiagnostics_.type)),
            **kwargs,
        )

    def _setSubTypeAttributes(self, **kwargs):
        self.population_name = None
        if "population_name" not in kwargs and kwargs["quantity"] == "flux":
            raise ValueError("Error: missing population_name")
        elif "population_name" in kwargs:
            self.population_name = kwargs["population_name"]

        if "population_name" not in kwargs and kwargs["quantity"] == "density":
            raise ValueError("Error: cannot use density without population name")

        if kwargs["quantity"] not in FluidDiagnostics_.fluid_quantities:
            error_msg = "Error: '{}' not a valid fluid diagnostics : " + ", ".join(
                FluidDiagnostics_.fluid_quantities
            )
            raise ValueError(error_msg.format(kwargs["quantity"]))
        elif kwargs["quantity"] == "flux" and kwargs["population_name"] == "ions":
            raise ValueError(
                "'flux' is only available for specific populations, try 'bulkVelocity"
            )
        elif kwargs["quantity"] == "pressure_tensor":
            raise ValueError("'pressure_tensor' is not available yet")
        else:
            self.quantity = kwargs["quantity"]

        if self.population_name is None:
            self.quantity = "/ions/" + self.quantity
        else:
            if not population_in_model(self.population_name):
                raise ValueError(
                    "Error: population '{}' not in simulation initial model".format(
                        self.population_name
                    )
                )
            self.quantity = "/ions/pop/" + self.population_name + "/" + self.quantity

    def to_dict(self):
        return super().to_dict(type(self).type, population_name=self.population_name)


def for_total_ions(**kwargs):
    return "population_name" not in kwargs


class FluidDiagnostics:
    """
    Writes ion fluid moments (density, velocity, pressure, ...) to disk.
    These are moments of the particle distribution(s), either for one named
    population or, for some quantities, for the ions as a whole.

    **Usage example:**

    .. code-block:: python

        from pyphare.pharein import FluidDiagnostics

        # a per-population quantity: population_name required
        FluidDiagnostics(
            quantity="density",
            population_name="protons",
            write_timestamps=timestamps,
        )

        # a quantity that also works for the ions as a whole: population_name omitted
        FluidDiagnostics(quantity="mass_density", write_timestamps=timestamps)

    **quantity**, always requires **population_name**:

        * "density" - particle density of the named population
        * "flux" - flux of the named population

    **quantity**, works either for the ions as a whole (omit
    `population_name`) or for one named population (give `population_name`):

        * "charge_density" - charge density
        * "mass_density" - mass density
        * "bulkVelocity" - bulk velocity
        * "momentum_tensor" - momentum tensor
        * "pressure_tensor" - pressure tensor. This is a convenience quantity:
          declaring it registers the underlying moments it is computed from
          (`mass_density`, `bulkVelocity` and `momentum_tensor` for the ions
          as a whole; `density`, `flux` and `momentum_tensor` for a named
          population) as separate diagnostics under the hood.

    **population_name** (``str``), name of the ion population this
    diagnostics applies to, as declared in
    :class:`~pyphare.pharein.MaxwellianFluidModel`.
    """

    def __init__(self, **kwargs):
        if kwargs["quantity"] == "pressure_tensor":
            if for_total_ions(**kwargs):
                needed_quantities = ["mass_density", "bulkVelocity", "momentum_tensor"]
            else:
                needed_quantities = ["density", "flux", "momentum_tensor"]

            for quantity in needed_quantities:
                kwargs["quantity"] = quantity
                FluidDiagnostics_(**kwargs)
        else:
            FluidDiagnostics_(**kwargs)


# ------------------------------------------------------------------------------


class ParticleDiagnostics(Diagnostics):
    """
    Writes raw macro-particles of one ion population to disk. Typically much
    heavier (in disk size and I/O time) than the other diagnostics types.

    **Usage example:**

    .. code-block:: python

        from pyphare.pharein import ParticleDiagnostics

        ParticleDiagnostics(
            quantity="domain",
            population_name="protons",
            write_timestamps=timestamps,
        )

    **quantity**, one of:

        * "domain" (default), particles living within the interior of the simulation patches - the right default for most uses.
        * "levelGhost", particles living in the ghost region between two AMR levels - mainly useful for debugging refinement/coarsening.
        * "space_box", particles living within a user-given spatial region - requires the additional **extent** parameter (the box, in the same units as the domain).

    **population_name** (``str``), mandatory, name of the ion population to
    write particles for, as declared in
    :class:`~pyphare.pharein.MaxwellianFluidModel`.
    """

    particle_quantities = ["space_box", "domain", "levelGhost"]
    type = "particle"

    def __init__(self, **kwargs):
        super(ParticleDiagnostics, self).__init__(
            ParticleDiagnostics.type
            + str(global_vars.sim.count_diagnostics(ParticleDiagnostics.type)),
            **kwargs,
        )

    def _setSubTypeAttributes(self, **kwargs):
        # domain is good default for users who should not worry about what that means
        # even less about ghosts...
        kwargs["quantity"] = kwargs.get("quantity", "domain")

        if kwargs["quantity"] not in ParticleDiagnostics.particle_quantities:
            error_msg = "Error: '{}' not a valid particle diagnostics : " + ", ".join(
                ParticleDiagnostics.particle_quantities
            )
            raise ValueError(error_msg.format(kwargs["quantity"]))

        self.quantity = kwargs["quantity"]

        self.space_box(**kwargs)

        if "population_name" not in kwargs:
            raise ValueError("Error: missing population_name")
        else:
            self.population_name = kwargs["population_name"]

        if not population_in_model(self.population_name):
            raise ValueError(
                "Error: population '{}' not in simulation initial model".format(
                    self.population_name
                )
            )

        self.quantity = "/ions/pop/" + self.population_name + "/" + self.quantity

    def space_box(self, **kwargs):
        if "extent" not in kwargs and self.quantity == "space_box":
            raise ValueError(
                "Error: missing 'extent' parameter required by 'space_box' the ParticleDiagnostics type"
            )
        elif "extent" in kwargs:
            self.extent = kwargs["extent"]

    def to_dict(self):
        return super().to_dict(
            type(self).type,
            extent=", ".join([str(x) for x in self.extent]),
            population_name=self.population_name,
        )


# ------------------------------------------------------------------------------


class MetaDiagnostics(Diagnostics):
    """
    Writes bookkeeping information about the AMR hierarchy itself, rather
    than a physical quantity.

    **Usage example:**

    .. code-block:: python

        from pyphare.pharein import MetaDiagnostics

        MetaDiagnostics(quantity="tags", write_timestamps=timestamps)

    **quantity**, one of:

        * "tags" - the per-cell refinement tags produced by
          ``refinement="tagging"``, i.e. which cells were flagged for
          refinement at each write time. Useful to inspect/tune
          `tagging_threshold` and `tag_buffer`.
    """

    meta_quantities = ["tags"]
    type = "meta"

    def __init__(self, **kwargs):
        super(MetaDiagnostics, self).__init__(
            MetaDiagnostics.type
            + str(global_vars.sim.count_diagnostics(MetaDiagnostics.type)),
            **kwargs,
        )

    def _setSubTypeAttributes(self, **kwargs):
        if kwargs["quantity"] not in MetaDiagnostics.meta_quantities:
            error_msg = "Error: '{}' not a valid meta diagnostics : " + ", ".join(
                MetaDiagnostics.meta_quantities
            )
            raise ValueError(error_msg.format(kwargs["quantity"]))

        self.quantity = f"/{kwargs['quantity']}"

    def to_dict(self):
        return super().to_dict(type(self).type)


# ------------------------------------------------------------------------------


class InfoDiagnostics(Diagnostics):
    """
    Writes lightweight run-monitoring information to disk, cheap enough to
    dump at every time step.

    **Usage example:**

    .. code-block:: python

        from pyphare.pharein import InfoDiagnostics

        InfoDiagnostics(quantity="particle_count")

    **quantity**, one of:

        * "particle_count" - number of particles per patch/level, useful to
          monitor AMR/load-balancing behavior over the run.

    Unlike other diagnostics, **write_timestamps** defaults to every
    simulation time step if not given explicitly.
    """

    info_quantities = ["particle_count"]
    type = "info"

    @classmethod
    def default_kwargs(cls, **kwargs):
        if "write_timestamps" not in kwargs:
            kwargs["write_timestamps"] = all_timestamps(global_vars.sim)
        return kwargs

    def __init__(self, **kwargs):
        super(InfoDiagnostics, self).__init__(
            InfoDiagnostics.type
            + str(global_vars.sim.count_diagnostics(InfoDiagnostics.type)),
            **InfoDiagnostics.default_kwargs(**kwargs),
        )

    def _setSubTypeAttributes(self, **kwargs):
        if kwargs["quantity"] not in InfoDiagnostics.info_quantities:
            error_msg = "Error: '{}' not a valid info diagnostics : " + ", ".join(
                InfoDiagnostics.info_quantities
            )
            raise ValueError(error_msg.format(kwargs["quantity"]))

        self.quantity = f"/{kwargs['quantity']}"

    def to_dict(self):
        return super().to_dict(type(self).type)
