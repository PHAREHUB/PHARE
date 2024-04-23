import numpy as np

from ..core import phare_utilities
from . import global_vars


def all_timestamps(sim):
    nbr_dump_step = int(sim.final_time / sim.time_step) + 1
    return sim.time_step * np.arange(nbr_dump_step)


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

        accepted_keywords = [
            "path",
            "compute_timestamps",
            "population_name",
            "flush_every",
        ]
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
        raise RuntimeError(
            f"Error: timestamp({sim.time_step_nbr}) cannot be less than simulation.init_time({init_time}))"
        )
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


# ------------------------------------------------------------------------------


def try_cpp_dep_vers():
    try:
        from pyphare.cpp import cpp_etc_lib

        return cpp_etc_lib().phare_deps()
    except ImportError:
        return {}


class Diagnostics(object):
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
        self.compute_timestamps = validate_timestamps(
            self.__class__.__name__, "compute_timestamps", **kwargs
        )

        self.attributes = kwargs.get("attributes", {})
        self.attributes["git_hash"] = phare_utilities.top_git_hash()

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
                f"{self.__class__.__name__,}.flush_every cannot be negative"
            )

        self.__extent = None

        # if a diag already is registered we just contactenate the timestamps
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


# ------------------------------------------------------------------------------


class ElectromagDiagnostics(Diagnostics):
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
        return {
            "name": self.name,
            "type": ElectromagDiagnostics.type,
            "quantity": self.quantity,
            "write_timestamps": self.write_timestamps,
            "compute_timestamps": self.compute_timestamps,
            "path": self.path,
        }


# ------------------------------------------------------------------------------


def population_in_model(population):
    return population in [p for p in global_vars.sim.model.populations]


class FluidDiagnostics_(Diagnostics):
    fluid_quantities = [
        "density",
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
        return {
            "name": self.name,
            "type": FluidDiagnostics_.type,
            "quantity": self.quantity,
            "write_timestamps": self.write_timestamps,
            "compute_timestamps": self.compute_timestamps,
            "path": self.path,
            "population_name": self.population_name,
        }


def for_total_ions(**kwargs):
    return "population_name" not in kwargs


class FluidDiagnostics:
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
    particle_quantities = ["space_box", "domain", "levelGhost", "patchGhost"]
    type = "particle"

    def __init__(self, **kwargs):
        super(ParticleDiagnostics, self).__init__(
            ParticleDiagnostics.type
            + str(global_vars.sim.count_diagnostics(ParticleDiagnostics.type)),
            **kwargs,
        )

    def _setSubTypeAttributes(self, **kwargs):
        if kwargs["quantity"] not in ParticleDiagnostics.particle_quantities:
            error_msg = "Error: '{}' not a valid particle diagnostics : " + ", ".join(
                ParticleDiagnostics.particle_quantities
            )
            raise ValueError(error_msg.format(kwargs["quantity"]))

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
        return {
            "name": self.name,
            "type": ParticleDiagnostics.type,
            "quantity": self.quantity,
            "write_timestamps": self.write_timestamps,
            "compute_timestamps": self.compute_timestamps,
            "path": self.path,
            "extent": ", ".join([str(x) for x in self.extent]),
            "population_name": self.population_name,
        }


# ------------------------------------------------------------------------------


class MetaDiagnostics(Diagnostics):
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
        return {
            "name": self.name,
            "type": MetaDiagnostics.type,
            "quantity": self.quantity,
            "write_timestamps": self.write_timestamps,
            "compute_timestamps": self.compute_timestamps,
            "path": self.path,
        }


# ------------------------------------------------------------------------------


class InfoDiagnostics(Diagnostics):
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
        return {
            "name": self.name,
            "type": InfoDiagnostics.type,
            "quantity": self.quantity,
            "write_timestamps": self.write_timestamps,
            "path": self.path,
        }
