
from ..core import phare_utilities
from . import global_vars

# ------------------------------------------------------------------------------



def diagnostics_checker(func):
    def wrapper(diagnostics_object, name, **kwargs):

        mandatory_keywords = ['write_timestamps', 'quantity']

        # check if some mandatory keywords are not missing
        missing_mandatory_kwds = phare_utilities.check_mandatory_keywords(mandatory_keywords, **kwargs)
        if len(missing_mandatory_kwds) > 0:
            raise RuntimeError("Error: missing mandatory parameters : " + ', '.join(missing_mandatory_kwds))

        accepted_keywords = ['path', 'compute_timestamps', 'population_name', 'flush_every']
        accepted_keywords += mandatory_keywords

        # check that all passed keywords are in the accepted keyword list
        wrong_kwds = phare_utilities.not_in_keywords_list(accepted_keywords, **kwargs)
        if len(wrong_kwds) > 0:
            raise RuntimeError("Error: invalid arguments - " + " ".join(wrong_kwds))

        try:
            # just take mandatory arguments from the dict
            # since if we arrived here we are sure they are there

            kwargs['path'] = kwargs.get("path", './')

            return func(diagnostics_object, name, **kwargs)

        except ValueError as msg:
            print(msg)

    return wrapper


import numpy as np
# ------------------------------------------------------------------------------
def validate_timestamps(clazz, **kwargs):
    sim = global_vars.sim

    for key in ["write_timestamps", "compute_timestamps"]:
        timestamps = kwargs[key]

        if np.any(timestamps < sim.init_time):
            raise RuntimeError(f"Error: timestamp({sim.time_step_nbr}) cannot be less than simulation.init_time({sim.init_time}))")
        if np.any(timestamps > sim.final_time):
            raise RuntimeError(f"Error: timestamp({sim.time_step_nbr}) cannot be greater than simulation.final_time({sim.final_time}))")
        if not np.all(np.diff(timestamps) >= 0):
            raise RuntimeError(f"Error: {clazz}.{key} not in ascending order)")
        if not np.all(np.abs(timestamps / sim.time_step - np.rint(timestamps/sim.time_step) < 1e-9)):
            raise RuntimeError(f"Error: {clazz}.{key} is inconsistent with simulation.time_step")





# ------------------------------------------------------------------------------

def try_cpp_dep_vers():
    try:
        from pyphare.cpp import cpp_lib
        return cpp_lib().phare_deps()
    except ImportError:
        return {}

class Diagnostics(object):

    h5_flush_never = 0
    cpp_dep_vers = try_cpp_dep_vers()

    @diagnostics_checker
    def __init__(self,name, **kwargs):

        if global_vars.sim is None:
            raise RuntimeError("A simulation must be created before adding diagnostics")

        self.name = name
        self.path = kwargs['path']

        validate_timestamps(self.__class__.__name__, **kwargs)
        self.write_timestamps = kwargs['write_timestamps'] #[0, 1, 2]
        self.compute_timestamps = kwargs['compute_timestamps']

        self.attributes = kwargs.get("attributes", {})
        self.attributes["git_hash"]  = phare_utilities.top_git_hash()

        for dep, dep_ver in Diagnostics.cpp_dep_vers.items():
            self.attributes[f"{dep}_version"] = dep_ver

        for key in self.attributes:
            self.attributes[key] = self.attributes[key]

        self._setSubTypeAttributes(**kwargs)
        self.flush_every = kwargs.get("flush_every", 1) # flushes every dump, safe, but costly

        if self.flush_every < 0:
            raise RuntimeError(f"{self.__class__.__name__,}.flush_every cannot be negative")

        if any([self.quantity == diagnostic.quantity for diagnostic in global_vars.sim.diagnostics]):
            raise RuntimeError(f"Error: Diagnostic ({kwargs['quantity']}) already registered")

        self.__extent = None
        global_vars.sim.add_diagnostics(self)


    def extent(self):
        return self.__extent

# ------------------------------------------------------------------------------


class ElectromagDiagnostics(Diagnostics):

    em_quantities = ['E', 'B']
    type = "electromag"

    def __init__(self, **kwargs):
        super(ElectromagDiagnostics, self).__init__(ElectromagDiagnostics.type \
                                                    + str(global_vars.sim.count_diagnostics(ElectromagDiagnostics.type))
                                                    , **kwargs)

    def _setSubTypeAttributes(self, **kwargs):

        if kwargs['quantity'] not in ElectromagDiagnostics.em_quantities:
            error_msg = "Error: '{}' not a valid electromag diagnostics : " + ', '.join(ElectromagDiagnostics.em_quantities)
            raise ValueError(error_msg.format(kwargs['quantity']))
        else:
            self.quantity = "/EM_" + kwargs['quantity']

    def to_dict(self):
        return {"name": self.name,
                "type": ElectromagDiagnostics.type,
                "quantity": self.quantity,
                "write_timestamps": self.write_timestamps,
                "compute_timestamps": self.compute_timestamps,
                "path": self.path}

# ------------------------------------------------------------------------------

def population_in_model(population):
    return population in [p for p in global_vars.sim.model.populations] + ["all","ions"]




class FluidDiagnostics (Diagnostics):

    fluid_quantities = ['density', 'flux', 'bulkVelocity']
    type = "fluid"

    def __init__(self, **kwargs):
        super(FluidDiagnostics, self).__init__(FluidDiagnostics.type \
                                               + str(global_vars.sim.count_diagnostics(FluidDiagnostics.type)),
                                               **kwargs)

    def _setSubTypeAttributes(self, **kwargs):

        self.population_name = None
        if 'population_name' not in kwargs and kwargs['quantity'] == "flux":
            raise ValueError("Error: missing population_name")
        elif 'population_name' in kwargs:
            self.population_name = kwargs['population_name']

        if kwargs['quantity'] not in FluidDiagnostics.fluid_quantities:
            error_msg = "Error: '{}' not a valid fluid diagnostics : " + ', '.join(FluidDiagnostics.fluid_quantities)
            raise ValueError(error_msg.format(kwargs['quantity']))
        elif kwargs['quantity'] == 'flux' and kwargs['population_name'] == "ions":
            raise ValueError("'flux' is only available for specific populations, try 'bulkVelocity")
        else:
            self.quantity = kwargs['quantity']

        if self.population_name is None:
            self.quantity = "/ions/" + self.quantity
        else:
            if not population_in_model(self.population_name):
                raise ValueError("Error: population '{}' not in simulation initial model".format(self.population_name))
            self.quantity = "/ions/pop/" + self.population_name + "/" + self.quantity


    def to_dict(self):
        return {"name": self.name,
                "type": FluidDiagnostics.type,
                "quantity": self.quantity,
                "write_timestamps": self.write_timestamps,
                "compute_timestamps": self.compute_timestamps,
                "path": self.path,
                "population_name": self.population_name}



# ------------------------------------------------------------------------------


class ParticleDiagnostics(Diagnostics):

    particle_quantities = ['space_box', 'domain', 'levelGhost', 'patchGhost']
    type = "particle"

    def __init__(self, **kwargs):
        super(ParticleDiagnostics, self).__init__(ParticleDiagnostics.type \
                                                  + str(global_vars.sim.count_diagnostics(ParticleDiagnostics.type)),
                                                  **kwargs)

    def _setSubTypeAttributes(self, **kwargs):

        if kwargs['quantity'] not in ParticleDiagnostics.particle_quantities:
            error_msg = "Error: '{}' not a valid particle diagnostics : " + ', '.join(ParticleDiagnostics.particle_quantities)
            raise ValueError(error_msg.format(kwargs['quantity']))

        self.quantity = kwargs['quantity']

        self.space_box(**kwargs)

        if 'population_name' not in kwargs:
            raise ValueError("Error: missing population_name")
        else:
            self.population_name = kwargs['population_name']

        if not population_in_model(self.population_name):
            raise ValueError("Error: population '{}' not in simulation initial model".format(self.population_name))

        self.quantity = "/ions/pop/" + self.population_name + "/" + self.quantity

    def space_box(self, **kwargs):

        if 'extent' not in kwargs and self.quantity == 'space_box':
            raise ValueError("Error: missing 'extent' parameter required by 'space_box' the ParticleDiagnostics type")
        elif 'extent' in kwargs:
            self.extent = kwargs['extent']

    def to_dict(self):
        return {"name": self.name,
                "type": ParticleDiagnostics.type,
                "quantity": self.quantity,
                "write_timestamps": self.write_timestamps,
                "compute_timestamps": self.compute_timestamps,
                "path": self.path,
                "extent": ", ".join([str(x) for x in self.extent]),
                "population_name":self.population_name}
