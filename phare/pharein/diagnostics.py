
from . import phare_utilities
#from .globals import objects
from . import globals

# ------------------------------------------------------------------------------

def check_last_iteration(**kwargs):
    if kwargs['last_iteration'] < kwargs['start_iteration']:
        raise ValueError("Error: 'last_iteration' should be equal or larger than 'start_iteration'")

# ------------------------------------------------------------------------------


def diagnostics_checker(func):
    def wrapper(diagnostics_object, name, **kwargs):

        accepted_keywords = ['write_timestamps',
                             'start_iteration', 'last_iteration',
                             'path', 'compute_timestamps']

        mandatory_keywords = ['write_timestamps', 'diag_type']

        # check that all passed keywords are in the accepted keyword list
        # wrong_kwds = phare_utilities.not_in_keywords_list(accepted_keywords, **kwargs)
        # if len(wrong_kwds) > 0:
        #     raise ValueError("Error: invalid arguments - " + " ".join(wrong_kwds))

        # check if some mandatory keywords are not missing
        missing_mandatory_kwds = phare_utilities.check_mandatory_keywords(mandatory_keywords, **kwargs)
        if len(missing_mandatory_kwds) > 0:
            raise ValueError("Error: missing mandatory parameters : " + ', '.join(missing_mandatory_kwds))

        try:
            # just take mandatory arguments from the dict
            # since if we arrived here we are sure they are there
            kwargs['start_iteration'] = kwargs.get("start_iteration", 0)
            kwargs['last_iteration'] = kwargs.get("last_iteration", 0)
            check_last_iteration(**kwargs)

            kwargs['path'] = kwargs.get("path", './')

            return func(diagnostics_object, name, **kwargs)

        except ValueError as msg:
            print(msg)

    return wrapper


# ------------------------------------------------------------------------------


class Diagnostics(object):

    @diagnostics_checker
    def __init__(self,name, **kwargs):

        self.name = name

        self.start_iteration = kwargs['start_iteration']
        self.last_iteration = kwargs['last_iteration']
        self.path = kwargs['path']

        self.write_timestamps = kwargs['write_timestamps'] #[0, 1, 2]
        self.compute_timestamps = kwargs['compute_timestamps']

        self.__extent = None

        if globals.sim is None:
            raise RuntimeError("A simulation must be created before adding diagnostics")

        globals.sim.add_diagnostics(self)

    def iteration_interval(self):
        return self.start_iteration, self.last_iteration

    def extent(self):
        return self.__extent

# ------------------------------------------------------------------------------


class ElectromagDiagnostics(Diagnostics):

    em_diag_types = ['E', 'B']
    category = "electromag"

    def __init__(self, **kwargs):
        super(ElectromagDiagnostics, self).__init__(ElectromagDiagnostics.category \
                                                    + str(globals.sim.count_diagnostics(ElectromagDiagnostics.category))
                                                    , **kwargs)


        if 'diag_type' not in kwargs:
            raise ValueError("Error: missing diag_type parameter")
        else:
            if kwargs['diag_type'] not in ElectromagDiagnostics.em_diag_types:
                error_msg = "Error: '{}' not a valid electromag diagnostics : " + ', '.join(ElectromagDiagnostics.em_diag_types)
                raise ValueError(error_msg.format(kwargs['diag_type']))
            else:
                self.diag_type = "/EM_" + kwargs['diag_type']

    def to_dict(self):
        return {"name": self.name,
                "diag_category": ElectromagDiagnostics.category,
                "diag_type": self.diag_type,
                "write_timestamps": self.write_timestamps,
                "compute_timestamps": self.compute_timestamps,
                "start_iteration": self.start_iteration,
                "last_iteration": self.last_iteration,
                "path": self.path}

# ------------------------------------------------------------------------------

def population_in_model(population):
    return population in [p for p in globals.sim.model.populations] + ["all","ions"]




class FluidDiagnostics (Diagnostics):

    fluid_diag_types = ['density', 'flux', 'bulkVelocity']
    category = "fluid"

    def __init__(self, **kwargs):
        super(FluidDiagnostics, self).__init__(FluidDiagnostics.category \
                                               + str(globals.sim.count_diagnostics(FluidDiagnostics.category)),
                                               **kwargs)
        if 'diag_type' not in kwargs:
            raise ValueError("Error: missing diag_type parameter")

        self.population_name = None
        if 'population_name' not in kwargs and kwargs['diag_type'] == "flux":
            raise ValueError("Error: missing population_name")
        elif 'population_name' in kwargs:
            self.population_name = kwargs['population_name']

        if kwargs['diag_type'] not in FluidDiagnostics.fluid_diag_types:
            error_msg = "Error: '{}' not a valid fluid diagnostics : " + ', '.join(FluidDiagnostics.fluid_diag_types)
            raise ValueError(error_msg.format(kwargs['diag_type']))
        elif kwargs['diag_type'] == 'flux' and kwargs['population_name'] == "ions":
            raise ValueError("'flux' is only available for specific populations, try 'bulkVelocity")
        else:
            self.diag_type = kwargs['diag_type']

        if self.population_name is None:
            self.diag_type = "/ions/" + self.diag_type
        else:
            if not population_in_model(self.population_name):
                raise ValueError("Error: population '{}' not in simulation initial model".format(self.population_name))
            self.diag_type = "/ions/pop/ions_" + self.population_name + "/" + self.diag_type


    def to_dict(self):
        return {"name": self.name,
                "diag_category": FluidDiagnostics.category,
                "diag_type": self.diag_type,
                "write_timestamps": self.write_timestamps,
                "compute_timestamps": self.compute_timestamps,
                "start_iteration": self.start_iteration,
                "last_iteration": self.last_iteration,
                "path": self.path,
                "population_name": self.population_name}



# ------------------------------------------------------------------------------


class ParticleDiagnostics(Diagnostics):

    particle_diag_types = ['space_box', 'domain', 'levelGhost', 'patchGhost']
    category = "particle"

    def __init__(self, **kwargs):
        super(ParticleDiagnostics, self).__init__(ParticleDiagnostics.category \
                                                  + str(globals.sim.count_diagnostics(ParticleDiagnostics.category)),
                                                  **kwargs)

        if 'diag_type' not in kwargs:
            raise ValueError("Error: missing diag_type parameter")

        if kwargs['diag_type'] not in ParticleDiagnostics.particle_diag_types:
            error_msg = "Error: '{}' not a valid particle diagnostics : " + ', '.join(ParticleDiagnostics.particle_diag_types)
            raise ValueError(error_msg.format(kwargs['diag_type']))

        self.diag_type = kwargs['diag_type']

        self.space_box(**kwargs)

        if 'population_name' not in kwargs:
            raise ValueError("Error: missing population_name")
        else:
            self.population_name = kwargs['population_name']

        if not population_in_model(self.population_name):
            raise ValueError("Error: population '{}' not in simulation initial model".format(self.population_name))

        self.diag_type = "/ions/pop/ions_" + self.population_name + "/" + self.diag_type

    def space_box(self, **kwargs):

        if 'extent' not in kwargs and self.diag_type == 'space_box':
            raise ValueError("Error: missing 'extent' parameter required by 'space_box' the ParticleDiagnostics type")
        elif 'extent' in kwargs:
            self.extent = kwargs['extent']

    def to_dict(self):
        return {"name": self.name,
                "diag_category": ParticleDiagnostics.category,
                "diag_type": self.diag_type,
                "write_timestamps": self.write_timestamps,
                "compute_timestamps": self.compute_timestamps,
                "start_iteration": self.start_iteration,
                "last_iteration": self.last_iteration,
                "path": self.path,
                "extent": ", ".join([str(x) for x in self.extent]),
                "population_name":self.population_name}
