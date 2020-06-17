
from ..core import phare_utilities
from . import globals

class MaxwellianFluidModel(object):

    def defaulter(self, input, value):
        if input is not None:
            import inspect
            argspec = inspect.getargspec(input)
            has_args_but_not_kwargs = argspec.varargs != None and argspec.keywords == None
            sig = inspect.signature(input)
            params = sig.parameters
            param_per_dim = len(params) == self.dim
            param_is_varargs = len(params) == 1 and has_args_but_not_kwargs
            assert param_per_dim or param_is_varargs
            return input
        if self.dim == 1:
            return lambda x:value + x*0
        if self.dim == 2:
            return lambda x,y:value
        if self.dim == 3:
            return lambda x,y,z:value


    def __init__(self, bx = None,
                       by = None,
                       bz = None,
                       **kwargs):

        if globals.sim is None:
            raise RuntimeError("A simulation must be declared before a model")

        if globals.sim.model is not None:
            raise RuntimeError("A model is already created")

        self.dim = globals.sim.dims
        bx = self.defaulter(bx, 1.)
        by = self.defaulter(by, 0.)
        bz = self.defaulter(bz, 0.)


        self.model_dict = {"model": "model", "model_name": "custom"}



        self.model_dict.update({"bx": bx,
                                "by": by,
                                "bz": bz})

        self.populations = kwargs.keys()
        for population in self.populations:
            self.add_population(population, **kwargs[population])

        self.validate(globals.sim)

        globals.sim.set_model(self)


# ------------------------------------------------------------------------------

    def nbr_populations(self):
        """
        returns the number of species currently registered in the model
        """
        return len(self.populations)

#------------------------------------------------------------------------------

    def add_population(self, name,
                        charge=1.,
                        mass=1.,
                        nbr_part_per_cell=100,
                        density= None,
                        vbulkx = None,
                        vbulky = None,
                        vbulkz = None,
                        vthx   = None,
                        vthy   = None,
                        vthz   = None,
                        init   = {}):
        """
        add a particle population to the current model

        add_population(name,charge=1, mass=1, nbrPartCell=100, density=1, vbulk=(0,0,0), beta=1, anisotropy=1)

        Parameters:
        -----------
        name        : name of the species, str

        Optional Parameters:
        -------------------
        charge      : charge of the species particles, float (default = 1.)
        nbrPartCell : number of particles per cell, int (default = 100)
        density     : particle density, float (default = 1.)
        vbulk       : bulk velocity, tuple of size 3  (default = (0,0,0))
        beta        : beta of the species, float (default = 1)
        anisotropy  : Pperp/Ppara of the species, float (default = 1)
        """

        init_keys = ['seed']
        wrong_keys = phare_utilities.not_in_keywords_list(init_keys, **init)
        if len(wrong_keys) > 0:
            raise ValueError("Model Error: invalid init arguments - " + " ".join(wrong_keys))
        init["seed"] = init["seed"] if "seed" in init else None

        density = self.defaulter(density, 1.)

        vbulkx = self.defaulter(vbulkx, 0.)
        vbulky = self.defaulter(vbulky, 0.)
        vbulkz = self.defaulter(vbulkz, 0.)

        vthx   = self.defaulter(vthx, 1.)
        vthy   = self.defaulter(vthy, 1.)
        vthz   = self.defaulter(vthz, 1.)

        idx = str(self.nbr_populations())

        new_population = {name: {
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
                          "init": init}}

        keys = self.model_dict.keys()
        if name in keys:
            raise ValueError("population already registered")

        self.model_dict.update(new_population)
#------------------------------------------------------------------------------

    def to_dict(self):
        self.model_dict['nbr_ion_populations'] = self.nbr_populations()
        return self.model_dict

#------------------------------------------------------------------------------

    def validate(self, sim):
        import math

        def periodic_function_check(vec_field, dic):
            for xyz in ["x", "y", "z"]:
                fn = dic[vec_field + xyz]
                if not math.isclose(fn(*sim.origin), fn(*sim.simulation_domain()), rel_tol=1e-5):
                    return False
            return True

        if sim.boundary_types[0] == "periodic":
            model_dict = self.model_dict
            valid = True
            for pop_index, pop in enumerate(self.populations):
                for v in ["vth", "v"]:
                    valid &= periodic_function_check(v, model_dict[pop])
            for em in ["b"]:
                valid &= periodic_function_check(em, model_dict)
            if not valid:
                print("Warning: Simulation is periodic but some functions are not")
