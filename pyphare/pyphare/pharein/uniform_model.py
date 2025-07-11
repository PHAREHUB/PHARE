from . import global_vars


class UniformModel(object):
    def __init__(self, b=(1.0, 0.0, 0.0), e=(0.0, 0.0, 0.0), **kwargs):
        self.model_dict = {"model": "model", "model_name": "uniform"}

        if global_vars.sim is None:
            raise RuntimeError("A simulation must be declared before a model")

        if global_vars.sim.model is not None:
            raise RuntimeError("A model is already created")

        global_vars.sim.set_uniform_model(self)

        if len(b) != 3 or (not isinstance(b, tuple) and not isinstance(b, list)):
            raise ValueError("invalid B")
        if len(e) != 3 or (not isinstance(e, tuple) and not isinstance(e, list)):
            raise ValueError("invalid E")

        self.model_dict.update(
            {
                "bx": lambda x: b[0],
                "by": lambda x: b[1],
                "bz": lambda x: b[2],
                "ex": lambda x: e[0],
                "ey": lambda x: e[1],
                "ez": lambda x: e[2],
            }
        )

        self.populations = kwargs.keys()
        for population in self.populations:
            self.add_population(population, **kwargs[population])

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
        density=1.0,
        vbulk=(0.0, 0.0, 0.0),
        beta=1.0,
        anisotropy=1.0,
    ):
        """
        add a particle population to the current model

        add_population(name,charge=1, mass=1, nbrPartCell=100, density=1, vbulk=(0,0,0), beta=1, anisotropy=1)

        Parameters
        ----------
        name        : name of the species, str

        Other Parameters
        ----------------
        charge      : charge of the species particles, float (default = 1.)
        nbrPartCell : number of particles per cell, int (default = 100)
        density     : particle density, float (default = 1.)
        vbulk       : bulk velocity, tuple of size 3  (default = (0,0,0))
        beta        : beta of the species, float (default = 1)
        anisotropy  : Pperp/Ppara of the species, float (default = 1)
        """
        new_population = {
            name: {
                "charge": charge,
                "mass": mass,
                "density": lambda x: density,
                "vx": lambda x: vbulk[0],
                "vy": lambda x: vbulk[1],
                "vz": lambda x: vbulk[2],
                "beta": lambda x: beta,
                "anisotropy": lambda x: anisotropy,
                "nbrParticlesPerCell": nbr_part_per_cell,
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
