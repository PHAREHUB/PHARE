

from ..core import phare_utilities
from ..core.gridlayout import yee_element_is_primal
from . import global_vars


class MaxwellianFluidModel(object):

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
            return lambda x:value + x*0
        if self.dim == 2:
            return lambda x,y:value
        if self.dim == 3:
            return lambda x,y,z:value


    def __init__(self, bx = None,
                       by = None,
                       bz = None,
                       **kwargs):

        if global_vars.sim is None:
            raise RuntimeError("A simulation must be declared before a model")

        if global_vars.sim.model is not None:
            raise RuntimeError("A model is already created")

        self.dim = global_vars.sim.ndim
        bx = self.defaulter(bx, 1.)
        by = self.defaulter(by, 0.)
        bz = self.defaulter(bz, 0.)


        self.model_dict = {"model": "model", "model_name": "custom"}



        self.model_dict.update({"bx": bx,
                                "by": by,
                                "bz": bz})

        self.populations = list(kwargs.keys())
        for population in self.populations:
            self.add_population(population, **kwargs[population])

        self.validate(global_vars.sim)

        global_vars.sim.set_model(self)


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
        import numpy as np
        from ..core.gridlayout import GridLayout, directions as all_directions
        from ..core.box import Box

        directions = all_directions[:sim.ndim]
        domain_box = Box([0] * sim.ndim, np.asarray(sim.cells) - 1)
        assert len(sim.origin) == domain_box.ndim
        ndim = domain_box.ndim

        def _primal_directions(qty):
            return [1 if yee_element_is_primal(qty, direction) else 0 for direction in directions]

        def _nbrGhosts(layout, qty, primal_directions):
            return [
                layout.nbrGhosts(sim.interp_order, "primal" if primal_directions[dim] else "dual") for dim in range(ndim)
            ]

        def split_to_columns(ar):
            ar = ar.reshape(-1, ndim) # drop outer arrays if any
            return [ar[:,dim] for dim in range(ndim)]

        def _nd_check(layout, qty, nbrGhosts, fn):
            from pyphare.pharesee.geometry import hierarchy_overlaps
            from pyphare.pharesee.hierarchy import hierarchy_from_box
            from ..core import box as boxm

            hier = hierarchy_from_box(domain_box, nbrGhosts)
            coords = layout.meshCoords(qty)

            for ilvl, overlaps in hierarchy_overlaps(hier).items():
                for overlap in overlaps:
                    (pd1, pd2), box, offsets = overlap["pdatas"], overlap["box"], overlap["offset"]
                    slice1 = boxm.select(coords, boxm.amr_to_local(box, boxm.shift(pd1.ghost_box, offsets[0])))
                    slice2 = boxm.select(coords, boxm.amr_to_local(box, boxm.shift(pd2.ghost_box, offsets[1])))
                    if not np.allclose(fn(*split_to_columns(slice1)) , fn(*split_to_columns(slice2)), atol=1e-15):
                        return False
            return True

        def periodic_function_check(vec_field, dic):
            valid = True
            for xyz in ["x", "y", "z"]:
                qty = vec_field + xyz
                layout = GridLayout(domain_box, sim.origin, sim.dl, sim.interp_order)
                nbrGhosts = _nbrGhosts(layout, qty, _primal_directions(qty))
                valid &= _nd_check(layout, qty, nbrGhosts, fn=dic[qty])
            return valid

        if sim.boundary_types[0] == "periodic":
            model_dict = self.model_dict
            valid = True
            for pop in self.populations:
                for v in ["vth", "v"]:
                    valid &= periodic_function_check(v, model_dict[pop])
            valid &= periodic_function_check("b", model_dict)
            if not valid:
                print("Warning: Simulation is periodic but some functions are not")
                if sim.strict:
                    raise RuntimeError("Simulation is not periodic")
