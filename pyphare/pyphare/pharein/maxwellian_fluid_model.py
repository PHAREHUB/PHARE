

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

        self.populations = kwargs.keys()
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
        import itertools
        import numpy as np
        from ..core.gridlayout import GridLayout, directions as all_directions
        from ..core.box import Box
        from pyphare.pharein import py_fn_wrapper

        dl = np.array(sim.dl)
        directions = all_directions[:sim.ndim]

        domain_box = Box([0] * sim.ndim, np.asarray(sim.cells) - 1)
        assert len(sim.origin) == domain_box.ndim

        def _primal_directions(qty):
            return [1 if yee_element_is_primal(qty, direction) else 0 for direction in directions]

        def _nbrGhosts(qty, layout, primal_directions):
            return [
                layout.nbrGhosts(sim.interp_order, "primal" if primal_directions[dim] else "dual") for dim in range(sim.ndim)
            ]

        def qty_shifts(qty, nbrGhosts):
          primal_shifts = np.arange(-5, 6)
          dual_shifts   = np.arange(-5, 5)
          shifts = []
          for dim in range(sim.ndim):
              shifts += [primal_shifts * dl[dim] if yee_element_is_primal(qty, directions[dim]) else dual_shifts *  dl[dim] + (.5 * dl[dim])]
          return shifts

        def _1d_check(layout, nbrGhosts, primal_directions, qty, fn):
            fnX = py_fn_wrapper(fn)(layout.yeeCoordsFor(qty, "x"))
            include_border = primal_directions[0]
            np.testing.assert_allclose(fnX[:nbrGhosts * 2 + include_border] , fnX[-(nbrGhosts * 2) - include_border:], atol=1e-15)
            return (
                np.allclose(fnX[:nbrGhosts * 2 + include_border] , fnX[-(nbrGhosts * 2) - include_border:], atol=1e-15)
            )


        def split_to_columns(ar):
            ar = ar.reshape(-1, sim.ndim) # drop outer arrays if any
            return [ar[:,dim] for dim in range(sim.ndim)]


        def _2d_check(nbrGhosts, primal_directions, qty, fn, shift):
            layout = GridLayout(domain_box, sim.origin + (shift * sim.dl), sim.dl, sim.interp_order)
            x = layout.yeeCoordsFor(qty, "x")[nbrGhosts[0]:-(nbrGhosts[0]-1)-(primal_directions[0])]
            y = layout.yeeCoordsFor(qty, "y")[nbrGhosts[1]:-(nbrGhosts[1]-1)-(primal_directions[1])]
            X,Y = np.meshgrid(x,y,indexing="ij")
            coords = np.array([X.flatten(),Y.flatten()]).T
            top = coords[:len(x)]
            bot = coords[-len(x):]
            left  = coords[0 :: len(x)]
            right = coords[len(x) - 1 :: len(x)]
            return (
                np.allclose(fn(*split_to_columns(left)) , fn(*split_to_columns(right)), atol=1e-15) and
                np.allclose(fn(*split_to_columns(top))  , fn(*split_to_columns(bot)), atol=1e-15)
            )

        def _3d_check(nbrGhosts, primal_directions, qty, fn, shift):
            def get_3d_faces(M):
                def get_face(M, dim, front_side):
                    return M[tuple((0 if front_side else -1) if i == dim else slice(None) for i in range(M.ndim))]
                return [get_face(M, dim, front_side) for dim in range(sim.ndim) for front_side in (True, False)]

            layout = GridLayout(domain_box, sim.origin + (shift * sim.dl), sim.dl, sim.interp_order)
            x = layout.yeeCoordsFor(qty, "x")[nbrGhosts[0]:-(nbrGhosts[0]-1)-(primal_directions[0])]
            y = layout.yeeCoordsFor(qty, "y")[nbrGhosts[1]:-(nbrGhosts[1]-1)-(primal_directions[1])]
            z = layout.yeeCoordsFor(qty, "z")[nbrGhosts[2]:-(nbrGhosts[2]-1)-(primal_directions[2])]
            coords = np.array([X.flatten(), Y.flatten(), Z.flatten()]).T.reshape((len(x), len(y), len(z), sim.ndim))
            left, right, top, bot, front, back = get_3d_faces(coords)
            return (
                np.allclose(fn(split_to_columns(left)) , fn(split_to_columns(right)), atol=1e-15) and
                np.allclose(fn(split_to_columns(top))  , fn(split_to_columns(bot)), atol=1e-15) and
                np.allclose(fn(split_to_columns(front)), fn(split_to_columns(back)), atol=1e-15)
            )

        def periodic_function_check(vec_field, dic):
            valid = True
            for xyz in ["x", "y", "z"]:
                qty = vec_field + xyz
                fn = dic[qty]

                layout = GridLayout(domain_box, sim.origin, sim.dl, sim.interp_order)
                primal_directions = _primal_directions(qty)
                nbrGhosts = _nbrGhosts(qty, layout, primal_directions)

                if sim.ndim == 1:
                    valid &= _1d_check(layout, nbrGhosts[0], primal_directions, qty, fn)

                if sim.ndim == 2:
                    permutations = itertools.product(*qty_shifts(qty, nbrGhosts))
                    valid &= all([_2d_check(nbrGhosts, primal_directions, qty, fn, np.asarray(perm)) for perm in permutations])

                if sim.ndim == 3:
                    # untested block - add in during 3d/splitting PR https://github.com/PHAREHUB/PHARE/pull/314
                    permutations = itertools.product(*qty_shifts(qty, nbrGhosts))
                    valid &= all([_3d_check(nbrGhosts, primal_directions, qty, fn, np.asarray(perm)) for perm in permutations])

            return valid

        if sim.boundary_types[0] == "periodic":
            model_dict = self.model_dict
            valid = True
            for pop_index, pop in enumerate(self.populations):
                for v in ["vth", "v"]:
                    valid &= periodic_function_check(v, model_dict[pop])
            valid &= periodic_function_check("b", model_dict)
            if not valid:
                print("Warning: Simulation is periodic but some functions are not")
                if sim.strict:
                    raise RuntimeError("Simulation is not periodic")
