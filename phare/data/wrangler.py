class DataWrangler:
    is_primal = {"bx": True, "by": False, "bz": False, "ex": False, "ey": True, "ez": True}

    def __init__(self, sim, hier):
        import phare.pharein as ph, phare_lib.data_wrangler as data_wrangler

        self.dim = ph.globals.sim.dims
        self.interp = ph.globals.sim.interp_order
        self.cpp = getattr(
            data_wrangler,
            "DataWrangler_" + str(self.dim) + "_" + str(self.interp),
        )(sim, hier)

    def getPatchLevel(self, lvl):
        return self.cpp.getPatchLevel(lvl)

    def _lvl0FullContiguous(self, input, is_primal=True):
        return self.cpp.sync_merge(input, is_primal)

    def lvl0IonDensity(self):
        return self._lvl0FullContiguous(self.getPatchLevel(0).getDensity())

    def lvl0BulkVelocity(self):
        return {
            xyz: self._lvl0FullContiguous(bv)
            for xyz, bv in self.getPatchLevel(0).getBulkVelocity().items()
        }

    def lvl0PopDensity(self):
        return {
            pop: self._lvl0FullContiguous(density)
            for pop, density in self.getPatchLevel(0).getPopDensities().items()
        }

    def lvl0PopFluxs(self):
        return {
            pop: {xyz: self._lvl0FullContiguous(data) for xyz, data in flux.items()}
            for pop, flux in self.getPatchLevel(0).getPopFluxs().items()
        }

    def extract_is_primal_key_from(self, em_xyz):
        """ extract "ex" from "EM_E_x"  """
        return "".join(em_xyz.lower().split("_"))[2:]

    def lvl0EM(self):
        return {
            em: {
                em_xyz: self.cpp.sync_merge(
                    data, DataWrangler.is_primal[self.extract_is_primal_key_from(em_xyz)]
                )
                for em_xyz, data in xyz_map.items()
            }
            for em, xyz_map in self.getPatchLevel(0).getEM().items()
        }


# for pop, particles in dw.getPatchLevel(0).getParticles().items():
#     print("pop :", pop)
#     for key, patches in particles.items():
#         print("\tkey :", key)
#         for patch in patches:
#             print("\t\t", patch.patchID, "size:", patch.data.size())
