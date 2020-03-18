class DataWrangler:
    is_primal = {"bx": 1, "by": 0, "bz": 0, "ex": 0, "ey": 1, "ez": 1}

    def __init__(self, sim, hier):
        import phare.pharein as ph, phare.data.data_wrangler

        self.dim = len(ph.globals.sim.dl)
        self.interp = ph.globals.sim.interp_order
        self.cpp = getattr(
            phare.data.data_wrangler,
            "DataWrangler_" + str(self.dim) + "_" + str(self.interp),
        )(sim, hier)

    def getPatchLevel(self, lvl):
        return self.cpp.getPatchLevel(lvl)

    def _lvl0FullContigous(self, input, is_primal=True):
        return self.cpp.sync_merge(input, is_primal)

    def lvl0IonDensity(self):
        return self._lvl0FullContigous(self.getPatchLevel(0).getDensity())

    def lvl0BulkVelocity(self):
        return {
            xyz: self._lvl0FullContigous(bv)
            for xyz, bv in self.getPatchLevel(0).getBulkVelocity().items()
        }

    def lvl0PopDensity(self):
        return {
            pop: self._lvl0FullContigous(density)
            for pop, density in self.getPatchLevel(0).getPopDensities().items()
        }

    def lvl0PopFluxs(self):
        return {
            pop: {xyz: self._lvl0FullContigous(data) for xyz, data in flux.items()}
            for pop, flux in self.getPatchLevel(0).getPopFluxs().items()
        }

    def lvl0EM(self):
        return {
            em: {
                xyz: self.cpp.sync_merge(
                    data, DataWrangler.is_primal["".join(xyz.lower().split("_"))[2:]]
                )
                for xyz, data in xyz_map.items()
            }
            for em, xyz_map in self.getPatchLevel(0).getEM().items()
        }


# for pop, particles in dw.getPatchLevel(0).getParticles().items():
#     print("pop :", pop)
#     for key, patches in particles.items():
#         print("\tkey :", key)
#         for patch in patches:
#             print("\t\t", patch.patchID, "size:", patch.data.size())
