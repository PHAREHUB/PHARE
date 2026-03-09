from pyphare.core import gridlayout


class DataWrangler:
    def __init__(self, simulator):
        from pyphare import cpp

        sim = simulator.simulation
        self.dim = sim.ndim
        self.interp = sim.interp_order
        self.refined_particle_nbr = sim.refined_particle_nbr
        self.cpp = getattr(cpp.cpp_lib(sim), "DataWrangler")(
            simulator.cpp_sim, simulator.cpp_hier
        )
        # todo mhd
        self.modelsPerLevel = ["Hybrid" for lvl in range(self.cpp.getNumberOfLevels())]

    def kill(self):
        del self.cpp

    def getNumberOfLevels(self):
        return self.cpp.getNumberOfLevels()

    def getMHDPatchLevel(self, lvl):
        return self.cpp.getMHDPatchLevel(lvl)

    def getHybridPatchLevel(self, lvl):
        return self.cpp.getHybridPatchLevel(lvl)

    def getPatchLevel(self, lvl):
        lvlModel = self.modelsPerLevel[lvl]
        assert lvlModel in ["Hybrid", "MHD"]
        if lvlModel == "MHD":
            return self.getMHDPatchLevel(lvl)
        return self.getHybridPatchLevel(lvl)

    def sync(self, patch_datas):
        return self.cpp.sync(patch_datas)
