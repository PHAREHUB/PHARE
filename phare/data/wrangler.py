primal_datasets = {
    "bx": True,
    "by": False,
    "bz": False,
    "ex": False,
    "ey": True,
    "ez": True,
}


def is_primal_dataset(em_xyz):
    return primal_datasets[_extract_primal_key(em_xyz)]


def _extract_primal_key(em_xyz):
    """ extract "ex" from "EM_E_x"  """
    return "".join(em_xyz.lower().split("_"))[2:]


class DataWrangler:
    def __init__(self, sim_cpp, hier_cpp):
        import phare.pharein as ph, phare.data.data_wrangler

        self.dim = ph.globals.sim.dims
        self.interp = ph.globals.sim.interp_order
        self.cpp_class = "DataWrangler_" + str(self.dim) + "_" + str(self.interp)
        self.cpp = getattr(phare.data.data_wrangler, self.cpp_class)(sim_cpp, hier_cpp)

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

    def lvl0EM(self):
        return {
            em: {
                em_xyz: self.cpp.sync_merge(data, is_primal_dataset(em_xyz))
                for em_xyz, data in xyz_map.items()
            }
            for em, xyz_map in self.getPatchLevel(0).getEM().items()
        }


## we have no need for exposing particles currently, but this is the structure built from C++
# for pop_name, particles in dw.getPatchLevel(0).getParticles().items():
#     print("pop_name :", pop_name)
#     for key, patches in particles.items():
#         print("\tkey :", key)
#         for patch in patches:
#             print("\t\t", patch.patchID, "size:", patch.data.size())
