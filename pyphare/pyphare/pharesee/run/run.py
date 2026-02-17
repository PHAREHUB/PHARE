#
#
#

import os
import glob
from pathlib import Path

from pyphare.pharesee.hierarchy import all_times_from
from pyphare.pharesee.hierarchy import default_time_from
from pyphare.pharesee.hierarchy import hierarchy_from
from pyphare.pharesee.hierarchy import hierarchy_compute as hc
from pyphare.pharesee.hierarchy import ScalarField, VectorField

from pyphare.core import phare_utilities as phut
from pyphare.pharesee.hierarchy.hierarchy_utils import compute_hier_from
from pyphare.pharesee.hierarchy.hierarchy_utils import flat_finest_field

from pyphare.logger import getLogger

from .man import RunMan

from .utils import (
    _compute_to_primal,
    _compute_pop_pressure,
    _compute_pressure,
    _compute_current,
    _compute_divB,
    _get_rank,
    make_interpolator,
)

logger = getLogger(__name__)


class Run:
    def __init__(self, path, default_time=None):
        self.path = path
        self.default_time_ = default_time
        self.available_diags = self._available_diags()

    def GetTags(self, time, merged=False, **kwargs):
        hier = self._get_hierarchy(time, "tags.h5")
        return self._get(hier, time, merged, "nearest")

    def GetB(self, time, merged=False, interp="nearest", all_primal=True, **kwargs):
        if merged:
            all_primal = False
        hier = self._get_hier_for(time, "EM_B", **kwargs)
        if not all_primal:
            return self._get(hier, time, merged, interp)

        return VectorField(compute_hier_from(_compute_to_primal, hier))

    def GetE(self, time, merged=False, interp="nearest", all_primal=True, **kwargs):
        if merged:
            all_primal = False
        hier = self._get_hier_for(time, "EM_E", **kwargs)
        if not all_primal:
            return self._get(hier, time, merged, interp)

        return VectorField(compute_hier_from(_compute_to_primal, hier))

    def GetMassDensity(self, time, merged=False, interp="nearest", **kwargs):
        hier = self._get_hier_for(time, "ions_mass_density", **kwargs)
        return ScalarField(self._get(hier, time, merged, interp))

    def GetNi(self, time, merged=False, interp="nearest", **kwargs):
        hier = self._get_hier_for(time, "ions_charge_density", **kwargs)
        return ScalarField(self._get(hier, time, merged, interp, drop_ghosts=True))

    def GetN(self, time, pop_name, merged=False, interp="nearest", **kwargs):
        hier = self._get_hier_for(time, f"ions_pop_{pop_name}_density", **kwargs)
        return ScalarField(self._get(hier, time, merged, interp))

    def GetVi(self, time, merged=False, interp="nearest", **kwargs):
        hier = self._get_hier_for(time, "ions_bulkVelocity", **kwargs)
        return VectorField(self._get(hier, time, merged, interp, drop_ghosts=True))
        return RunMan(self).GetVi(time, merged, interp, **kwargs)

    def GetFlux(self, time, pop_name, merged=False, interp="nearest", **kwargs):
        hier = self._get_hier_for(time, f"ions_pop_{pop_name}_flux", **kwargs)
        return VectorField(self._get(hier, time, merged, interp))

    def GetPressure(self, time, pop_name, merged=False, interp="nearest", **kwargs):
        # return RunMan(self).GetPressure(time, pop_name, merged, interp, **kwargs)
        M = self._get_hierarchy(
            time, f"ions_pop_{pop_name}_momentum_tensor.h5", **kwargs
        )
        V = self.GetFlux(time, pop_name, **kwargs)
        N = self.GetN(time, pop_name, **kwargs)
        P = compute_hier_from(
            _compute_pop_pressure,
            (M, V, N),
            popname=pop_name,
            mass=self.GetMass(pop_name, **kwargs),
        )
        return self._get(P, time, merged, interp)  # should later be a TensorField
        return RunMan(self).GetPressure(time, pop_name, merged, interp, **kwargs)

    def GetPi(self, time, merged=False, interp="nearest", **kwargs):
        M = self._get_hierarchy(time, "ions_momentum_tensor.h5", **kwargs)
        massDensity = self.GetMassDensity(time, **kwargs)
        Vi = self._get_hier_for(time, "ions_bulkVelocity", **kwargs)
        Pi = compute_hier_from(_compute_pressure, (M, massDensity, Vi))
        return self._get(Pi, time, merged, interp)  # should later be a TensorField

    def GetPe(self, time, merged=False, interp="nearest", all_primal=True):
        hier = self._get_hier_for(time, "ions_charge_density")
        Te = hier.sim.electrons.closure.Te
        if not all_primal:
            return Te * self._get(hier, time, merged, interp)
        h = compute_hier_from(hc.drop_ghosts, hier)
        return ScalarField(h) * Te

    def GetJ(self, time, merged=False, interp="nearest", all_primal=True, **kwargs):
        if merged:
            all_primal = False
        B = self.GetB(time, all_primal=False, **kwargs)
        J = compute_hier_from(_compute_current, B)
        if not all_primal:
            return self._get(J, time, merged, interp)
        J = hc.rename(compute_hier_from(_compute_to_primal, J), ["Jx", "Jy", "Jz"])
        return VectorField(J)

    def GetDivB(self, time, merged=False, interp="nearest", **kwargs):
        B = self.GetB(time, all_primal=False, **kwargs)
        db = compute_hier_from(_compute_divB, B)
        return ScalarField(self._get(db, time, merged, interp))

    def GetRanks(self, time, merged=False, interp="nearest", **kwargs):
        """
        returns a hierarchy of MPI ranks
        takes the information from magnetic field diagnostics arbitrarily
        this fails if the magnetic field is not written and could at some point
        be replace by a search of any available diag at the requested time.
        """
        B = self.GetB(time, all_primal=False, **kwargs)
        ranks = compute_hier_from(_get_rank, B)
        return ScalarField(self._get(ranks, time, merged, interp))

    def GetParticles(self, time, pop_name, hier=None, **kwargs):
        def filename(name):
            return f"ions_pop_{name}_domain.h5"

        if isinstance(pop_name, (list, tuple)):
            for pop in pop_name:
                hier = self._get_hierarchy(time, filename(pop), hier=hier, **kwargs)
            return hier
        return self._get_hierarchy(time, filename(pop_name), hier=hier, **kwargs)

    def GetParticleCount(self, time, **kwargs):
        c = self._get_hierarchy(time, "particle_count.h5", **kwargs)
        return c

    def GetMass(self, pop_name, **kwargs):
        list_of_qty = ["density", "flux", "domain", "levelGhost"]
        list_of_mass = []

        import h5py

        for qty in list_of_qty:
            file = os.path.join(self.path, "ions_pop_{}_{}.h5".format(pop_name, qty))
            if os.path.isfile(file):
                h5_file = h5py.File(file, "r")
                list_of_mass.append(h5_file.attrs["pop_mass"])

        assert all(m == list_of_mass[0] for m in list_of_mass)

        return list_of_mass[0]

    def GetDomainSize(self, **kwargs):
        hier = self._get_any_hierarchy(self.default_time)
        root_cell_width = hier.level(0).cell_width
        domain_box = hier.domain_box
        return (domain_box.upper + 1) * root_cell_width

    def GetDl(self, level="finest", time=None):
        """
        gives the ndarray containing the grid sizes at the given time
        for the hierarchy defined in the given run, and for the given level
        (default is 'finest', but can also be a int)

        :param level: the level at which get the associated grid size
        :param time: the time because level depends on it
        """

        time = time if time is not None else self.default_time
        hier = self._get_any_hierarchy(time)
        level = hier.finest_level(time) if level == "finest" else level
        return hier.level(level).cell_width

    def all_times(self):
        return {Path(file).stem: all_times_from(file) for file in self.available_diags}

    def times(self, qty):
        return self.all_times()[qty]

    def _available_diags(self):
        files = glob.glob(os.path.join(self.path, "*.h5"))
        if files:
            return files
        files = glob.glob(os.path.join(self.path, "*.vtkhdf"))
        if files:
            return files
        raise RuntimeError(f"No HDF5 files found at: {self.path}")

    def _get_hier_for(self, time, qty, **kwargs):
        path = Path(self.path) / f"{qty}.h5"
        if path.exists():
            return self._get_hierarchy(time, f"{qty}.h5", **kwargs)
        path = Path(self.path) / f"{qty}.vtkhdf"
        if path.exists():
            return self._get_hierarchy(time, f"{qty}.vtkhdf", **kwargs)
        raise RuntimeError(f"No HDF5 file found for: {qty}")

    def _get_any_hierarchy(self, time):
        ref_file = Path(self.available_diags[0]).stem
        return self._get_hier_for(time, ref_file)

    def _get_hierarchy(self, times, filename, hier=None, **kwargs):
        from pyphare.core.box import Box

        times = phut.listify(times)
        times = [f"{t:.10f}" for t in times]
        if "selection_box" in kwargs:
            if isinstance(kwargs["selection_box"], tuple):
                lower = kwargs["selection_box"][:2]
                upper = kwargs["selection_box"][2:]
                kwargs["selection_box"] = Box(lower, upper)

        return hierarchy_from(
            h5_filename=os.path.join(self.path, filename),
            times=times,
            hier=hier,
            **kwargs,
        )

    # TODO maybe transform that so multiple times can be accepted
    def _get(self, hierarchy, time, merged, interp, drop_ghosts=False):
        """
        if merged=True, will return an interpolator and a tuple of 1d arrays
        with the coordinates of the finest grid where the interpolator
        can be calculated (that is the return of flat_finest_field)
        """
        if merged:
            domain = self.GetDomainSize()
            dl = self.GetDl(time=time)

            # assumes all qties in the hierarchy have the same ghost width
            # so take the first patch data of the first patch of the first level....
            nbrGhosts = list(hierarchy.level(0).patches[0].patch_datas.values())[
                0
            ].ghosts_nbr
            merged_qties = {}
            for qty in hierarchy.quantities():
                data, coords = flat_finest_field(hierarchy, qty, time=time)
                merged_qties[qty] = make_interpolator(
                    data, coords, interp, domain, dl, qty, nbrGhosts
                )
            return merged_qties
        else:
            return (
                compute_hier_from(hc.drop_ghosts, hierarchy)
                if drop_ghosts
                else hierarchy
            )

    @property
    def default_time(self):
        if self.default_time_ is None:
            ref_file = self.available_diags[0]
            self.default_time_ = default_time_from(ref_file)
        return self.default_time_
