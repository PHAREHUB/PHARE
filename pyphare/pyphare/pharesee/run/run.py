#
#
#

import os
import glob
from pathlib import Path

from pyphare.pharesee.hierarchy import all_times_from
from pyphare.pharesee.hierarchy import default_time_from
from pyphare.pharesee.hierarchy import hierarchy_from
from pyphare.pharesee.hierarchy import compute as hc
from pyphare.pharesee.hierarchy import ScalarField, VectorField

from pyphare.core import phare_utilities as phut
from pyphare.pharesee.hierarchy.hierarchy_utils import compute_hier_from

from pyphare.logger import getLogger

from .man import RunMan

from .utils import (
    _compute_to_primal,
    _compute_pop_pressure,
    _compute_pressure,
    _compute_current,
    _compute_divB,
    _get_rank,
)

logger = getLogger(__name__)


class Run:
    def __init__(self, path, default_time=None):
        self.path = path
        self.default_time_ = default_time
        self.available_diags = self._available_diags()

    def GetTags(self, time, **kwargs):
        hier = self._get_hierarchy(time, "tags.h5")
        return self._get(hier)

    def GetB(self, time, all_primal=True, **kwargs):
        return RunMan(self).GetB(time, all_primal=all_primal, **kwargs)

    def GetE(self, time, all_primal=True, **kwargs):
        hier = self._get_hier_for(time, "EM_E", **kwargs)
        if not all_primal:
            return self._get(hier)
        return VectorField.FROM(compute_hier_from(_compute_to_primal, hier))

    def GetMassDensity(self, time, **kwargs):
        hier = self._get_hier_for(time, "ions_mass_density", **kwargs)
        return ScalarField.FROM(self._get(hier))

    def GetNi(self, time, **kwargs):
        hier = self._get_hier_for(time, "ions_charge_density", **kwargs)
        return ScalarField.FROM(self._get(hier, drop_ghosts=True))

    def GetN(self, time, pop_name, **kwargs):
        hier = self._get_hier_for(time, f"ions_pop_{pop_name}_density", **kwargs)
        return ScalarField.FROM(self._get(hier))

    def GetVi(self, time, **kwargs):
        return RunMan(self).GetVi(time, **kwargs)

    def GetFlux(self, time, pop_name, **kwargs):
        hier = self._get_hier_for(time, f"ions_pop_{pop_name}_flux", **kwargs)
        return VectorField.FROM(self._get(hier))

    def GetPressure(self, time, pop_name, **kwargs):
        M = self._get_hier_for(time, f"ions_pop_{pop_name}_momentum_tensor", **kwargs)
        V = self.GetFlux(time, pop_name, **kwargs)
        N = self.GetN(time, pop_name, **kwargs)
        P = compute_hier_from(
            _compute_pop_pressure,
            (M, V, N),
            popname=pop_name,
            mass=self.GetMass(pop_name, **kwargs),
        )
        return self._get(P)  # should later be a TensorField

    def GetPi(self, time, **kwargs):
        M = self._get_hier_for(time, "ions_momentum_tensor", **kwargs)
        massDensity = self.GetMassDensity(time, **kwargs)
        Vi = self._get_hier_for(time, "ions_bulkVelocity", **kwargs)
        Pi = compute_hier_from(_compute_pressure, (M, massDensity, Vi))
        return self._get(Pi)  # should later be a TensorField

    def GetPe(self, time, all_primal=True):
        hier = self._get_hier_for(time, "ions_charge_density")
        Te = hier.sim.electrons.closure.Te
        if not all_primal:
            return Te * self._get(hier)
        h = compute_hier_from(hc.drop_ghosts, hier)
        return ScalarField.FROM(h) * Te

    def GetJ(self, time, all_primal=True, **kwargs):
        B = self.GetB(time, all_primal=False, **kwargs)
        J = compute_hier_from(_compute_current, B)
        if not all_primal:
            return self._get(J)
        return VectorField.FROM(compute_hier_from(_compute_to_primal, J))

    def GetDivB(self, time, **kwargs):
        B = self.GetB(time, all_primal=False, **kwargs)
        db = compute_hier_from(_compute_divB, B)
        return ScalarField.FROM(self._get(db))

    def GetRanks(self, time, **kwargs):
        """
        returns a hierarchy of MPI ranks
        takes the information from magnetic field diagnostics arbitrarily
        this fails if the magnetic field is not written and could at some point
        be replace by a search of any available diag at the requested time.
        """
        B = self.GetB(time, all_primal=False, **kwargs)
        ranks = compute_hier_from(_get_rank, B)
        return ScalarField.FROM(self._get(ranks))

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

    def _get(self, hierarchy, drop_ghosts=False):
        return (
            compute_hier_from(hc.drop_ghosts, hierarchy) if drop_ghosts else hierarchy
        )

    @property
    def default_time(self):
        if self.default_time_ is None:
            ref_file = self.available_diags[0]
            self.default_time_ = default_time_from(ref_file)
        return self.default_time_
