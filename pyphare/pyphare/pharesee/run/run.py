import os
import glob
import numpy as np

from pyphare.pharesee.hierarchy import hierarchy_from
from pyphare.pharesee.hierarchy import ScalarField, VectorField

from pyphare.pharesee.hierarchy.hierarchy_utils import compute_hier_from
from pyphare.pharesee.hierarchy.hierarchy_utils import flat_finest_field
from pyphare.pharesee.hierarchy import hierarchy_compute as hc
from pyphare.core.phare_utilities import listify

from pyphare.logger import getLogger
from .utils import (
    _compute_to_primal,
    _compute_pop_pressure,
    _compute_pressure,
    _compute_current,
    _compute_divB,
    _get_rank,
    make_interpolator,
    finest_coords_for,
)


logger = getLogger(__name__)

quantities_per_file = {
    "EM_B": "B",
    "EM_E": "E",
    "ions_bulkVelocity": "Vi",
    "ions_charge_density": "Ni",
    "particle_count": "nppc",
    "mhd_rho": "mhd_rho",
    "mhd_V": "mhd_V",
    "mhd_P": "mhd_P",
    "mhd_rhoV": "mhd_rhoV",
    "mhd_Etot": "mhd_Etot",
}


class Run:
    def __init__(self, path, default_time=None):
        self.path = path
        self.default_time_ = default_time
        self.available_diags = glob.glob(os.path.join(self.path, "*.h5"))

    def _get_hierarchy(self, times, filename, hier=None, **kwargs):
        from pyphare.core.box import Box

        times = listify(times)
        times = [f"{t:.10f}" for t in times]
        if "selection_box" in kwargs:
            if isinstance(kwargs["selection_box"], tuple):
                lower = kwargs["selection_box"][:2]
                upper = kwargs["selection_box"][2:]
                kwargs["selection_box"] = Box(lower, upper)

        def _get_hier(h):
            return hierarchy_from(
                h5_filename=os.path.join(self.path, filename),
                times=times,
                hier=h,
                **kwargs,
            )

        return _get_hier(hier)

    # TODO maybe transform that so multiple times can be accepted
    def _get(self, hierarchy, time, merged, interp, drop_ghosts=False):
        """
        if merged=True, returns {qty: (interpolator, finest_coords)} where
        interpolator is a callable scipy interpolator and finest_coords is a
        tuple of 1d coordinate arrays at the finest grid resolution.
        """
        if merged:
            domain = self.GetDomainSize()
            dl = self.GetDl(time=time)
            nbrGhosts = list(hierarchy.level(0).patches[0].patch_datas.values())[
                0
            ].ghosts_nbr
            merged_qties = {}
            for qty in hierarchy.quantities():
                data, coords = flat_finest_field(hierarchy, qty, time=time)
                merged_qties[qty] = (
                    make_interpolator(data, coords, interp),
                    finest_coords_for(domain, dl, qty, nbrGhosts),
                )
            return merged_qties
        else:
            return compute_hier_from(hc.drop_ghosts, hierarchy) if drop_ghosts else hierarchy

    def GetTags(self, time, merged=False, **kwargs):
        hier = self._get_hierarchy(time, "tags.h5")
        return self._get(hier, time, merged, "nearest")

    def GetB(self, time, merged=False, interp="nearest", all_primal=True, **kwargs):
        if merged:
            all_primal = False
        hier = self._get_hierarchy(time, "EM_B.h5", **kwargs)
        if not all_primal:
            return self._get(hier, time, merged, interp)

        h = compute_hier_from(_compute_to_primal, hier, x="Bx", y="By", z="Bz")
        return VectorField.FROM(h)

    def GetE(self, time, merged=False, interp="nearest", all_primal=True, **kwargs):
        if merged:
            all_primal = False
        hier = self._get_hierarchy(time, "EM_E.h5", **kwargs)
        if not all_primal:
            return self._get(hier, time, merged, interp)

        h = compute_hier_from(_compute_to_primal, hier, x="Ex", y="Ey", z="Ez")
        return VectorField.FROM(h)

    def GetMassDensity(self, time, merged=False, interp="nearest", **kwargs):
        hier = self._get_hierarchy(time, "ions_mass_density.h5", **kwargs)
        return ScalarField.FROM(self._get(hier, time, merged, interp))

    def GetNi(self, time, merged=False, interp="nearest", **kwargs):
        hier = self._get_hierarchy(time, "ions_charge_density.h5", **kwargs)
        return ScalarField.FROM(self._get(hier, time, merged, interp, drop_ghosts=True))

    def GetN(self, time, pop_name, merged=False, interp="nearest", **kwargs):
        hier = self._get_hierarchy(time, f"ions_pop_{pop_name}_density.h5", **kwargs)
        return ScalarField.FROM(self._get(hier, time, merged, interp))

    def GetVi(self, time, merged=False, interp="nearest", **kwargs):
        hier = self._get_hierarchy(time, "ions_bulkVelocity.h5", **kwargs)
        return VectorField.FROM(self._get(hier, time, merged, interp, drop_ghosts=True))

    def GetFlux(self, time, pop_name, merged=False, interp="nearest", **kwargs):
        hier = self._get_hierarchy(time, f"ions_pop_{pop_name}_flux.h5", **kwargs)
        return VectorField.FROM(self._get(hier, time, merged, interp))

    def GetPressure(self, time, pop_name, merged=False, interp="nearest", **kwargs):
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

    def GetPi(self, time, merged=False, interp="nearest", **kwargs):
        M = self._get_hierarchy(time, "ions_momentum_tensor.h5", **kwargs)
        massDensity = self.GetMassDensity(time, **kwargs)
        Vi = self._get_hierarchy(time, "ions_bulkVelocity.h5", **kwargs)
        Pi = compute_hier_from(_compute_pressure, (M, massDensity, Vi))
        return self._get(Pi, time, merged, interp)  # should later be a TensorField

    def GetPe(self, time, merged=False, interp="nearest", all_primal=True):
        hier = self._get_hierarchy(time, "ions_charge_density.h5")

        Te = hier.sim.electrons.closure.Te

        if not all_primal:
            return Te * self._get(hier, time, merged, interp)

        h = compute_hier_from(hc.drop_ghosts, hier)
        return ScalarField.FROM(h) * Te

    def GetJ(self, time, merged=False, interp="nearest", all_primal=True, **kwargs):
        if merged:
            all_primal = False
        B = self.GetB(time, all_primal=False, **kwargs)
        J = compute_hier_from(_compute_current, B)
        if not all_primal:
            return self._get(J, time, merged, interp)
        h = compute_hier_from(_compute_to_primal, J, x="Jx", y="Jy", z="Jz")
        return VectorField.FROM(h)

    def GetDivB(self, time, merged=False, interp="nearest", **kwargs):
        B = self.GetB(time, all_primal=False, **kwargs)
        db = compute_hier_from(_compute_divB, B)
        return ScalarField.FROM(self._get(db, time, merged, interp))

    def GetMHDrho(
        self, time, merged=False, interp="nearest", all_primal=True, **kwargs
    ):
        if merged:
            all_primal = False
        hier = self._get_hierarchy(time, "mhd_rho.h5", **kwargs)
        if not all_primal:
            return self._get(hier, time, merged, interp)

        h = compute_hier_from(_compute_to_primal, hier, value="mhdRho")
        return ScalarField.FROM(h)

    def GetMHDV(self, time, merged=False, interp="nearest", all_primal=True, **kwargs):
        if merged:
            all_primal = False
        hier = self._get_hierarchy(time, "mhd_V.h5", **kwargs)
        if not all_primal:
            return self._get(hier, time, merged, interp)

        h = compute_hier_from(_compute_to_primal, hier, x="mhdVx", y="mhdVy", z="mhdVz")
        return VectorField.FROM(h)

    def GetMHDP(self, time, merged=False, interp="nearest", all_primal=True, **kwargs):
        if merged:
            all_primal = False
        hier = self._get_hierarchy(time, "mhd_P.h5", **kwargs)
        if not all_primal:
            return self._get(hier, time, merged, interp)

        h = compute_hier_from(_compute_to_primal, hier, value="mhdP")
        return ScalarField.FROM(h)

    def GetMHDrhoV(
        self, time, merged=False, interp="nearest", all_primal=True, **kwargs
    ):
        if merged:
            all_primal = False
        hier = self._get_hierarchy(time, "mhd_rhoV.h5", **kwargs)
        if not all_primal:
            return self._get(hier, time, merged, interp)

        h = compute_hier_from(
            _compute_to_primal, hier, x="mhdRhoVx", y="mhdRhoVy", z="mhdRhoVz"
        )
        return VectorField.FROM(h)

    def GetMHDEtot(
        self, time, merged=False, interp="nearest", all_primal=True, **kwargs
    ):
        if merged:
            all_primal = False
        hier = self._get_hierarchy(time, "mhd_Etot.h5", **kwargs)
        if not all_primal:
            return self._get(hier, time, merged, interp)

        h = compute_hier_from(_compute_to_primal, hier, value="mhdEtot")
        return ScalarField.FROM(h)

    def GetMagneticFlux(
        self, time, interp="nearest", xn=None, yn=None, Xn=None, Yn=None
    ):
        # Reuse grids if provided, otherwise generate them
        if xn is None or yn is None or Xn is None or Yn is None:
            domain = self.GetDomainSize()
            dl = self.GetDl(level="finest", time=time)
            xn = np.arange(0, domain[0] + dl[0], dl[0])
            yn = np.arange(0, domain[1] + dl[1], dl[1])
            Xn, Yn = np.meshgrid(xn, yn, indexing="ij")

        merged_B = self.GetB(time, merged=True, interp=interp)
        bx_interp = merged_B["Bx"][0]
        by_interp = merged_B["By"][0]

        bx = bx_interp(Xn, Yn)
        by = by_interp(Xn, Yn)

        from scipy.integrate import cumulative_trapezoid

        Az_x0 = -cumulative_trapezoid(by[:, 0], xn, initial=0)
        Az = cumulative_trapezoid(bx, yn, axis=1, initial=0)
        Az += Az_x0[:, np.newaxis]

        return Az, (xn, yn)

    def FindPrimaryXPoint(self, Az, xn, yn):
        dAz_dx = np.gradient(Az, xn, axis=0)
        dAz_dy = np.gradient(Az, yn, axis=1)

        grad_mag = np.sqrt(dAz_dx**2 + dAz_dy**2)
        threshold = np.percentile(grad_mag, 5)

        d2Az_dx2 = np.gradient(dAz_dx, xn, axis=0)
        d2Az_dy2 = np.gradient(dAz_dy, yn, axis=1)
        d2Az_dxdy = np.gradient(dAz_dx, yn, axis=1)

        det_hessian = d2Az_dx2 * d2Az_dy2 - d2Az_dxdy**2

        candidates = (grad_mag < threshold) & (det_hessian < 0)

        det_hessian_masked = np.where(candidates, det_hessian, np.inf)
        idx = np.unravel_index(np.argmin(det_hessian_masked), Az.shape)

        return xn[idx[0]], yn[idx[1]], idx

    def GetReconnectionRate(self, times, interp="nearest"):
        domain = self.GetDomainSize()
        dl = self.GetDl(level="finest", time=times[0])
        xn = np.arange(0, domain[0] + dl[0], dl[0])
        yn = np.arange(0, domain[1] + dl[1], dl[1])
        Xn, Yn = np.meshgrid(xn, yn, indexing="ij")

        flux_at_xpoint = []
        xpoint_trajectory = []

        for t in times:
            Az, _ = self.GetMagneticFlux(t, interp=interp, xn=xn, yn=yn, Xn=Xn, Yn=Yn)

            x_xp, y_xp, idx = self.FindPrimaryXPoint(Az, xn, yn)

            xpoint_trajectory.append([x_xp, y_xp])

            flux_at_xpoint.append(Az[idx])

        flux_at_xpoint = np.array(flux_at_xpoint)
        xpoint_trajectory = np.array(xpoint_trajectory)

        dt = np.diff(times)
        rates = np.diff(flux_at_xpoint) / dt
        times_centered = (times[:-1] + times[1:]) / 2

        return times_centered, rates, flux_at_xpoint, xpoint_trajectory

    def GetRanks(self, time, merged=False, interp="nearest", **kwargs):
        """
        returns a hierarchy of MPI ranks
        takes the information from magnetic field diagnostics arbitrarily
        this fails if the magnetic field is not written and could at some point
        be replace by a search of any available diag at the requested time.
        """
        B = self.GetB(time, all_primal=False, **kwargs)
        ranks = compute_hier_from(_get_rank, B)
        return ScalarField.FROM(self._get(ranks, time, merged, interp))

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
        import h5py

        # That is the first file in the available diags
        data_file = h5py.File(self.available_diags[0], "r")
        root_cell_width = np.asarray(data_file.attrs["cell_width"])

        return (data_file.attrs["domain_box"] + 1) * root_cell_width

    def GetDl(self, level="finest", time=None):
        """
        gives the ndarray containing the grid sizes at the given time
        for the hierarchy defined in the given run, and for the given level
        (default is 'finest', but can also be a int)

        :param level: the level at which get the associated grid size
        :param time: the time because level depends on it
        """

        import h5py

        def _get_time():
            if time:
                return time
            if self.default_time_:
                return self.default_time_
            self.default_time_ = float(list(data_file[h5_time_grp_key].keys())[0])
            return self.default_time_

        h5_time_grp_key = "t"
        files = self.available_diags

        for h5_filename in files:
            data_file = h5py.File(h5_filename, "r")

            time = _get_time()

            try:
                hier = self._get_hierarchy(time, h5_filename.split("/")[-1])

                if level == "finest":
                    level = hier.finest_level(time)
                fac = np.power(hier.refinement_ratio, level)

                root_cell_width = np.asarray(data_file.attrs["cell_width"])

                return root_cell_width / fac

            except KeyError:
                ...  # time may not be available for given quantity

        raise RuntimeError("Unable toGetDl")

    def all_times(self):
        import h5py

        files = self.available_diags
        ts = {}
        for file in files:
            basename = os.path.basename(file).split(".")[0]
            ff = h5py.File(file)
            time_keys = ff["t"].keys()
            time = np.zeros(len(time_keys))
            for it, t in enumerate(time_keys):
                time[it] = float(t)
            if basename in quantities_per_file:
                ts[quantities_per_file[basename]] = time
            ff.close()
        return ts

    def times(self, qty):
        return self.all_times()[qty]
