import glob
import os
from pathlib import Path

import numpy as np

from pyphare.core import phare_utilities as phut
from pyphare.logger import getLogger
from pyphare.pharesee.hierarchy import (
    ScalarField,
    VectorField,
    all_times_from,
    default_time_from,
    hierarchy_from,
)
from pyphare.pharesee.hierarchy import hierarchy_compute as hc
from pyphare.pharesee.hierarchy.fromh5 import pop_attribute_from_h5
from pyphare.pharesee.hierarchy.hierarchy_utils import compute_hier_from

from .utils import (
    _compute_current,
    _compute_divB,
    _compute_pop_pressure,
    _compute_pressure,
    _compute_to_primal,
    _get_rank,
)

logger = getLogger(__name__)


class Run:
    def __init__(self, path, default_time=None):
        self.path = path
        self.default_time_ = default_time
        self.available_diags = self._available_diags()

    def GetTags(self, time, **kwargs):
        return self._get_hierarchy(time, "tags.h5")

    def GetB(self, time, all_primal=True, **kwargs):
        hier = self._get_hierarchy(time, "EM_B.h5", **kwargs)
        if not all_primal:
            return hier

        return VectorField.FROM(compute_hier_from(_compute_to_primal, hier))

    def GetE(self, time, all_primal=True, **kwargs):
        hier = self._get_hierarchy(time, "EM_E.h5", **kwargs)
        if not all_primal:
            return hier

        return VectorField.FROM(compute_hier_from(_compute_to_primal, hier))

    def GetMassDensity(self, time, **kwargs):
        hier = self._get_hierarchy(time, "ions_mass_density.h5", **kwargs)
        return ScalarField.FROM(compute_hier_from(hc.drop_ghosts, hier))

    def GetNi(self, time, **kwargs):
        hier = self._get_hierarchy(time, "ions_charge_density.h5", **kwargs)
        return ScalarField.FROM(compute_hier_from(hc.drop_ghosts, hier))

    def GetN(self, time, pop_name, **kwargs):
        hier = self._get_hierarchy(time, f"ions_pop_{pop_name}_density.h5", **kwargs)
        return ScalarField.FROM(compute_hier_from(hc.drop_ghosts, hier))

    def GetVi(self, time, **kwargs):
        hier = self._get_hierarchy(time, "ions_bulkVelocity.h5", **kwargs)
        return VectorField.FROM(compute_hier_from(hc.drop_ghosts, hier))

    def GetChargeFlux(self, time, **kwargs):
        hier = self._get_hierarchy(time, "ions_charge_flux.h5", **kwargs)
        return VectorField.FROM(hier)

    def GetFlux(self, time, pop_name, **kwargs):
        hier = self._get_hierarchy(time, f"ions_pop_{pop_name}_flux.h5", **kwargs)
        return VectorField.FROM(compute_hier_from(hc.drop_ghosts, hier))

    def GetPressure(self, time, pop_name, **kwargs):
        M = self._get_hierarchy(
            time, f"ions_pop_{pop_name}_momentum_tensor.h5", **kwargs
        )
        M = compute_hier_from(hc.drop_ghosts, M)
        V = self.GetFlux(time, pop_name, **kwargs)
        N = self.GetN(time, pop_name, **kwargs)
        P = compute_hier_from(
            _compute_pop_pressure,
            (M, V, N),
            popname=pop_name,
            mass=self.GetMass(pop_name, **kwargs),
        )
        return P  # should later be a TensorField

    def GetPi(self, time, **kwargs):
        M = self._get_hierarchy(time, "ions_momentum_tensor.h5", **kwargs)
        M = compute_hier_from(hc.drop_ghosts, M)
        massDensity = self.GetMassDensity(time, **kwargs)
        Vi = self._get_hierarchy(time, "ions_bulkVelocity.h5", **kwargs)
        Vi = compute_hier_from(hc.drop_ghosts, Vi)
        Pi = compute_hier_from(_compute_pressure, (M, massDensity, Vi))
        return Pi  # should later be a TensorField

    def GetPe(self, time, all_primal=True):
        hier = self._get_hierarchy(time, "ions_charge_density.h5")
        Te = hier.sim.electrons.closure.Te
        if not all_primal:
            return Te * hier

        h = compute_hier_from(hc.drop_ghosts, hier)
        return ScalarField.FROM(h) * Te

    def GetJ(self, time, all_primal=True, **kwargs):
        B = self.GetB(time, all_primal=False, **kwargs)
        J = compute_hier_from(_compute_current, B)
        if not all_primal:
            return J
        h = compute_hier_from(_compute_to_primal, J)
        return VectorField.FROM(h)

    def GetDivB(self, time, **kwargs):
        B = self.GetB(time, all_primal=False, **kwargs)
        db = compute_hier_from(_compute_divB, B)
        return ScalarField.FROM(db)

    def GetMHDrho(self, time, all_primal=True, **kwargs):
        hier = self._get_hierarchy(time, "mhd_rho.h5", **kwargs)
        if not all_primal:
            return hier

        h = compute_hier_from(_compute_to_primal, hier, value="mhdRho")
        return ScalarField.FROM(h)

    def GetMHDV(self, time, all_primal=True, **kwargs):
        hier = self._get_hierarchy(time, "mhd_V.h5", **kwargs)
        if not all_primal:
            return hier

        h = compute_hier_from(_compute_to_primal, hier, x="mhdVx", y="mhdVy", z="mhdVz")
        return VectorField.FROM(h)

    def GetMHDP(self, time, all_primal=True, **kwargs):
        hier = self._get_hierarchy(time, "mhd_P.h5", **kwargs)
        if not all_primal:
            return hier

        h = compute_hier_from(_compute_to_primal, hier, value="mhdP")
        return ScalarField.FROM(h)

    def GetMHDrhoV(self, time, all_primal=True, **kwargs):
        hier = self._get_hierarchy(time, "mhd_rhoV.h5", **kwargs)
        if not all_primal:
            return hier

        h = compute_hier_from(
            _compute_to_primal, hier, x="mhdRhoVx", y="mhdRhoVy", z="mhdRhoVz"
        )
        return VectorField.FROM(h)

    def GetMHDEtot(self, time, all_primal=True, **kwargs):
        hier = self._get_hierarchy(time, "mhd_Etot.h5", **kwargs)
        if not all_primal:
            return hier

        h = compute_hier_from(_compute_to_primal, hier, value="mhdEtot")
        return ScalarField.FROM(h)

    def GetMagneticFlux(self, time, interp="nearest"):
        B = self.GetB(time)  # all primal, on the nodes of the finest grid
        bx = B.finest(time, "x", interp)
        xn, yn = bx.x, bx.y
        bx = bx[:]
        by = B.finest(time, "y", interp)[:]

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
        flux_at_xpoint = []
        xpoint_trajectory = []

        for t in times:
            Az, (xn, yn) = self.GetMagneticFlux(t, interp=interp)

            x_xp, y_xp, idx = self.FindPrimaryXPoint(Az, xn, yn)

            xpoint_trajectory.append([x_xp, y_xp])

            flux_at_xpoint.append(Az[idx])

        flux_at_xpoint = np.array(flux_at_xpoint)
        xpoint_trajectory = np.array(xpoint_trajectory)

        dt = np.diff(times)
        rates = np.diff(flux_at_xpoint) / dt
        times_centered = (times[:-1] + times[1:]) / 2

        return times_centered, rates, flux_at_xpoint, xpoint_trajectory

    def GetRanks(self, time, **kwargs):
        """
        returns a hierarchy of MPI ranks
        takes the information from magnetic field diagnostics arbitrarily
        this fails if the magnetic field is not written and could at some point
        be replace by a search of any available diag at the requested time.
        """
        B = self.GetB(time, all_primal=False, **kwargs)
        ranks = compute_hier_from(_get_rank, B)
        return ScalarField.FROM(ranks)

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
        return pop_attribute_from_h5(self.path, pop_name, "pop_mass")

    def GetCharge(self, pop_name, **kwargs):
        return pop_attribute_from_h5(self.path, pop_name, "pop_charge")

    def GetDomainSize(self, **kwargs):
        hier = self._get_any_hierarchy(self.default_time)
        root_cell_width = hier.level(0).cell_width
        domain_box = hier.domain_box
        return (domain_box.upper + 1) * root_cell_width

    def GetDl(self, level=0):
        """
        gives the ndarray containing the grid sizes for the given level

        :param level: the level at which get the associated grid size
        """
        if level < 0:
            raise ValueError(f"Run::GetDl: level must be >= 0, got {level}")

        hier = self._get_any_hierarchy(self.default_time)
        return hier.level(0).cell_width / phut.refinement_ratio**level

    def all_times(self):
        return {Path(file).stem: all_times_from(file) for file in self.available_diags}

    def times(self, qty):
        return self.all_times()[qty]

    def _available_diags(self):
        files = glob.glob(os.path.join(self.path, "*.h5"))
        if not files:
            raise RuntimeError(f"No HDF5 files found at: {self.path}")
        return files

    def _get_any_hierarchy(self, time):
        for file in self.available_diags:
            if f"{time:.10f}" in (f"{t:.10f}" for t in all_times_from(file)):
                return self._get_hierarchy(time, os.path.basename(file))
        raise RuntimeError(f"No diagnostic found at time {time} in {self.path}")

    def _get_hierarchy(self, times, filename, hier=None, **kwargs):
        from pyphare.core.box import Box

        times = phut.listify(times)
        times = [f"{t:.10f}" for t in times]
        if "selection_box" in kwargs and isinstance(kwargs["selection_box"], tuple):
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

    @property
    def default_time(self):
        if self.default_time_ is None:
            ref_file = self.available_diags[0]
            self.default_time_ = default_time_from(ref_file)
        return self.default_time_
