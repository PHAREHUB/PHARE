"""Composite AMR error norm for space-time convergence tests.

Error norm in the style of Gardiner & Stone (JCP 227, 2008) / Stone et al.
(ApJS 178, 2008), on the Berger-Colella composite grid as in McCorquodale &
Colella (CAMCoS 6, 2011): for each conserved variable q in {rho, rhoVx, rhoVy,
rhoVz, Etot, Bx, By, Bz}, compute the volume-weighted L1 error over the
composite grid (fine level everywhere it exists, coarse level only on its
uncovered region), then report eps = sqrt(sum_q dq^2) ("L2-of-L1").

The comparison is t=final_time against the t=0 dump of the SAME run, patch by
patch (identical discrete representation -> no staggered-sampling error); this
requires the exact solution at final_time to equal the initial condition (wave
period / traversal time). All variables are brought to cell centers (B via
2-point face average) so covered-cell masking is exact and shapes are uniform.
Runs must dump the conserved set natively (rho, rhoV, Etot + B).
"""

import numpy as np

CONSERVED = ["rho", "rhoVx", "rhoVy", "rhoVz", "Etot", "Bx", "By", "Bz"]


def interior(pd):
    ng = pd.ghosts_nbr
    sl = tuple(slice(g, -g if g else None) for g in ng)
    return np.asarray(pd.dataset[:])[sl]


def to_center(arr, axis_primal):
    """Average a face-centered array to cell centers along its primal axis."""
    if axis_primal is None:
        return arr
    lo = [slice(None)] * arr.ndim
    hi = [slice(None)] * arr.ndim
    lo[axis_primal] = slice(None, -1)
    hi[axis_primal] = slice(1, None)
    return 0.5 * (arr[tuple(lo)] + arr[tuple(hi)])


def patchwise(hier, time, names):
    """{key: {name: interior array}} with key = (ilvl, lower, upper)."""
    out = {}
    for ilvl, level in hier.levels(time).items():
        for patch in level.patches:
            key = (ilvl, tuple(patch.box.lower), tuple(patch.box.upper))
            d = out.setdefault(key, {})
            for name in names:
                d[name] = interior(patch.patch_datas[name])
    return out


def conserved_fields(run, time):
    """{key: {var: cell-centered array}} for the 8 conserved variables."""
    b = patchwise(run.GetB(time, all_primal=False), time, ["Bx", "By", "Bz"])
    rho = patchwise(run.GetMHDrho(time, all_primal=False), time, ["mhdRho"])
    rhoV = patchwise(
        run.GetMHDrhoV(time, all_primal=False), time,
        ["mhdRhoVx", "mhdRhoVy", "mhdRhoVz"],
    )
    etot = patchwise(run.GetMHDEtot(time, all_primal=False), time, ["mhdEtot"])

    fields = {}
    for key in rho:
        bc = {
            "Bx": to_center(b[key]["Bx"], 0),
            "By": to_center(b[key]["By"], 1),
            "Bz": b[key]["Bz"],
        }
        f = {"rho": rho[key]["mhdRho"], **bc}
        for c, name in zip("xyz", ["mhdRhoVx", "mhdRhoVy", "mhdRhoVz"]):
            f[f"rhoV{c}"] = rhoV[key][name]
        f["Etot"] = etot[key]["mhdEtot"]
        fields[key] = f
    return fields


def valid_mask(key, fine_box, ndim=2):
    """1 where this patch's cells belong to the composite grid, 0 where a coarse
    cell is covered by the fine level. fine_box = [lower, upper] in L0 cells."""
    ilvl, lower, upper = key
    shape = tuple(up - lo + 1 for lo, up in zip(lower, upper))
    mask = np.ones(shape)
    if ilvl == 0:
        sl = []
        for d in range(ndim):
            lo = max(fine_box[0][d], lower[d]) - lower[d]
            hi = min(fine_box[1][d], upper[d]) - lower[d]
            if hi < lo:
                return mask
            sl.append(slice(lo, hi + 1))
        mask[tuple(sl)] = 0.0
    return mask


def composite_errors(run, final_time, fine_box, refinement_ratio=2, ndim=2):
    """(per_level eps, composite eps) -- L2-of-L1 over conserved vars."""
    f1 = conserved_fields(run, final_time)
    f0 = conserved_fields(run, 0.0)
    if f1.keys() != f0.keys():
        raise RuntimeError(
            f"patch layout changed between t=0 and t={final_time}: "
            f"{sorted(f1.keys() ^ f0.keys())}"
        )
    comp_num = {v: 0.0 for v in CONSERVED}
    comp_den = 0.0
    lvl_num, lvl_den = {}, {}
    for key in f1:
        ilvl = key[0]
        mask = valid_mask(key, fine_box, ndim)
        w = float(refinement_ratio) ** (-ndim * ilvl)  # relative cell volume
        ln = lvl_num.setdefault(ilvl, {v: 0.0 for v in CONSERVED})
        for v in CONSERVED:
            diff = np.abs(f1[key][v] - f0[key][v])
            if diff.shape != mask.shape:
                raise RuntimeError(f"{v}: field shape {diff.shape} != mask {mask.shape}")
            comp_num[v] += (diff * mask).sum() * w
            ln[v] += diff.sum()
        comp_den += mask.sum() * w
        lvl_den[ilvl] = lvl_den.get(ilvl, 0) + mask.size

    def eps(num, den):
        return float(np.sqrt(sum((num[v] / den) ** 2 for v in CONSERVED)))

    per_level = {ilvl: eps(lvl_num[ilvl], lvl_den[ilvl]) for ilvl in lvl_num}
    return per_level, eps(comp_num, comp_den)


def uniform_error(run, final_time):
    """eps for a single-level uniform run (no masking, no volume weights)."""
    f1 = conserved_fields(run, final_time)
    f0 = conserved_fields(run, 0.0)
    if f1.keys() != f0.keys():
        raise RuntimeError(f"patch layout changed between t=0 and t={final_time}")

    num = {v: 0.0 for v in CONSERVED}
    den = 0.0
    for key in f1:
        if key[0] != 0:
            raise RuntimeError(f"uniform run has a fine level: {key}")
        for v in CONSERVED:
            num[v] += np.abs(f1[key][v] - f0[key][v]).sum()
        den += f1[key][CONSERVED[0]].size
    return float(np.sqrt(sum((num[v] / den) ** 2 for v in CONSERVED)))
