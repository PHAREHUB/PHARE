"""
Magnetic initial condition given as a vector potential A instead of B.

B is computed in C++ as the discrete (Yee) curl of A sampled at the edges, so the
initial discrete div B is at round-off.
"""

import numpy as np

from pyphare.core.box import Box
from pyphare.core.gridlayout import HybridGridLayoutFor, MHDGridLayoutFor
from pyphare.core.gridlayout import mhdGhostNbrFromReconstruction

_components = ("x", "y", "z")
_directions = ("X", "Y", "Z")


def resolve_magnetic_init(ndim, defaulter, bx, by, bz, ax, ay, az):
    """
    Decide between the B mode (``bx, by, bz``) and the A mode (``ax, ay, az``).

    B mode, when no ``a*`` is given: returns ``{"mode": "b", "bx", "by", "bz"}``
    with the defaults (1, 0, 0) for missing components.

    A mode: only components of B along non-simulated directions may be given directly
    (``bz`` in 2D, none in 3D), and ``bz`` excludes ``ax``/``ay`` in 2D since those
    feed Bz. Missing ``a*`` default to 0. Returns
    ``{"mode": "a", "ax", "ay", "az", "direct": {"z": bz}}``, ``direct`` being
    non-empty only for a 2D ``bz``.
    """
    a_given = {c: f is not None for c, f in zip(_components, (ax, ay, az))}
    b_given = {c: f is not None for c, f in zip(_components, (bx, by, bz))}

    if not any(a_given.values()):
        return {
            "mode": "b",
            "bx": defaulter(bx, 1.0),
            "by": defaulter(by, 0.0),
            "bz": defaulter(bz, 0.0),
        }

    if ndim == 1:
        raise ValueError(
            "Vector potential (ax, ay, az) is not supported in 1D, give bx, by, bz"
        )

    in_plane = [c for c in _components[:ndim] if b_given[c]]
    if in_plane:
        raise ValueError(
            f"b{', b'.join(in_plane)} cannot be given with a vector potential in "
            f"{ndim}D: B along simulated directions comes from the curl of A"
        )

    if ndim == 2 and b_given["z"] and (a_given["x"] or a_given["y"]):
        raise ValueError(
            "bz cannot be given together with ax or ay in 2D: both define Bz"
        )

    resolved = {
        "mode": "a",
        "ax": defaulter(ax, 0.0),
        "ay": defaulter(ay, 0.0),
        "az": defaulter(az, 0.0),
        "direct": {},
    }
    if ndim == 2 and b_given["z"]:
        resolved["direct"]["z"] = defaulter(bz, 0.0)
    return resolved


def _validation_layout(sim):
    domain_box = Box([0] * sim.ndim, sim.cells)
    if sim.interp_order:
        return HybridGridLayoutFor(
            domain_box, domain_box.lower, sim.dl, sim.interp_order
        )
    ghosts = mhdGhostNbrFromReconstruction((sim.reconstruction or "").lower()) or 2
    return MHDGridLayoutFor(
        domain_box, domain_box.lower, sim.dl, sim.reconstruction, [ghosts] * sim.ndim
    )


def _plane_coords(L, R, idir, ndim):
    """coordinates along idir, zero along the other directions"""

    def coords(values):
        return tuple(values if d == idir else np.zeros_like(values) for d in range(ndim))

    return coords(L), coords(R)


def _curl(ndim, dl, ax, ay, az, b_name, coords):
    """centred-difference curl component of A at coords, half-steps on the edges"""
    A = {"x": ax, "y": ay, "z": az}

    def d(fn, idir):  # 0 along non-simulated directions
        if idir >= ndim:
            return 0.0
        h = 0.5 * dl[idir]
        plus = tuple(c + h if d == idir else c for d, c in enumerate(coords))
        minus = tuple(c - h if d == idir else c for d, c in enumerate(coords))
        return (np.asarray(fn(*plus)) - np.asarray(fn(*minus))) / dl[idir]

    if b_name == "Bx":
        return d(A["z"], 1) - d(A["y"], 2)
    if b_name == "By":
        return d(A["x"], 2) - d(A["z"], 0)
    return d(A["y"], 0) - d(A["x"], 1)


def _ghost_ranges(sim, layout, idir):
    nbrGhosts = layout.nbrGhostsPrimal(sim.interp_order)
    length = sim.simulation_domain()[idir]
    dual_left = (np.arange(-nbrGhosts, nbrGhosts) + 0.5) * sim.dl[idir]
    primal_left = np.arange(-nbrGhosts, nbrGhosts) * sim.dl[idir]
    return (dual_left, dual_left + length), (primal_left, primal_left + length)


def curl_non_periodic(sim, ax, ay, az, atol, curl_rtol=1e-10):
    """
    For each periodic direction, compare the centred-difference curl of A (step
    ``sim.dl``) on the lower and upper boundary planes, at the B face coordinates and
    over the ghost ranges the B periodicity check uses. A itself may be non-periodic
    (e.g. Harris Az, which gains a constant per period): only its curl is compared.

    The tolerance is ``atol`` plus the round-off of differencing A,
    ``eps * max|A| / dl``, since a non-periodic A can be large at the boundary, plus
    ``curl_rtol * max|B|``, so that profile tails reaching into the ghost range at
    the 1e-12 level (e.g. a tanh current sheet) are not reported.

    Returns the list of ``(B component name, direction index)`` that mismatch.
    """
    layout = _validation_layout(sim)
    not_periodic = []

    for idir in range(sim.ndim):
        if sim.boundary_types[idir] != "periodic":
            continue
        dual, primal = _ghost_ranges(sim, layout, idir)

        for b_name in ("Bx", "By", "Bz"):
            L, R = dual if layout.qtyIsDual(b_name, _directions[idir]) else primal
            coordsL, coordsR = _plane_coords(L, R, idir, sim.ndim)
            curlL = _curl(sim.ndim, sim.dl, ax, ay, az, b_name, coordsL)
            curlR = _curl(sim.ndim, sim.dl, ax, ay, az, b_name, coordsR)

            a_max = max(
                np.max(np.abs(np.asarray(fn(*c))))
                for fn in (ax, ay, az)
                for c in (coordsL, coordsR)
            )
            b_max = max(np.max(np.abs(curlL)), np.max(np.abs(curlR)))
            tol = (
                atol
                + 16 * np.finfo(float).eps * a_max / min(sim.dl)
                + curl_rtol * b_max
            )

            if not np.allclose(curlL, curlR, atol=tol, rtol=0):
                not_periodic += [(b_name, idir)]

    return not_periodic


def direct_non_periodic(sim, direct, atol):
    """
    periodicity check of the B components given directly in A mode (``{"z": bz}``),
    same as the B-mode check. Returns ``[(B component name, direction index)]``.
    """
    layout = _validation_layout(sim)
    not_periodic = []

    for idir in range(sim.ndim):
        if sim.boundary_types[idir] != "periodic":
            continue
        dual, primal = _ghost_ranges(sim, layout, idir)

        for c, b_i in direct.items():
            b_name = "B" + c
            L, R = dual if layout.qtyIsDual(b_name, _directions[idir]) else primal
            coordsL, coordsR = _plane_coords(L, R, idir, sim.ndim)
            if not np.allclose(b_i(*coordsL), b_i(*coordsR), atol=atol, rtol=0):
                not_periodic += [(b_name, idir)]

    return not_periodic


def magnetic_non_periodic(sim, model_dict, atol):
    """A-mode periodicity check: curl of A, then direct components"""
    vecpot = model_dict["vector_potential"]
    direct = model_dict["direct_b"]
    curl_bad = curl_non_periodic(
        sim, vecpot["ax"], vecpot["ay"], vecpot["az"], atol
    )
    overridden = {"B" + c for c in direct}
    return [nb for nb in curl_bad if nb[0] not in overridden] + direct_non_periodic(
        sim, direct, atol
    )
