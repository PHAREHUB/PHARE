"""Derived quantities served by ``Run`` on the new reader (DESIGN 9).

A ``Recipe`` says which files a derived quantity needs, how to compute it, and how far
around a requested box its stencil reaches:

    Recipe(name, inputs, kernel, halo, peak)

name    the quantity name given to the result ("J", "divB", "Pe", ...)
inputs  file stems read ("EM_B", "ions_charge_density")
kernel  ``kernel(ctx, **options) -> result``; ``ctx`` (a ``RecipeContext``) opens the inputs
        lazily on the read box at the requested times
halo    nodes of every level the stencil needs beyond a node: with ``box=``, the read box
        is the user's box grown by ``halo * dl_0`` (level-0 mesh size) on every side, so
        that every level gets at least ``halo`` nodes of margin; the result's extent is
        then set back to the user's box. Nodes strictly inside the user's box get the
        same values as without ``box=``.
peak    memory of the result as a multiple of the inputs' bytes (for ``Run.plan``)

The registry ``RECIPES`` holds: ``J`` (node stencil, ``curl`` of node B), ``J_yee`` (h5
only: master's Yee stencil moved to nodes), ``divB`` (node formula; on h5 the Yee divB
averaged to nodes, equal to the node formula away from coarse-fine faces), ``Pe`` (Ni * Te,
h5 only), ``Az`` (magnetic flux function, 2D, from level 0).

Derived quantities are computed when the getter is called (eager, in memory): pass ``box=``
to compute them on a region only. Their inputs are opened through ``Run``, which reuses a
whole-domain input field the caller still holds (e.g. ``GetJ`` after ``GetB(t).load()`` reads
nothing).
"""

from __future__ import annotations

import math
from collections.abc import Callable
from dataclasses import dataclass

import numpy as np

from .boxes import PhysBox
from .store import LevelData, Snapshot

VECTOR_NAMES = ("x", "y", "z")


@dataclass(frozen=True)
class Recipe:
    """A derived quantity (see the module docstring for the fields)."""

    name: str
    inputs: tuple
    kernel: Callable
    halo: int
    peak: float
    doc: str = ""


class RecipeContext:
    """What a recipe kernel may use.

    opener    : callable ``opener(stem, times, extent, **kw) -> AMRField`` (lazy, no field
                bytes read); ``kw`` are ``AMRField.from_backend`` keywords (derived, names,
                dtype, quantity, stencil)
    times     : the time(s) asked for (number or list of numbers)
    read_box  : PhysBox the inputs are opened on (the user's box widened), or None
    formats   : ``{stem: "h5" | "vtkhdf"}`` of the input files
    sim       : callable returning the simulation object, or None when the files have none
    """

    def __init__(self, opener, times, read_box, formats, sim=None):
        self.opener = opener
        self.times = times
        self.read_box = read_box
        self.formats = dict(formats)
        self._sim = sim

    def field(self, stem, **kw):
        """Input ``stem`` opened lazily at ``self.times`` on ``self.read_box``."""
        return self.opener(stem, self.times, self.read_box, **kw)

    def format(self, stem):
        """"h5" or "vtkhdf": the format of the input file ``stem``."""
        return self.formats[stem]

    def sim(self):
        """The simulation object of the run (ValueError when the files do not have one)."""
        if self._sim is None:
            raise ValueError("the files of this run carry no simulation parameters")
        return self._sim()


# --------------------------------------------------------------------------- helpers


def widen(box, halo, dl0):
    """``box`` (PhysBox) grown by ``halo * dl0[a]`` on both sides of every bounded axis.

    Unbounded sides stay unbounded; the result may stick out of the domain (selection
    clips it). None stays None.
    """
    if box is None or halo == 0:
        return box
    lo = tuple(v - halo * float(d) if math.isfinite(v) else v for v, d in zip(box.lo, dl0))
    hi = tuple(v + halo * float(d) if math.isfinite(v) else v for v, d in zip(box.hi, dl0))
    return PhysBox(lo, hi)


def relabel(field, quantity, stencil=None, extent=False):
    """Same data with another quantity name / stencil, and optionally another extent.

    extent : False keeps the extent; a PhysBox or None replaces it (shares every buffer).
    """
    snaps = field.snapshots()
    if extent is not False:
        snaps = [Snapshot(s.hier.with_extent(extent), s.levels) for s in snaps]
    kw = dict(quantity=quantity)
    if stencil is not None:
        kw["stencil"] = stencil
    return field._like(snaps, **kw)


def _per_time(field, func):
    """``[func(snapshot) for every time]``."""
    return [func(s) for s in field.snapshots()]


# DESIGN-NOTE: the patches selected by the widened box are kept in the result (only the extent
# is set back to the user's box): they hold valid values and integrals/finest/plot follow the
# extent. Periodic wrap is not widened across the domain edge (PhysBox clips): nodes of a box
# touching x=0 on a periodic axis may use one-sided derivatives there.
def evaluate(recipe, opener, times, box, dl0, formats, sim=None, **options):
    """Run ``recipe`` at ``times`` on ``box`` (PhysBox or None) and return its result.

    The inputs are opened on ``widen(box, recipe.halo, dl0)``; a field result gets the
    user's ``box`` back as its extent. Raises what the kernel raises.
    """
    read_box = widen(box, recipe.halo, dl0)
    ctx = RecipeContext(opener, times, read_box, formats, sim)
    out = recipe.kernel(ctx, **options)
    if box is not None and hasattr(out, "snapshots"):
        out = relabel(out, out.quantity, extent=box)
    return out


# --------------------------------------------------------------------------- kernels


def _current_node(ctx):
    from .ops import deriv

    B = ctx.field("EM_B")
    return B._like(_per_time(B, deriv.curl), quantity="J", stencil="node")


def _current_yee(ctx):
    if ctx.format("EM_B") != "h5":
        raise ValueError(
            "GetJ(stencil='yee') needs the Yee-centred B of a phareh5 file (EM_B.h5); a "
            "vtkhdf file holds node values only: use stencil='node' (the default)"
        )
    return ctx.field("EM_B", derived="J_yee", names=VECTOR_NAMES, quantity="J",
                     stencil="yee")


def _mask_one_sided(field):
    """In-memory copy of a scalar field with NaN where a one-sided derivative was used."""
    from .ops.deriv import one_sided_mask

    snaps = []
    for s in field.snapshots():
        levels = []
        for ld in s.levels:
            v = np.array(ld.buffer("value"), dtype=ld.dtype, copy=True)
            v[one_sided_mask(ld.geom)] = np.nan
            levels.append(LevelData.from_arrays(ld.geom, ld.time, {"value": v}))
        snaps.append(Snapshot(s.hier, tuple(levels)))
    return field._like(snaps)


def default_mask_faces(fmt):
    """``GetDivB``'s ``mask_faces`` default for an input file format ("h5" / "vtkhdf").

    vtkhdf: True. A node file has no data across a coarse-fine face, so divB there comes from
    one-sided derivatives and is dominated by their truncation error (76-119x the h5 value
    on the sample runs, all of it at those nodes): the file cannot give divB there, NaN says
    so. h5: False. The Yee divB comes from the file's ghost values and is meaningful there.
    """
    return fmt != "h5"


def _divergence(ctx, mask_faces=None):
    from .ops import deriv

    if mask_faces is None:
        mask_faces = default_mask_faces(ctx.format("EM_B"))
    if ctx.format("EM_B") == "h5":
        div = ctx.field("EM_B", derived="divB_yee", names=("value",), quantity="divB",
                        stencil="yee")
        return _mask_one_sided(div) if mask_faces else div
    B = ctx.field("EM_B")
    snaps = _per_time(B, lambda s: deriv.div(s, mask_faces=mask_faces))
    return B._like(snaps, quantity="divB", stencil="node")


def electron_temperature(sim):
    """``Te`` of an isothermal closure, read as master does (run.py:167)."""
    return sim.electrons.closure.Te


def _electron_pressure(ctx):
    if ctx.format("ions_charge_density") != "h5":
        raise ValueError(
            "GetPe needs the electron temperature, which is stored in the simulation "
            "parameters of phareh5 files; a vtkhdf file has no simulation parameters. Compute "
            "it yourself: run.GetNi(t) * Te"
        )
    Te = electron_temperature(ctx.sim())
    Ni = ctx.field("ions_charge_density")
    return relabel(Ni * Te, "Pe")


def magnetic_flux(B, interp="nearest", xn=None, yn=None, Xn=None, Yn=None, time=None):
    """Magnetic flux function ``Az`` of a 2D B field, master's algorithm (run.py:248-274).

    B      : VectorField (one time, or ``time=``), read at level 0 only (``finest(max_level=0)``
             paints level 0 on its node grid; finer levels are never read)
    xn, yn : 1D node coordinates of the output grid. None -> the level-0 nodes
             ``0, dl0, ..., L`` (master used the finest level's spacing; see DESIGN-NOTE).
             Given -> the level-0 values are interpolated there (``interp="nearest"`` or
             "linear"/"bilinear") on ``Xn, Yn`` (default ``meshgrid(xn, yn, indexing="ij")``).
    Returns ``(Az, (xn, yn))`` with ``Az[i, j]`` at ``(xn[i], yn[j])``, float64:
    ``Az = -int_0^x By(x', 0) dx' + int_0^y Bx(x, y') dy'`` (cumulative trapezoids).
    Raises ValueError for a field that is not 2D.
    """
    # DESIGN-NOTE: DESIGN 9 feeds this recipe with finest(max_level=0); master integrated the
    # merged finest data on the finest level's grid. The default grid here is level 0's.
    from scipy.integrate import cumulative_trapezoid

    if B.ndim != 2:
        raise ValueError(f"the magnetic flux function needs a 2D field, this one is {B.ndim}D")
    g = B.finest(time=time, max_level=0, names=("x", "y"))
    bx, by = (np.asarray(g[c], dtype=np.float64) for c in ("x", "y"))
    if xn is None or yn is None:
        xn, yn = np.asarray(g.x, dtype=np.float64), np.asarray(g.y, dtype=np.float64)
    else:
        from scipy.interpolate import RegularGridInterpolator

        xn, yn = np.asarray(xn, dtype=np.float64), np.asarray(yn, dtype=np.float64)
        if Xn is None or Yn is None:
            Xn, Yn = np.meshgrid(xn, yn, indexing="ij")
        method = "nearest" if interp == "nearest" else "linear"
        pts = np.stack([np.asarray(Xn).ravel(), np.asarray(Yn).ravel()], axis=-1)

        def at(v):
            f = RegularGridInterpolator((g.x, g.y), v, method=method, bounds_error=False,
                                        fill_value=None)
            return f(pts).reshape(np.shape(Xn))

        bx, by = at(bx), at(by)
    Az_x0 = -cumulative_trapezoid(by[:, 0], xn, initial=0)
    Az = cumulative_trapezoid(bx, yn, axis=1, initial=0)
    Az += Az_x0[:, np.newaxis]
    return Az, (xn, yn)


def _flux_kernel(ctx, interp="nearest", xn=None, yn=None, Xn=None, Yn=None):
    return magnetic_flux(ctx.field("EM_B"), interp, xn, yn, Xn, Yn)


RECIPES = {
    r.name: r
    for r in (
        Recipe("J", ("EM_B",), _current_node, halo=1, peak=1.0,
               doc="curl of node B: PHARE's Yee J with a (1,2,1)/4 smoothing per derivative; "
                   "its peak is 7-50 % below the Yee J on under-resolved sheets"),
        Recipe("J_yee", ("EM_B",), _current_yee, halo=0, peak=1.0,
               doc="master's Yee J per patch (ghosts), averaged to nodes; phareh5 only"),
        Recipe("divB", ("EM_B",), _divergence, halo=1, peak=1.0 / 3.0,
               doc="sum of centred derivatives (h5: Yee divB averaged to nodes, the same "
                   "away from coarse-fine faces); vtkhdf: NaN at one-sided nodes by default"),
        Recipe("Pe", ("ions_charge_density",), _electron_pressure, halo=0, peak=1.0,
               doc="Ni * Te (isothermal closure); phareh5 only"),
        Recipe("Az", ("EM_B",), _flux_kernel, halo=0, peak=1.0,
               doc="magnetic flux function from level 0, 2D"),
    )
}


__all__ = ["RECIPES", "Recipe", "RecipeContext", "evaluate", "magnetic_flux", "relabel",
           "widen"]
