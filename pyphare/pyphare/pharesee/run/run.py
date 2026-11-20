import functools
import glob
import inspect
import numbers
import os
import warnings
import weakref
from typing import NamedTuple

import numpy as np

from pyphare.core.phare_utilities import listify, refinement_ratio
from pyphare.logger import getLogger
from pyphare.pharesee.hierarchy import ScalarField, VectorField, hierarchy_from
from pyphare.pharesee.hierarchy.fromh5 import pop_attribute_from_h5
from pyphare.pharesee.hierarchy.hierarchy_utils import (
    compute_hier_from,
    flat_finest_field,
)

from .utils import (
    _compute_current,
    _compute_divB,
    _compute_pop_pressure,
    _compute_pressure,
    _compute_to_primal,
    _get_rank,
    make_interpolator,
)

logger = getLogger(__name__)

quantities_per_file = {
    "EM_B": "B",
    "EM_E": "E",
    "ions_bulkVelocity": "Vi",
    "ions_charge_flux": "CF",
    "ions_charge_density": "Ni",
    "ions_mass_density": "rho",
    "particle_count": "nppc",
    "mhd_rho": "mhd_rho",
    "mhd_V": "mhd_V",
    "mhd_P": "mhd_P",
    "mhd_rhoV": "mhd_rhoV",
    "mhd_Etot": "mhd_Etot",
}


class LegacyRun:
    """Master's reader of phareh5 runs, unchanged (``Run`` routes to it on the legacy path)."""

    def __init__(self, path, default_time=None):
        self.path = path
        self.default_time_ = default_time
        self.available_diags = glob.glob(os.path.join(self.path, "*.h5"))

    def _get_hierarchy(self, times, filename, hier=None, **kwargs):
        from pyphare.core.box import Box

        times = listify(times)
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

    # TODO maybe transform that so multiple times can be accepted
    def _get(self, hierarchy, time, merged, interp):
        """
        if merged=True, will return an interpolator and a tuple of 1d arrays
        with the coordinates of the finest grid where the interpolator
        can be calculated (that is the return of flat_finest_field)
        """
        if merged:
            domain = self.GetDomainSize()
            dl = self.GetDl(hierarchy.finest_level(time))

            # assumes all qties in the hierarchy have the same ghost width
            # so take the first patch data of the first patch of the first level....
            nbrGhosts = next(
                iter(hierarchy.level(0).patches[0].patch_datas.values())
            ).ghosts_nbr
            merged_qties = {}
            for qty in hierarchy.quantities():
                data, coords = flat_finest_field(hierarchy, qty, time=time)
                merged_qties[qty] = make_interpolator(
                    data, coords, interp, domain, dl, qty, nbrGhosts
                )
            return merged_qties
        else:
            return hierarchy

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
        return ScalarField.FROM(self._get(hier, time, merged, interp))

    def GetN(self, time, pop_name, merged=False, interp="nearest", **kwargs):
        hier = self._get_hierarchy(time, f"ions_pop_{pop_name}_density.h5", **kwargs)
        return ScalarField.FROM(self._get(hier, time, merged, interp))

    def GetVi(self, time, merged=False, interp="nearest", **kwargs):
        hier = self._get_hierarchy(time, "ions_bulkVelocity.h5", **kwargs)
        return VectorField.FROM(self._get(hier, time, merged, interp))

    def GetChargeFlux(self, time, merged=False, interp="nearest", **kwargs):
        hier = self._get_hierarchy(time, "ions_charge_flux.h5", **kwargs)
        return VectorField.FROM(self._get(hier, time, merged, interp))

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

        h = compute_hier_from(_compute_to_primal, hier, scalar="rho")
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
        hier = self._get_hierarchy(time, "EM_B.h5")

        # Reuse grids if provided, otherwise generate them
        if xn is None or yn is None or Xn is None or Yn is None:
            domain = self.GetDomainSize()
            dl = self.GetDl(hier.finest_level(time))
            xn = np.arange(0, domain[0] + dl[0], dl[0])
            yn = np.arange(0, domain[1] + dl[1], dl[1])
            Xn, Yn = np.meshgrid(xn, yn, indexing="ij")

        merged_B = self._get(hier, time, merged=True, interp=interp)
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
        hier = self._get_hierarchy(times[0], "EM_B.h5")
        dl = self.GetDl(hier.finest_level(times[0]))
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
        return pop_attribute_from_h5(self.path, pop_name, "pop_mass")

    def GetCharge(self, pop_name, **kwargs):
        return pop_attribute_from_h5(self.path, pop_name, "pop_charge")

    def GetDomainSize(self, **kwargs):
        import h5py

        data_file = h5py.File(
            self.available_diags[0], "r"
        )  # That is the first file in th available diags
        root_cell_width = np.asarray(data_file.attrs["cell_width"])

        return (data_file.attrs["domain_box"] + 1) * root_cell_width

    def GetDl(self, level=0):
        """
        gives the ndarray containing the grid sizes for the given level

        :param level: the level at which get the associated grid size
        """
        import h5py

        if level < 0:
            raise ValueError(f"Run::GetDl: level must be >= 0, got {level}")

        with h5py.File(self.available_diags[0], "r") as data_file:
            root_cell_width = np.asarray(data_file.attrs["cell_width"])
        return root_cell_width / refinement_ratio**level

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
            ts[quantities_per_file.get(basename, basename)] = time
            ff.close()
        return ts

    def times(self, qty):
        return self.all_times()[qty]


# ============================================================================ routing
#
# ``Run`` decides ONCE, in ``__init__``, which reader serves the files of a directory:
#   * ``.vtkhdf`` files: always the new reader (``pyphare.pharesee.amr``); the legacy
#     reader cannot read them.
#   * ``.h5`` files: the new reader by default. What it does not provide is served by the
#     legacy reader (``LegacyRun``, master's code, unchanged) automatically: ``merged=True``,
#     particles, pressure tensors, tags, ranks, particle counts, ``all_primal=False``.
#     ``Run(path, legacy=True)`` or ``PHARESEE_LEGACY=1`` sends every h5 file to master's
#     reader (escape hatch).
# When a directory holds only legacy-served files, every method is master's, byte for byte.

_LEGACY_ENV = "PHARESEE_LEGACY"
_BOX_KEYWORDS = ("box", "selection_box")

_MERGED_MSG = (
    "{getter}(merged=True) is not available on the new reader: use field.finest(box=...) "
    "for a regular grid of the finest data, or probes field(x=..., y=...) for values at "
    "points"
)


class _Route(NamedTuple):
    legacy: bool  # True: master's code serves the quantity
    path: str  # the file, or None on an all-legacy run (master finds it itself)
    fmt: str  # "h5" or "vtkhdf"


# fields of the new h5 path carry a lazy proxy as ``sim`` (deserialised on first
# attribute access) so that GetB stays at index cost; ``field.sim is None`` only on vtkhdf.
class _LazySim:
    """The simulation object of a phareh5 file, deserialised on first attribute access.

    Master reads it from ``/py_attrs`` (hierarchy.py:104-126); doing so costs ~40 ms, so the
    new reader defers it until a getter (GetPe) or the user needs it.
    """

    __slots__ = ("_path", "_sim")

    def __init__(self, path):
        self._path = path
        self._sim = None

    def load(self):
        """The ``pharein.Simulation`` (deserialised once). ValueError if the file has none."""
        if self._sim is None:
            import h5py

            from pyphare.pharein.simulation import deserialize

            with h5py.File(self._path, "r") as f:
                if "py_attrs" not in f:
                    raise ValueError(
                        f"{self._path}: Simulation is not available for deserialization"
                    )
                text = f["py_attrs"].attrs["serialized_simulation"]
            self._sim = deserialize(text)
        return self._sim

    def __getattr__(self, name):
        if name.startswith("__"):
            raise AttributeError(name)
        return getattr(self.load(), name)

    def __reduce__(self):
        # the loaded Simulation holds the user's (often local) functions: pickle the path only
        return (_LazySim, (self._path,))

    def __repr__(self):
        state = "loaded" if self._sim is not None else "not loaded yet"
        return f"<simulation of {self._path} ({state})>"


class RunPlan:
    """What a read of a quantity would touch, computed without reading field bytes.

    nbytes    bytes read from the file (all levels and times); not the memory the
              computation will use, which is larger
    patches   ``{(time, level): number of patches read}`` (patches already in memory, e.g.
              of a field the caller holds, are not counted)
    """

    def __init__(self, what, nbytes, patches):
        self.what = what
        self.nbytes = int(nbytes)
        self.patches = dict(patches)

    def __repr__(self):
        from pyphare.pharesee.amr.plan import format_bytes

        npatch = sum(self.patches.values())
        return (
            f"RunPlan({self.what}: {format_bytes(self.nbytes)} read from the file, "
            f"in {npatch} patches)"
        )


@functools.cache
def _master_params(getter):
    """Names of the positional parameters of master's ``getter`` (LegacyRun), without self."""
    params = list(inspect.signature(getattr(LegacyRun, getter)).parameters.values())[1:]
    return tuple(p.name for p in params if p.kind is inspect.Parameter.POSITIONAL_OR_KEYWORD)


def _master_arguments(getter, args):
    """``{name: value}`` of master's positional parameters of ``getter`` given in ``args``;
    TypeError if ``args`` has more values than ``getter`` has parameters."""
    params = _master_params(getter)
    if len(args) > len(params):
        raise TypeError(f"master's {getter}() takes {len(params)} positional arguments "
                        f"({', '.join(params)}), {len(args)} given")
    return dict(zip(params, args))


def _check_keywords(getter, kwargs, accepted, signature):
    """TypeError naming every keyword of ``kwargs`` not in ``accepted``."""
    unknown = [k for k in kwargs if k not in accepted]
    if unknown:
        allowed = ", ".join(list(signature) + list(accepted))
        raise TypeError(
            f"{getter}() got unexpected keyword argument(s) {unknown}; on the new reader it "
            f"accepts: {allowed}"
        )


def _to_legacy(route, merged):
    """True when master's reader serves this call: the file is routed to it, or
    ``merged=True`` on a phareh5 file (the new reader has no merged interpolators)."""
    return route.legacy or (merged and route.fmt == "h5")


def _check_merged_interp(getter, merged, interp, h5):
    if merged:
        # only vtkhdf gets here: on h5, merged=True goes to master's reader (_to_legacy)
        raise ValueError(_MERGED_MSG.format(getter=getter))
    if interp != "nearest":
        raise TypeError(
            f"{getter}(interp={interp!r}): interp only applies to merged=True, which the new "
            "reader does not have; drop it (finest(fill='cubic'|'nearest') and probes choose "
            "their own interpolation)"
            + ("; on phareh5 files merged=True is served by master's reader" if h5 else "")
        )


# a dict selection_box is not master's (master took a Box or a 2D tuple): it is
# read as box= (PhysBox rule) with the same deprecation warning.
def _selection_to_box(selection_box, ndim):
    """master's ``selection_box`` -> ``(extent PhysBox, MasterSelection or None)``.

    Box or ``(lower..., upper...)`` tuple: the patches are those master selects
    (fromh5.py:235-247, ``field.MasterSelection``): the bounds are truncated to integers
    as master's ``Box`` does (core/box.py:15-16, issue #1451) and a patch whose truncated
    physical box touches the selection is kept. The field's extent is the box as given
    (floats kept). A dict: the PhysBox rule of ``box=``.
    """
    from pyphare.core.box import Box
    from pyphare.pharesee.amr.boxes import PhysBox
    from pyphare.pharesee.amr.field import MasterSelection

    if isinstance(selection_box, Box):
        lo, hi = np.asarray(selection_box.lower), np.asarray(selection_box.upper)
    elif isinstance(selection_box, dict):
        return PhysBox.parse(selection_box, ndim), None
    elif isinstance(selection_box, (tuple, list)) and len(selection_box) == 2 * ndim and all(
        isinstance(v, numbers.Number) for v in selection_box
    ):
        lo, hi = selection_box[:ndim], selection_box[ndim:]
    else:
        raise TypeError(
            f"selection_box={selection_box!r}: expected a Box or a tuple (lower..., upper...) "
            f"of {2 * ndim} numbers; better: box=dict(x=(x0, x1), ...)"
        )
    if len(lo) != ndim:
        raise ValueError(f"selection_box has {len(lo)} axes, the data has {ndim}")
    extent = PhysBox(tuple(float(v) for v in lo), tuple(float(v) for v in hi))
    # master: Box(lower, upper) keeps lower.astype(int), upper.astype(int)
    selection = MasterSelection(np.asarray(lo, dtype=np.float64).astype(int),
                                np.asarray(hi, dtype=np.float64).astype(int))
    return extent, selection


def _parse_legacy_env(value):
    """``PHARESEE_LEGACY`` -> True (legacy), False (new reader) or None (unset/empty).

    Accepts 0/1/true/false/yes/no in any case; ValueError naming the variable otherwise.
    """
    text = (value or "").strip().lower()
    if not text:
        return None
    if text in ("1", "true", "yes"):
        return True
    if text in ("0", "false", "no"):
        return False
    raise ValueError(
        f"{_LEGACY_ENV}={value!r} is not understood: use 1/true/yes (master's reader for .h5 "
        f"files) or 0/false/no (the new reader), or unset it"
    )


def _times_argument(time, default_time, getter):
    if time is None:
        if default_time is None:
            raise TypeError(f"{getter}: give a time (run.times(qty) lists them)")
        return default_time
    if not isinstance(time, (str, bytes)) and np.ndim(time) > 0 and np.size(time) == 0:
        raise ValueError(f"{getter}: no time given (an empty list); give a time or a list "
                         "of times (run.times(qty) lists them)")
    return time


class Run(LegacyRun):
    """Entry point to the output of a PHARE run (one directory of diagnostic files).

    ``Run(path, default_time=None, *, legacy=None, format=None)``

    path
        the diagnostics directory
    legacy
        None: the environment variable ``PHARESEE_LEGACY`` if set, else False;
        False: the new reader (``pyphare.pharesee.amr``) reads ``.h5`` files, and
        what it does not provide (``merged=True``, particles, pressure tensors, tags,
        ranks, particle counts, ``all_primal=False``) is served by master's reader
        (``LegacyRun``) automatically; True: master's reader for every ``.h5`` file
        (raises ValueError when the selected files are ``.vtkhdf``, which it cannot read)
    format
        None, "h5" or "vtkhdf": which files to use when a quantity exists in both
        formats (without it such a quantity raises ValueError asking for it)

    The reader is chosen once, here; ``repr(run)`` says which. ``.vtkhdf`` files always go
    to the new reader.

    New reader: every getter keeps master's signature and also takes the keywords
    ``box=dict(x=(x0, x1), ...)`` (physical units; only the patches meeting the box are
    read, the box becomes the field's extent) and ``selection_box=`` (master's meaning,
    deprecated). Fields are lazy (bytes are read at first use), node centred, without
    ghosts, in the file's dtype. ``time`` is a number or a list of numbers. ``merged=True``
    goes to master's reader on h5 and raises on vtkhdf; unknown keywords raise; see each
    getter for the rest.
    ``selection_box=`` selects exactly the patches master selects (integer bounds, touching
    patches kept); ``box=`` uses the physical rule of ``PhysBox``.
    Reuse: the Run keeps WEAK references to the whole-domain fields it returned, per (file,
    time). Asking again for the same quantity and time while the first field is alive
    shares its data (no second read); derived quantities (GetJ, GetDivB on vtkhdf, GetPe,
    GetMagneticFlux) use such a held B or Ni; a ``box=`` request whose whole-domain field is
    alive and fully loaded is cut from it in memory. Dropping the field frees the memory.
    Fields sharing data share loading too: ``unload()`` on one unloads the others.
    """

    def __init__(self, path, default_time=None, *, legacy=None, format=None):
        from pyphare.pharesee.amr.io.files import normalize_format, scan

        if legacy is not None and not isinstance(legacy, bool):
            raise TypeError(
                f"Run(legacy={legacy!r}): legacy must be True, False or None (the default: "
                f"the environment variable {_LEGACY_ENV}, else the new reader)"
            )
        super().__init__(path, default_time)
        self._format = normalize_format(format)
        found = scan(path) if os.path.isdir(path) else {}
        self._has = {
            fmt: sorted(f[fmt] for f in found.values() if fmt in f) for fmt in ("h5", "vtkhdf")
        }
        if legacy is None:
            env = _parse_legacy_env(os.environ.get(_LEGACY_ENV))
            h5_legacy = False if env is None else env
        else:
            h5_legacy = legacy
        only_vtk = self._format == "vtkhdf" or (
            self._format is None and self._has["vtkhdf"] and not self._has["h5"]
        )
        if legacy is True and only_vtk:
            raise ValueError(
                f"{path}: the legacy reader cannot read .vtkhdf files; drop legacy=True "
                "(vtkhdf files always use the new reader), or point to phareh5 output"
            )
        self._h5_legacy = h5_legacy
        # every quantity served by master's code: nothing of the new reader runs
        self._all_legacy = h5_legacy and (self._format == "h5" or not self._has["vtkhdf"])
        self._backends = {}
        self._backend_sigs = {}  # path -> (inode, mtime_ns, size) when its backend was opened
        self._sims = {}
        # created when a file is first routed to the new reader: the legacy path never reads
        # the new reader's configuration (PHARESEE_GEOM_CACHE_MB ...)
        self._geom_cache = None
        # (file path, time) -> weak references to whole-domain fields this Run returned
        self._held = {}

    # ------------------------------------------------------------------ description
    @property
    def reader(self):
        """"legacy" (master's code for every file), "new", or "mixed"."""
        if self._all_legacy:
            return "legacy"
        uses_h5 = self._format != "vtkhdf" and self._has["h5"]
        if self._h5_legacy and uses_h5:
            return "mixed"
        return "new"

    def __repr__(self):
        fmts = ", ".join(f"{fmt} ({len(files)} files)" for fmt, files in self._has.items()
                         if files) or "no diagnostic file yet"
        h5 = "legacy reader" if self._h5_legacy else "new reader"
        if self._all_legacy:
            how = "legacy reader (master's code) for every file"
        else:
            parts = []
            if self._has["vtkhdf"] and self._format != "h5":
                parts.append("vtkhdf: new reader")
            if self._has["h5"] and self._format != "vtkhdf":
                parts.append(f"h5: {h5}")
            how = ", ".join(parts) or "new reader"
        chosen = f", format={self._format!r}" if self._format else ""
        return f"Run({self.path!r}{chosen}): {fmts}; {how}"

    # ------------------------------------------------------------------ plumbing
    def _route(self, stem):
        if self._all_legacy:
            return _Route(True, None, "h5")
        from pyphare.pharesee.amr.io.files import find_quantity_file, scan

        fmt_wanted = self._format
        # legacy=True (or PHARESEE_LEGACY=1) asked for master's reader, which reads h5 only:
        # a stem present in both formats is read from its h5 file
        if fmt_wanted is None and self._h5_legacy and len(scan(self.path).get(stem, {})) > 1:
            fmt_wanted = "h5"
        path = find_quantity_file(self.path, stem, fmt_wanted)
        fmt = "h5" if path.endswith(".h5") else "vtkhdf"
        return _Route(fmt == "h5" and self._h5_legacy, path, fmt)

    # master's getters call other getters on ``self`` (GetPressure calls GetFlux,
    # GetN, GetMass; GetPi calls GetMassDensity; GetJ/GetDivB/GetRanks call GetB). Called on
    # this Run they would get the new reader's node-centred fields; a plain LegacyRun with
    # the same state keeps every call inside master's code.
    def _master(self):
        """A LegacyRun with this Run's state: master's getters, all on master's reader."""
        master = LegacyRun.__new__(LegacyRun)
        master.path = self.path
        master.default_time_ = self.default_time_
        master.available_diags = self.available_diags
        return master

    def _backend(self, path):
        """The backend of ``path``, reopened when the file changed (inode, size or mtime)
        since it was opened (a run still writing appends dumps: the new times become
        visible; a replaced file is read anew)."""
        st = os.stat(path)
        signature = (st.st_ino, st.st_mtime_ns, st.st_size)
        backend = self._backends.get(path)
        if backend is None or self._backend_sigs.get(path) != signature:
            from pyphare.pharesee.amr.io import open_backend

            if backend is not None:
                self._forget_held_if_rewritten(path, backend)
                # the simulation is read again from the file as it is now; fields already
                # returned keep theirs
                self._sims.pop(path, None)
            backend = self._backends[path] = open_backend(path)
            self._backend_sigs[path] = signature
        return backend

    def _forget_held_if_rewritten(self, path, old):
        """Drop the held fields of ``path`` (``_reuse_held``) unless ``old``, the backend
        that read them, finds the times it served unchanged (an append keeps them; a
        replaced, shrunk or rewritten file does not: its data are another run's)."""
        from pyphare.pharesee.amr.errors import PhareseeFileChanged

        held = [k for k in self._held if k[0] == path]
        if not held:
            return
        try:
            old.index(held[0][1])  # checks every time this backend indexed (_revalidate)
        except PhareseeFileChanged:
            for key in held:
                del self._held[key]

    def _sim(self, path):
        sim = self._sims.get(path)
        if sim is None:
            sim = self._sims[path] = _LazySim(path)
        return sim

    def _geometry_cache(self):
        """The GeometryCache of this Run, created on first use (new reader only)."""
        if self._geom_cache is None:
            from pyphare.pharesee.amr.geometry import GeometryCache

            self._geom_cache = GeometryCache()
        return self._geom_cache

    def _open_path(self, path, times, extent, selection=None, **kw):
        """Lazy AMRField of the file at ``path`` (index tables only, no field bytes).

        Reuses what this Run already returned (``_reuse_held``) for a plain stored field
        (no ``derived``/``names``/``dtype``/``quantity``/``stencil`` keyword).
        """
        from pyphare.pharesee.amr.field import AMRField

        backend = self._backend(path)
        if backend.meta().format == "phareh5":
            kw.setdefault("sim", self._sim(path))
        plain = set(kw) <= {"sim"}
        field = AMRField.from_backend(backend, times=times, extent=extent,
                                      cache=self._geometry_cache(), selection=selection, **kw)
        if plain:
            field = self._reuse_held(path, field, whole=extent is None and selection is None)
        return field

    # only weak references are kept, so dropping the field frees its memory. A
    # boxed field cut from a held whole-domain field is in memory (copies of its patches);
    # pickling it still reopens the file (it keeps the file origin).
    def _reuse_held(self, path, field, whole):
        """``field`` sharing the data of the whole-domain fields of ``path`` still alive.

        whole=True : the snapshots of a held field at the same times replace the new ones
                     (same LevelData: what one read, the other has); ``field`` is then held
                     (weakly) for the next calls.
        whole=False: a time whose held whole-domain field is fully loaded is cut from it in
                     memory (``ops.arith.restrict``), so the file is not read again; other
                     times stay lazy on the file.
        Reads nothing.
        """
        snaps, changed = list(field.snapshots()), False
        for k, snap in enumerate(snaps):
            held = self._held_snapshot(path, snap.time)
            if held is None:
                continue
            if whole:
                snaps[k], changed = held, True
            elif all(bool(np.all(ld.loaded)) for ld in held.levels):
                from pyphare.pharesee.amr.ops.arith import restrict

                snaps[k], changed = restrict(held, snap.hier), True
        if changed:
            field = field._sharing(snaps)
        if whole:
            for t in field.times():
                key = (path, float(t))
                refs = [r for r in self._held.get(key, ()) if r() is not None]
                refs.append(weakref.ref(field))
                self._held[key] = refs
        return field

    def _held_snapshot(self, path, time):
        """Snapshot at ``time`` of a live whole-domain field of ``path``, or None."""
        key = (path, float(time))
        refs = self._held.get(key)
        if not refs:
            return None
        for ref in refs:
            f = ref()
            if f is not None:
                return f.snapshot(float(time))
        del self._held[key]  # all dead
        return None

    def _extent(self, getter, opts, ndim, stacklevel=4):
        """``(extent PhysBox or None, MasterSelection or None)`` from box= / selection_box=.

        stacklevel: of the DeprecationWarning, so that it points at the user's call (4: a
        getter calls ``_stored`` or ``_derived_getter``, which call this).
        """
        from pyphare.pharesee.amr.boxes import PhysBox

        box, sb = opts.get("box"), opts.get("selection_box")
        if box is not None and sb is not None:
            raise TypeError(f"{getter}: give box= or selection_box=, not both")
        if sb is not None:
            warnings.warn(
                f"{getter}(selection_box=...) is deprecated: use box=dict(x=(x0, x1), "
                "y=(y0, y1)) (physical units, float bounds kept as given; selection_box "
                "truncates them to integers and keeps patches touching the box, as master)",
                DeprecationWarning,
                stacklevel=stacklevel,
            )
            return _selection_to_box(sb, ndim)
        return PhysBox.parse(box, ndim), None

    # Every getter below keeps master's signature and hands its positional arguments, in
    # master's order, to one of three paths: _stored (a quantity in a file), _derived_getter
    # (computed from one), _unavailable (not on the new reader). Each path either serves the
    # call on the new reader or gives it to master's code through _legacy.
    def _legacy(self, getter, args, kwargs):
        """Master's ``getter`` (LegacyRun) called with ``args``, ``kwargs`` on master's reader."""
        if "box" in kwargs:
            # master's hierarchy_from has no box=: it would read the whole domain
            raise TypeError(
                f"{getter}(box=...) is a keyword of the new reader; this call is served by "
                "master's reader (merged=True, all_primal=False, legacy=True or a quantity the "
                "new reader does not provide), which takes selection_box= instead"
            )
        # by name, bound to master's own signature: the count is checked there, and a
        # static analyser does not guess it for every getter this one line may call
        return getattr(LegacyRun, getter)(self._master(), **_master_arguments(getter, args),
                                          **kwargs)

    def _new_route(self, getter, stem, args, kwargs, extra=(), check_interp=True):
        """Route of ``stem`` when the new reader serves this call, None when master's does.

        args: master's positional arguments of ``getter``, in order. Raises on keywords the
        new reader does not take (``box``, ``selection_box`` and ``extra`` are accepted) and,
        with ``check_interp``, on merged/interp, which the new reader does not have.
        """
        a = _master_arguments(getter, args)
        route = self._route(stem)
        if _to_legacy(route, a.get("merged", False)):
            return None
        _check_keywords(getter, kwargs, _BOX_KEYWORDS + tuple(extra), _master_params(getter))
        if check_interp:
            _check_merged_interp(getter, a.get("merged", False), a.get("interp", "nearest"),
                                 route.fmt == "h5")
        return route

    def _stored(self, getter, stem, args, kwargs):
        """A quantity stored in a file: legacy object or lazy AMRField (see the getters)."""
        route = self._new_route(getter, stem, args, kwargs)
        if route is None:
            return self._legacy(getter, args, kwargs)
        a = _master_arguments(getter, args)
        if not a.get("all_primal", True):
            return self._yee_on_legacy(getter, route, args, kwargs)
        extent, selection = self._extent(getter, kwargs, self._backend(route.path).meta().dim)
        return self._open_path(route.path, _times_argument(a["time"], self.default_time_, getter),
                               extent, selection)

    def _derived_getter(self, getter, recipe, stem, args, kwargs, **options):
        """A quantity computed by ``recipe`` (amr.recipes) on the new reader, from ``stem``."""
        a = _master_arguments(getter, args)
        route = self._route(stem)
        extent, selection = self._extent(getter, kwargs, self._backend(route.path).meta().dim)
        return self._derived(recipe, _times_argument(a["time"], self.default_time_, getter),
                             extent, selection, **options)

    def _yee_on_legacy(self, getter, route, args, kwargs):
        """``all_primal=False``: master's Yee object on h5, an error on vtkhdf."""
        if route.fmt == "vtkhdf":
            raise ValueError(
                f"{getter}(all_primal=False): a vtkhdf file holds node values only, there is "
                "no Yee-centred data to return; drop all_primal=False"
            )
        # master's object: it cannot be combined with fields of the new reader
        return self._legacy(getter, args, kwargs)

    def _derived(self, name, time, extent, selection=None, **options):
        """Result of recipe ``name`` (amr.recipes) on the new reader.

        selection: a MasterSelection (selection_box=) choosing the input patches, or None.
        """
        from pyphare.pharesee.amr.recipes import RECIPES, evaluate

        recipe = RECIPES[name]
        paths, formats = {}, {}
        for stem in recipe.inputs:
            route = self._route(stem)
            paths[stem], formats[stem] = route.path, route.fmt
        first = paths[recipe.inputs[0]]
        dl0 = self._backend(first).meta().dl0

        def opener(stem, times, ext, **kw):
            return self._open_path(paths[stem], times, ext, selection, **kw)

        sim = None
        if formats[recipe.inputs[0]] == "h5":
            def sim():
                return self._sim(first).load()

        return evaluate(recipe, opener, time, extent, dl0, formats, sim, **options)

    def _new_meta(self):
        """FileMeta of the first file the new reader serves (None if there is none)."""
        from pyphare.pharesee.amr.errors import FormatError

        order = (self._format,) if self._format else ("vtkhdf", "h5")
        for fmt in order:
            if fmt == "h5" and self._h5_legacy:
                continue
            for path in self._has[fmt]:
                try:
                    return self._backend(path).meta()
                except FormatError:
                    continue
        return None

    # the data is not in a vtkhdf file -> ValueError; on h5 the new reader does
    # not implement it -> master's reader serves it.
    def _unavailable(self, getter, stem, what, args, kwargs):
        """Getters the new reader does not provide: master's reader on h5, an error on vtkhdf."""
        if not self._all_legacy:
            try:
                fmt = self._route(stem).fmt
            except FileNotFoundError:
                if not (self._has["vtkhdf"] and self._format != "h5"):
                    raise
                fmt = "vtkhdf"
            if fmt == "vtkhdf":
                raise ValueError(
                    f"{getter}: the vtkhdf writer does not write {what}; read it from phareh5 "
                    "output"
                )
        return self._legacy(getter, args, kwargs)

    # ------------------------------------------------------------------ stored quantities
    def GetB(self, time, merged=False, interp="nearest", all_primal=True, **kwargs):
        """Magnetic field (EM_B) -> VectorField ("x", "y", "z").

        New reader: lazy node-centred field; keywords box=, selection_box=.
        all_primal=False: master's Yee hierarchy on h5 (legacy object), ValueError on vtkhdf.
        """
        return self._stored("GetB", "EM_B", (time, merged, interp, all_primal), kwargs)

    def GetE(self, time, merged=False, interp="nearest", all_primal=True, **kwargs):
        """Electric field (EM_E) -> VectorField; same rules as GetB."""
        return self._stored("GetE", "EM_E", (time, merged, interp, all_primal), kwargs)

    def GetMassDensity(self, time, merged=False, interp="nearest", **kwargs):
        """Ion mass density (ions_mass_density) -> ScalarField."""
        return self._stored("GetMassDensity", "ions_mass_density", (time, merged, interp),
                            kwargs)

    def GetNi(self, time, merged=False, interp="nearest", **kwargs):
        """Ion charge density (ions_charge_density) -> ScalarField."""
        return self._stored("GetNi", "ions_charge_density", (time, merged, interp), kwargs)

    def GetN(self, time, pop_name, merged=False, interp="nearest", **kwargs):
        """Density of population ``pop_name`` (ions_pop_<pop>_density) -> ScalarField."""
        return self._stored("GetN", f"ions_pop_{pop_name}_density",
                            (time, pop_name, merged, interp), kwargs)

    def GetVi(self, time, merged=False, interp="nearest", **kwargs):
        """Ion bulk velocity (ions_bulkVelocity) -> VectorField."""
        return self._stored("GetVi", "ions_bulkVelocity", (time, merged, interp), kwargs)

    def GetChargeFlux(self, time, merged=False, interp="nearest", **kwargs):
        """Ion charge flux (ions_charge_flux) -> VectorField."""
        return self._stored("GetChargeFlux", "ions_charge_flux", (time, merged, interp), kwargs)

    def GetFlux(self, time, pop_name, merged=False, interp="nearest", **kwargs):
        """Flux of population ``pop_name`` (ions_pop_<pop>_flux) -> VectorField."""
        return self._stored("GetFlux", f"ions_pop_{pop_name}_flux",
                            (time, pop_name, merged, interp), kwargs)

    def GetMHDrho(self, time, merged=False, interp="nearest", all_primal=True, **kwargs):
        """MHD density (mhd_rho) -> ScalarField; rules of GetB."""
        return self._stored("GetMHDrho", "mhd_rho", (time, merged, interp, all_primal), kwargs)

    def GetMHDV(self, time, merged=False, interp="nearest", all_primal=True, **kwargs):
        """MHD velocity (mhd_V) -> VectorField; rules of GetB."""
        return self._stored("GetMHDV", "mhd_V", (time, merged, interp, all_primal), kwargs)

    def GetMHDP(self, time, merged=False, interp="nearest", all_primal=True, **kwargs):
        """MHD pressure (mhd_P) -> ScalarField; rules of GetB."""
        return self._stored("GetMHDP", "mhd_P", (time, merged, interp, all_primal), kwargs)

    def GetMHDrhoV(self, time, merged=False, interp="nearest", all_primal=True, **kwargs):
        """MHD momentum (mhd_rhoV) -> VectorField; rules of GetB."""
        return self._stored("GetMHDrhoV", "mhd_rhoV", (time, merged, interp, all_primal),
                            kwargs)

    def GetMHDEtot(self, time, merged=False, interp="nearest", all_primal=True, **kwargs):
        """MHD total energy (mhd_Etot) -> ScalarField; rules of GetB."""
        return self._stored("GetMHDEtot", "mhd_Etot", (time, merged, interp, all_primal),
                            kwargs)

    # ------------------------------------------------------------------ derived quantities
    def GetJ(self, time, merged=False, interp="nearest", all_primal=True, **kwargs):
        """Current density curl(B) -> VectorField.

        Default stencil: ``"yee"`` on phareh5 files, ``"node"`` on vtkhdf; an explicit
        ``stencil=`` wins.
        ``stencil="yee"`` (phareh5 only): master's Yee stencil per patch, averaged to nodes,
        which reproduces master's (and the code's) J; vtkhdf files hold node values only and
        cannot provide it (ValueError).
        ``stencil="node"``: centred differences of the node-centred B, one-sided (3-point) at
        patch faces without a same-level neighbour, never coarse data. This is PHARE's (Yee)
        J smoothed by (1, 2, 1)/4 along each derivative direction: on under-resolved current
        sheets its peak is lower than the code's J (7-50 % lower on the sample runs: Jz peak
        ratio 0.83-0.93, Jy 0.50-0.62). It is the only J a vtkhdf file can give.
        ``J.stencil`` and ``repr(J)`` say which stencil was used, and ``J.plot()`` puts it in
        the default title.
        Computed at the call (in memory): pass ``box=`` for a region. With ``box=``, B is read
        on the box grown by one level-0 cell so that nodes strictly inside the box get the
        same values as without ``box=``. A whole-domain B this Run returned and the caller
        still holds is reused (no second read once it is loaded). all_primal=False on h5:
        master's Yee object.
        """
        args = (time, merged, interp, all_primal)
        route = self._new_route("GetJ", "EM_B", args, kwargs, extra=("stencil",))
        if route is None:
            return self._legacy("GetJ", args, kwargs)
        explicit_stencil = "stencil" in kwargs
        stencil = kwargs.pop("stencil", "yee" if route.fmt == "h5" else "node")
        if stencil not in ("node", "yee"):
            raise ValueError(f"GetJ(stencil={stencil!r}): use 'node' or 'yee'")
        if not all_primal:
            if explicit_stencil:
                raise TypeError("GetJ(all_primal=False) returns master's Yee J; drop stencil=")
            return self._yee_on_legacy("GetJ", route, args, kwargs)
        return self._derived_getter("GetJ", "J" if stencil == "node" else "J_yee", "EM_B",
                                    args, kwargs)

    def GetDivB(self, time, merged=False, interp="nearest", **kwargs):
        """Divergence of B -> ScalarField; ``stencil=`` alone chooses which one.

        Default stencil: ``"yee"`` on 2D phareh5 files, ``"node"`` otherwise (vtkhdf, and
        phareh5 in 1D/3D); an explicit ``stencil=`` wins.
        ``stencil="yee"`` (2D phareh5 only): master's per-cell Yee divB (cell-centred, with
        ghosts), the same object and values as master's GetDivB (served by master's reader,
        as ``GetB(all_primal=False)``). It takes master's keywords: ``selection_box=``, not
        ``box=`` nor ``mask_faces=`` (TypeError). vtkhdf files (node values only) and 1D/3D
        runs (master has no 1D/3D divB) raise ValueError. Its ghost cells at coarse-fine
        boundaries hold large values (interpolated coarse B is not divergence-free there):
        look at the interior of patches.
        ``stencil="node"`` (new reader, any dimension): a node-centred field
        (``d.stencil == "node"``), the sum of the centred node derivatives (vtkhdf), or
        PHARE's Yee divB averaged from the cells around each node (phareh5); both equal the
        node formula away from coarse-fine faces (nodes with same-level neighbours on every
        side). Keywords ``box=``, ``selection_box=``, ``mask_faces=``.
        ``mask_faces`` (node, default True, both formats): nodes where a one-sided derivative
        or ghost values across a coarse-fine face were needed (also non-periodic boundaries)
        become NaN. A node file has nothing across a coarse-fine face, and the ghost cells
        of an h5 file hold large divB there: divB at those nodes is not a property of B
        (76-119x the interior value on the sample runs).
        Reductions of a masked field: ``min()`` and ``max()`` ignore the NaN nodes
        (``abs(d).max()`` is max|divB| over the unmasked nodes); ``mean()`` and
        ``integral()`` return NaN, since part of the domain has no value; pass
        ``mask_faces=False`` for values everywhere. Computed at the call (pass ``box=`` for a
        region).
        """
        args = (time, merged, interp)
        route = self._new_route("GetDivB", "EM_B", args, kwargs,
                                extra=("mask_faces", "stencil"), check_interp=False)
        if route is None:
            return self._legacy("GetDivB", args, kwargs)
        ndim = self._backend(route.path).meta().dim
        h5_2d = route.fmt == "h5" and ndim == 2
        stencil = kwargs.pop("stencil", "yee" if h5_2d else "node")
        if stencil not in ("node", "yee"):
            raise ValueError(f"GetDivB(stencil={stencil!r}): use 'node' or 'yee'")
        if stencil == "yee":
            return self._master_divb(route, ndim, args, kwargs)
        _check_merged_interp("GetDivB", merged, interp, route.fmt == "h5")
        mask_faces = kwargs.pop("mask_faces", None)
        return self._derived_getter("GetDivB", "divB", "EM_B", args, kwargs,
                                    mask_faces=None if mask_faces is None else bool(mask_faces))

    def _master_divb(self, route, ndim, args, kwargs):
        """GetDivB(stencil="yee"): master's per-cell object, or why it cannot be had."""
        if route.fmt != "h5":
            raise ValueError(
                "GetDivB(stencil='yee') is master's per-cell Yee divB, which needs the "
                "Yee-centred B of a phareh5 file; a vtkhdf file holds node values only: use "
                "stencil='node' (the default on vtkhdf)"
            )
        if ndim != 2:
            raise ValueError(
                f"GetDivB(stencil='yee') is master's per-cell divB, which master computes in "
                f"2D only; this run is {ndim}D: use stencil='node' (the default in {ndim}D)"
            )
        if "box" in kwargs:
            raise TypeError(
                "GetDivB(box=...) with the Yee stencil (the 2D phareh5 default): master's "
                "per-cell divB is read by master's reader, which takes selection_box= "
                "instead; or pass stencil='node' for the node divB on the box"
            )
        if "mask_faces" in kwargs:
            raise TypeError(
                "GetDivB(mask_faces=...) applies to the node divB: pass stencil='node' "
                "(master's per-cell divB, the 2D phareh5 default, has no mask)"
            )
        # the new reader has no cell-centred derived fields (GetB(all_primal=False) is
        # master's too)
        return self._legacy("GetDivB", args, kwargs)

    def GetPe(self, time, merged=False, interp="nearest", all_primal=True, **kwargs):
        """Electron pressure Ni * Te -> ScalarField (isothermal closure).

        New reader on phareh5: Te from the simulation stored in the file, as master
        (run.py:164-173); box= accepted; computed at the call (from a held Ni if any).
        vtkhdf: ValueError (no simulation parameters in the file). Legacy path: master's
        GetPe, which takes no other keyword.
        """
        args = (time, merged, interp, all_primal)
        route = self._new_route("GetPe", "ions_charge_density", args, kwargs)
        if route is None:
            return self._legacy("GetPe", args, kwargs)
        if not all_primal:
            return self._yee_on_legacy("GetPe", route, args, kwargs)
        return self._derived_getter("GetPe", "Pe", "ions_charge_density", args, kwargs)

    def GetMagneticFlux(self, time, interp="nearest", xn=None, yn=None, Xn=None, Yn=None):
        """Magnetic flux function ``(Az, (xn, yn))`` of a 2D run (master's algorithm).

        New reader: B painted on the node grid of the finest level (``finest()``), as master's
        merged finest data; coarse-fine gaps are filled by ``finest()`` (cubic), not by
        master's interpolators. Default grid: the finest level's nodes, where nothing is
        interpolated, so ``interp`` has no effect there (Az within 0.6 % of max|Az| of
        master's ``interp="bilinear"`` and 1-2 % of master's ``"nearest"`` on the sample
        runs). ``interp`` ("nearest", "linear"/"bilinear"; other values raise ValueError, as
        on master) samples the finest-level values on given ``xn, yn`` grids.
        """
        if self._route("EM_B").legacy:
            return self._legacy("GetMagneticFlux", (time, interp, xn, yn, Xn, Yn), {})
        return self._derived("Az", time, None, interp=interp, xn=xn, yn=yn, Xn=Xn, Yn=Yn)

    def GetReconnectionRate(self, times, interp="nearest"):
        """``(t_centred, rates, flux_at_xpoint, xpoint_trajectory)``, master's algorithm.

        New reader: the flux function of every time from ``GetMagneticFlux``, sampled as on
        master on one grid, the finest level's nodes at the first time (a change of the
        finest level between dumps does not change the sampling).
        """
        if self._route("EM_B").legacy:
            return self._legacy("GetReconnectionRate", (times, interp), {})
        times = np.asarray(times, dtype=np.float64)
        flux_at_xpoint, xpoint_trajectory = [], []
        grid = {}
        for t in times:
            Az, (xn, yn) = self.GetMagneticFlux(t, interp=interp, **grid)
            if not grid:
                Xn, Yn = np.meshgrid(xn, yn, indexing="ij")
                grid = dict(xn=xn, yn=yn, Xn=Xn, Yn=Yn)
            x_xp, y_xp, idx = self.FindPrimaryXPoint(Az, xn, yn)
            xpoint_trajectory.append([x_xp, y_xp])
            flux_at_xpoint.append(Az[idx])
        flux_at_xpoint = np.array(flux_at_xpoint)
        xpoint_trajectory = np.array(xpoint_trajectory)
        rates = np.diff(flux_at_xpoint) / np.diff(times)
        times_centered = (times[:-1] + times[1:]) / 2
        return times_centered, rates, flux_at_xpoint, xpoint_trajectory

    # ------------------------------------------------------------------ not on the new reader
    def GetPressure(self, time, pop_name, merged=False, interp="nearest", **kwargs):
        """Pressure tensor of a population: master's reader (h5; ValueError on vtkhdf)."""
        return self._unavailable("GetPressure", f"ions_pop_{pop_name}_momentum_tensor",
                                 "momentum_tensor", (time, pop_name, merged, interp), kwargs)

    def GetPi(self, time, merged=False, interp="nearest", **kwargs):
        """Ion pressure tensor: master's reader (h5; ValueError on vtkhdf)."""
        return self._unavailable("GetPi", "ions_momentum_tensor", "momentum_tensor",
                                 (time, merged, interp), kwargs)

    def GetTags(self, time, merged=False, **kwargs):
        """Refinement tags: master's reader (h5; ValueError on vtkhdf)."""
        return self._unavailable("GetTags", "tags", "tags", (time, merged), kwargs)

    def GetRanks(self, time, merged=False, interp="nearest", **kwargs):
        """MPI rank of each patch: master's reader (h5; ValueError on vtkhdf)."""
        return self._unavailable("GetRanks", "EM_B", "MPI ranks", (time, merged, interp),
                                 kwargs)

    def GetParticles(self, time, pop_name, hier=None, **kwargs):
        """Particles: master's reader (h5; ValueError on vtkhdf)."""
        if isinstance(pop_name, (list, tuple)) and not pop_name:
            return hier  # nothing to read, as master
        first = pop_name[0] if isinstance(pop_name, (list, tuple)) else pop_name
        return self._unavailable("GetParticles", f"ions_pop_{first}_domain", "particles",
                                 (time, pop_name, hier), kwargs)

    def GetParticleCount(self, time, **kwargs):
        """Particle count per patch: master's reader (h5; ValueError on vtkhdf)."""
        return self._unavailable("GetParticleCount", "particle_count", "particle counts",
                                 (time,), kwargs)

    # ------------------------------------------------------------------ metadata
    def _pop_attribute(self, pop_name, attribute):
        from pyphare.pharesee.amr.errors import FormatError
        from pyphare.pharesee.amr.io.files import POP_KINDS, pop_stem

        values, sources = [], []
        for kind in POP_KINDS:
            try:
                route = self._route(pop_stem(pop_name, kind))
            except FileNotFoundError:
                continue
            try:
                meta = self._backend(route.path).meta()
                where, attrs = meta.path, meta.attrs
            except FormatError:
                if route.fmt != "h5":
                    continue
                # a file the new reader does not read (particles): its root attributes
                where, attrs = route.path, _h5_root_attrs(route.path)
            if attribute not in attrs:
                raise ValueError(f"{where} has no {attribute} attribute, it was probably "
                                 "written by an older PHARE version")
            values.append(attrs[attribute])
            sources.append(where)
        if not values:
            if self._has["h5"] and self._format != "vtkhdf":
                # only files the new reader does not read (particles): master's reading
                return pop_attribute_from_h5(self.path, pop_name, attribute)
            raise ValueError(f"no file of population {pop_name} found in {self.path}")
        if any(np.any(v != values[0]) for v in values):
            raise ValueError(f"inconsistent {attribute} across the files of population "
                             f"{pop_name}: {dict(zip(sources, values))}")
        return values[0]

    def GetMass(self, pop_name, **kwargs):
        """Mass of a population (``pop_mass`` attribute of its files, which must agree)."""
        if self._all_legacy:
            return LegacyRun.GetMass(self, pop_name, **kwargs)
        return self._pop_attribute(pop_name, "pop_mass")

    def GetCharge(self, pop_name, **kwargs):
        """Charge of a population (``pop_charge`` attribute of its files)."""
        if self._all_legacy:
            return LegacyRun.GetCharge(self, pop_name, **kwargs)
        return self._pop_attribute(pop_name, "pop_charge")

    def GetDomainSize(self, **kwargs):
        """Physical size of the domain, float64 ``(dim,)`` (cells of level 0 times dl)."""
        meta = None if self._all_legacy else self._new_meta()
        if meta is None:
            return LegacyRun.GetDomainSize(self, **kwargs)
        return np.asarray(meta.domain_cells0, dtype=np.float64) * np.asarray(meta.dl0)

    def GetDl(self, level=0):
        """Mesh size of ``level``, float64 ``(dim,)``."""
        meta = None if self._all_legacy else self._new_meta()
        if meta is None:
            return LegacyRun.GetDl(self, level)
        if level < 0:
            raise ValueError(f"Run::GetDl: level must be >= 0, got {level}")
        return np.asarray(meta.dl0, dtype=np.float64) / refinement_ratio**level

    def all_times(self):
        """``{key: sorted float64 times}`` for every quantity file (master's keys: "B", "E",
        "Ni", ... else the file stem). New reader: files written without data (empty
        ``*momentum_tensor.vtkhdf``) are left out; a file that cannot be opened is left out
        with a RuntimeWarning naming it and the error; a quantity present in both formats
        needs ``format=``."""
        if self._all_legacy:
            return LegacyRun.all_times(self)
        from pyphare.pharesee.amr.errors import FormatError
        from pyphare.pharesee.amr.io.files import scan, short_name

        # on a run that is not all-legacy, h5 files served by the legacy reader
        # list their times sorted (master: h5 group order); h5 files the new reader cannot
        # parse (particles, a tags.h5 without datasets) are listed the same way; vtkhdf files
        # without any dataset (empty *momentum_tensor.vtkhdf) are left out.
        # A file that cannot be opened (truncated, being written...) is skipped with a
        # warning naming it, so one bad file does not hide the times of the others.
        out = {}
        for stem in sorted(scan(self.path)):
            try:
                route = self._route(stem)
            except FileNotFoundError:  # format= chosen, and this stem has no such file
                continue
            try:
                if route.legacy:
                    out[short_name(stem)] = _times_of_file(route.path)
                else:
                    try:
                        out[short_name(stem)] = self._backend(route.path).times()
                    except FormatError:
                        if route.fmt != "h5":
                            raise
                        # particles, tags...: not read by the new reader; master lists them
                        out[short_name(stem)] = _times_of_file(route.path)
            except FormatError:
                continue
            except (OSError, ValueError, KeyError, RuntimeError) as e:
                self._backends.pop(route.path, None)
                warnings.warn(
                    f"all_times(): skipping {route.path}: it cannot be opened "
                    f"({type(e).__name__}: {e})",
                    RuntimeWarning,
                    stacklevel=2,
                )
        return out

    def times(self, qty):
        """Times of quantity ``qty`` ("B" or "EM_B", "Ni", ...): float64 array, sorted."""
        if self._all_legacy:
            return LegacyRun.times(self, qty)
        from pyphare.pharesee.amr.io.files import quantity_stem

        route = self._route(quantity_stem(qty))
        if route.legacy:
            return _times_of_file(route.path)
        return self._backend(route.path).times()

    # ------------------------------------------------------------------ planning
    def plan(self, qty, time, box=None):
        """What reading ``qty`` at ``time`` (on ``box``) would touch, without reading it.

        qty
            a stored quantity ("B", "EM_B", "Ni", "ions_pop_protons_density", ...) or "J",
            "divB" (B on the box widened by the stencil's halo: what ``GetJ(t)`` and
            ``GetDivB(t, stencil="node")`` read)

        Returns a ``RunPlan`` (bytes read from the file, patches per time and level). Reads only the index
        tables. ValueError when the quantity is served by the legacy reader.
        """
        from pyphare.pharesee.amr.boxes import PhysBox
        from pyphare.pharesee.amr.field import _missing_bytes
        from pyphare.pharesee.amr.io.files import quantity_stem
        from pyphare.pharesee.amr.recipes import RECIPES, VECTOR_NAMES, widen

        stem = "EM_B" if qty in ("J", "divB") else quantity_stem(qty)
        route = self._route(stem)
        if route.legacy:
            raise ValueError(f"plan(): {stem} is served by the legacy reader, which has no "
                             "plan; use Run(path, legacy=False)")
        meta = self._backend(route.path).meta()
        extent = PhysBox.parse(box, meta.dim)
        opening = {}
        if qty in ("J", "divB"):
            # the recipe GetJ / GetDivB(box=...) use: on h5 they read the Yee B through a
            # derived opening (amr.recipes._current_yee, _divergence), which does not share
            # the data of a held GetB
            name = "J_yee" if (qty == "J" and route.fmt == "h5") else qty
            extent = widen(extent, RECIPES[name].halo, meta.dl0)
            if route.fmt == "h5":
                opening = dict(derived=f"{qty}_yee", quantity=qty,
                               stencil="yee" if qty == "J" else "node",
                               names=VECTOR_NAMES if qty == "J" else ("value",))
        field = self._open_path(route.path, time, extent, **opening)
        nbytes, patches = 0, {}
        for snap in field.snapshots():
            for ld in snap.levels:
                nbytes += _missing_bytes(ld, None)
                patches[(float(ld.time), ld.geom.level)] = len(ld.missing()[0])
        return RunPlan(f"{qty} at t={time}", nbytes, patches)


def _h5_root_attrs(path):
    """Root attributes of the h5 file ``path``."""
    import h5py

    with h5py.File(path, "r") as f:
        return dict(f.attrs)


def _times_of_file(path):
    """master's way of listing the times of one h5 file (run.py:381-395), sorted."""
    import h5py

    with h5py.File(path, "r") as f:
        return np.sort(np.array([float(t) for t in f["t"]]))

