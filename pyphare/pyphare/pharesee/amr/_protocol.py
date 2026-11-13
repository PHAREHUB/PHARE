"""numpy protocol of AMR fields (DESIGN 7.5): operators, ufuncs, a whitelist of numpy functions.

``NumpyOps`` is a mixin of ``field.AMRField``. It only uses ``self.snapshots()``,
``self._like(snaps, ...)`` and the ``amr.ops`` modules (never ``field.py``).

What works
  * Python operators ``+ - * / // % ** << >> & | ^`` (and reflected forms), unary ``- + abs ~``,
    comparisons ``< <= > >=`` (elementwise, a field of bool). Operands: fields, Python/numpy
    scalars, 0-d arrays. Results are new in-memory fields; they chain.
  * numpy ufuncs called directly (``np.sqrt(B.x)``, ``np.add(a, b)``), one output only.
  * a whitelist of numpy functions mapped to field operations: ``np.mean``, ``np.nanmean``,
    ``np.min/amin/nanmin``, ``np.max/amax/nanmax`` (best-data reductions, see ``ops.reduce``),
    ``np.gradient`` (scalar -> ``grad``), ``np.copy``, ``np.clip``, ``np.where``,
    ``np.isclose``, ``np.allclose``, ``np.any``, ``np.all``.
What raises TypeError, with a message naming what to use instead
  * ``==`` / ``!=`` are NOT elementwise. Field vs field: identity, as on master's hierarchies
    (fields hash by identity, so sets and dicts of fields work). Field vs a number or an array
    (``N == 1.0``, ``N == N.max()``): TypeError, because a silent ``False`` would turn
    ``np.where(N == N.max(), 1, 0)`` into ``0``. Elementwise: ``np.equal(a, b)``,
    ``np.not_equal(a, b)``, ``np.isclose(a, b)`` (fields of bool); for a single number compare
    a reduction (``f.max() == x``). Other objects (None, strings, legacy hierarchies) compare
    by identity (False), so ``f in [None, f]`` keeps working.
    Numpy scalar on the LEFT (``np.float64(1.0) == f``): numpy itself turns this into
    ``np.equal(np.float64(1.0), f)`` before the field is asked, so the result is an
    elementwise field of bool (the same as ``np.equal``); the field cannot tell it apart from
    an explicit ``np.equal`` call without breaking numpy's protocol. Write ``np.equal(x, f)``
    to make the intent explicit.
  * ``np.sum`` (use ``integral()`` / ``mean()``), ``np.std``, ``np.var`` and every other numpy
    function (use ``finest()`` for a regular array); ``axis=``; ``out=``/``where=``;
    ufunc methods other than a plain call (``np.add.reduce``...); two-output ufuncs;
    ``bool(f)``, ``iter(f)``, ``len(f)``, ``np.asarray(f)``; arrays with ``ndim > 0``;
    objects from the legacy reader (``pyphare.pharesee.hierarchy``) mixed with new fields.
Several times: binary operations pair fields holding the same list of times (else
ValueError); two single-time fields combine whatever their times (the result has the first
operand's time); unary operations apply to every time; reductions need one time.
dtype: numpy rules (NEP 50) - a float32 field times a Python float stays float32, times a
``np.float64`` scalar becomes float64.
"""

from __future__ import annotations

import numbers

import numpy as np

_LEGACY_MODULE = "pyphare.pharesee.hierarchy"

_LEGACY_MSG = (
    "{what}: one object comes from the legacy reader (pyphare.pharesee.hierarchy, made by "
    "Run(path, legacy=True), a getter called with all_primal=False, or one the new reader "
    "does not provide, e.g. GetPressure) and the other from the new reader (pharesee.amr); "
    "they cannot be combined. Read both on the same path: e.g. drop all_primal=False, or "
    "drop legacy=True."
)

_FINEST_ADVICE = (
    "use field.finest() for a regular array (FinestGrid.data), or the field's methods "
    "integral(), mean(), min(), max()"
)


# --------------------------------------------------------------------------- operand checks


def _is_field(x):
    return isinstance(x, NumpyOps)


def _is_legacy(x):
    """True for objects of master's hierarchy classes (detected without importing them)."""
    module = getattr(type(x), "__module__", "") or ""
    return module.startswith(_LEGACY_MODULE)


def _is_scalar(x):
    if isinstance(x, (numbers.Number, np.generic)):
        return True
    return isinstance(x, np.ndarray) and x.ndim == 0


def _check_operand(x, what):
    """Raise TypeError for operands that can never combine with a field; True if usable."""
    if _is_field(x) or _is_scalar(x) or x is None:
        return True
    if _is_legacy(x):
        raise TypeError(_LEGACY_MSG.format(what=what))
    if isinstance(x, (np.ndarray, list, tuple)):
        shape = np.shape(x) if not isinstance(x, np.ndarray) else x.shape
        raise TypeError(
            f"{what}: an array of shape {shape} cannot be combined with an AMR field (its "
            "patches have no common array shape). Use a scalar, another field, or "
            "field.finest() to get a regular array first."
        )
    return False


def _ufunc_name(func):
    return getattr(func, "__name__", "operation")


def _time_tolerance(field, t):
    dtype = np.dtype(getattr(field, "_time_dtype", np.float64))
    return max(2.0 * float(np.spacing(np.abs(np.asarray(t, dtype=dtype)))), 5e-11)


def _common_count(fields, what):
    """Number of time pairs for these fields (rules in the module doc); ValueError otherwise."""
    # DESIGN-NOTE: two single-time fields combine even at different times (DESIGN 2
    # "arithmetic": B(t1) - B(t0)); only multi-time fields must hold the same time list.
    times = [f.times() for f in fields]
    if all(len(t) == 1 for t in times):
        return 1
    t0 = times[0]
    for f, t in zip(fields[1:], times[1:]):
        same = len(t) == len(t0) and all(
            abs(a - b) <= max(_time_tolerance(fields[0], a), _time_tolerance(f, b))
            for a, b in zip(t0, t)
        )
        if not same:
            raise ValueError(
                f"{what}: the fields hold different times {t0.tolist()} and {t.tolist()}; "
                "multi-time fields combine time by time and need the same times in the same "
                "order. Select one time with field.at(t), or read both at the same times."
            )
    return len(t0)


def _map(func, operands, what):
    """Apply ``func(*snapshots_or_scalars) -> Snapshot`` time by time -> new field.

    The result takes the first field's meta/stencil/sim (``_like``), quantity None.
    """
    # DESIGN-NOTE: results have quantity None (repr shows "(computed)"): "EM_B" on B + 1
    # would claim the file holds it.
    fields = [o for o in operands if _is_field(o)]
    n = _common_count(fields, what)
    snaps = []
    for k in range(n):
        args = [o.snapshots()[k] if _is_field(o) else o for o in operands]
        snaps.append(func(*args))
    return fields[0]._like(snaps, quantity=None)


def _apply(ufunc, operands):
    """Elementwise ``ufunc`` over operands (fields/scalars) -> field, or NotImplemented."""
    from .ops import arith

    what = _ufunc_name(ufunc)
    if not all(_check_operand(o, what) for o in operands):
        return NotImplemented
    if len(operands) == 1:
        return _map(lambda s: arith.unary(ufunc, s), operands, what)
    if len(operands) == 2:
        return _map(lambda a, b: arith.binary(ufunc, a, b), operands, what)
    raise TypeError(f"{what}: ufuncs with {len(operands)} inputs are not supported on AMR fields")


def _binop(ufunc):
    def forward(self, other):
        return _apply(ufunc, [self, other])

    def reflected(self, other):
        return _apply(ufunc, [other, self])

    name = ufunc.__name__
    forward.__doc__ = f"Elementwise np.{name}(self, other) -> new field."
    reflected.__doc__ = f"Elementwise np.{name}(other, self) -> new field."
    return forward, reflected


def _unop(ufunc):
    def op(self):
        return _apply(ufunc, [self])

    op.__doc__ = f"Elementwise np.{ufunc.__name__}(self) -> new field."
    return op


def _cmp(ufunc):
    def op(self, other):
        return _apply(ufunc, [self, other])

    op.__doc__ = f"Elementwise np.{ufunc.__name__}(self, other) -> field of bool."
    return op


# --------------------------------------------------------------------------- the mixin


_EQ_ADVICE = (
    "{op} between an AMR field and a {kind} is not elementwise (field == field is identity, "
    "as on master's hierarchies): use np.equal(f, x) / np.not_equal(f, x) or "
    "np.isclose(f, x) for a field of bool, or compare a reduction (f.max() == x)"
)


def _eq_operand_kind(other):
    """"number"/"array" for operands ``==`` refuses, None for those compared by identity."""
    if _is_scalar(other):
        return "number"
    if isinstance(other, (np.ndarray, list, tuple)):
        return "array"
    return None


class NumpyOps:
    """Operators and numpy dispatch for AMR fields (rules in the module docstring).

    ``==``/``!=``: identity between fields (``__hash__`` is identity too); TypeError against a
    number or an array (use ``np.equal``/``np.isclose``); identity (False) otherwise.
    """

    def __eq__(self, other):
        """Identity with another field; TypeError against a number or an array."""
        if _is_field(other):
            return self is other
        kind = _eq_operand_kind(other)
        if kind is not None:
            raise TypeError(_EQ_ADVICE.format(op="==", kind=kind))
        return NotImplemented  # Python falls back to identity: False

    def __ne__(self, other):
        """Negation of ``__eq__`` (TypeError against a number or an array)."""
        if _is_field(other):
            return self is not other
        kind = _eq_operand_kind(other)
        if kind is not None:
            raise TypeError(_EQ_ADVICE.format(op="!=", kind=kind))
        return NotImplemented

    __hash__ = object.__hash__

    __add__, __radd__ = _binop(np.add)
    __sub__, __rsub__ = _binop(np.subtract)
    __mul__, __rmul__ = _binop(np.multiply)
    __truediv__, __rtruediv__ = _binop(np.true_divide)
    __floordiv__, __rfloordiv__ = _binop(np.floor_divide)
    __mod__, __rmod__ = _binop(np.remainder)
    __pow__, __rpow__ = _binop(np.power)
    __lshift__, __rlshift__ = _binop(np.left_shift)
    __rshift__, __rrshift__ = _binop(np.right_shift)
    __and__, __rand__ = _binop(np.bitwise_and)
    __or__, __ror__ = _binop(np.bitwise_or)
    __xor__, __rxor__ = _binop(np.bitwise_xor)

    __lt__ = _cmp(np.less)
    __le__ = _cmp(np.less_equal)
    __gt__ = _cmp(np.greater)
    __ge__ = _cmp(np.greater_equal)

    __neg__ = _unop(np.negative)
    __pos__ = _unop(np.positive)
    __abs__ = _unop(np.absolute)
    __invert__ = _unop(np.invert)

    def __matmul__(self, other):
        raise TypeError("@ is not defined on AMR fields: use dot(a, b) for the scalar product "
                        "of two vector fields, cross(a, b) for the vector product")

    __rmatmul__ = __matmul__

    def __divmod__(self, other):
        raise TypeError("divmod() gives two results, not supported on AMR fields: "
                        "use a // b and a % b")

    __rdivmod__ = __divmod__

    # ---- conversions that have no meaning on an AMR field
    def __bool__(self):
        raise TypeError(
            "the truth value of an AMR field is ambiguous: use np.any(f) or np.all(f), or "
            "compare a reduction (f.max() > 0)"
        )

    def __iter__(self):
        raise TypeError(
            "AMR fields cannot be iterated: use field.times() for the times, "
            "field.levels_data() / field.levels() for the patches, or field.finest()"
        )

    def __len__(self):
        raise TypeError(
            "AMR fields have no length: use field.nbytes(), field.times(), or the patch "
            "counts of field.levels_data()"
        )

    def __array__(self, dtype=None, copy=None):
        raise TypeError(
            "an AMR field is not a numpy array (patches of several levels): use "
            "field.finest().data for a regular array, or field.levels_data()[k].buffer(name) "
            "for the flat values of one level"
        )

    # ---- numpy protocol
    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        name = _ufunc_name(ufunc)
        if method != "__call__":
            hint = ("use field.integral(), mean(), min() or max() for reductions"
                    if method == "reduce" else _FINEST_ADVICE)
            raise TypeError(f"np.{name}.{method} is not supported on AMR fields: {hint}")
        if "out" in kwargs:
            raise TypeError(f"np.{name}(..., out=) is not supported on AMR fields: results are "
                            "new fields, use the returned value (a = np.{name}(...))")
        if "where" in kwargs:
            raise TypeError(f"np.{name}(..., where=) is not supported on AMR fields: use "
                            "np.where(cond, np.{name}(...), other)")
        if kwargs:
            raise TypeError(f"np.{name}: keywords {sorted(kwargs)} are not supported on AMR "
                            "fields")
        if ufunc.nout != 1:
            raise TypeError(f"np.{name} has {ufunc.nout} outputs, not supported on AMR fields; "
                            "compute each output with a one-output ufunc")
        return _apply(ufunc, list(inputs))

    def __array_function__(self, func, types, args, kwargs):
        handler = _HANDLED.get(func)
        if handler is None:
            raise TypeError(_unsupported_message(func))
        return handler(*args, **kwargs)


# --------------------------------------------------------------------------- numpy functions


def _unsupported_message(func):
    name = getattr(func, "__name__", str(func))
    if name == "sum":
        return ("np.sum is not defined on AMR fields (levels overlap and cells differ in size): "
                "use field.integral() for the integral over the domain/box, field.mean() for "
                "the mean, or field.mean(per_level=True) for master's per-level means")
    if name in ("std", "var", "nanstd", "nanvar"):
        return (f"np.{name} is not supported on AMR fields: compute "
                "((f - f.mean())**2).mean() for the weighted variance, or use field.finest() "
                "for a regular array")
    return f"np.{name} is not supported on AMR fields: {_FINEST_ADVICE}"


def _require_field(x, fname):
    if _is_legacy(x):
        raise TypeError(_LEGACY_MSG.format(what=f"np.{fname}"))
    if not _is_field(x):
        raise TypeError(f"np.{fname}: expected an AMR field, got {type(x).__name__}")


def _no_axis(fname, axis, kwargs):
    if axis is not None:
        raise TypeError(
            f"np.{fname}(..., axis=) is not supported on AMR fields (patches have no array "
            "axes): use the field's own method, per_level=True, or field.finest()"
        )
    if kwargs:
        raise TypeError(f"np.{fname}: keywords {sorted(kwargs)} are not supported on AMR fields")


def _reducer(method, fname, nan=False):
    # DESIGN-NOTE: the nan* variants map to the same best-data reduction; they do not skip
    # NaN (the reduction kernels do not), so a NaN result raises with advice.
    def handler(a, axis=None, **kwargs):
        _require_field(a, fname)
        _no_axis(fname, axis, kwargs)
        value = getattr(a, method)()
        if nan and np.isnan(np.asarray(value, dtype=np.float64)).any():
            raise ValueError(
                f"np.{fname}: the field holds NaN and AMR reductions do not skip it; "
                "replace it first, e.g. np.where(np.isnan(f), 0.0, f), or use finest()"
            )
        return value

    return handler


def _gradient(f, *varargs, **kwargs):
    from .ops import deriv

    _require_field(f, "gradient")
    if varargs or kwargs:
        raise TypeError("np.gradient on an AMR field takes no spacing, axis or edge_order: "
                        "the mesh sizes come from the file (level-aware grad)")
    if tuple(f.names) != ("value",):
        raise TypeError(f"np.gradient needs a scalar field, got components {f.names}: use "
                        "grad() on a component (np.gradient(B.x)), or curl()/div() on a vector")
    return _map(deriv.grad, [f], "gradient")


def _copy(a, **kwargs):
    _require_field(a, "copy")
    if kwargs:
        raise TypeError(f"np.copy: keywords {sorted(kwargs)} are not supported on AMR fields")
    return a.copy()


_NO = object()


def _clip(a, a_min=_NO, a_max=_NO, out=None, *, min=_NO, max=_NO, **kwargs):
    from .ops import arith

    _require_field(a, "clip")
    if out is not None:
        raise TypeError("np.clip(..., out=) is not supported on AMR fields: use the result")
    if kwargs:
        raise TypeError(f"np.clip: keywords {sorted(kwargs)} are not supported on AMR fields")
    lo = a_min if a_min is not _NO else (min if min is not _NO else None)
    hi = a_max if a_max is not _NO else (max if max is not _NO else None)
    for bound in (lo, hi):
        _check_operand(bound, "clip")
    return _map(arith.clip, [a, lo, hi], "clip")


def _where(cond, x=_NO, y=_NO):
    from .ops import arith

    if x is _NO or y is _NO:
        raise TypeError("np.where(cond) with one argument returns array indices, which AMR "
                        "fields do not have: use np.where(cond, a, b), or field.finest()")
    for o in (cond, x, y):
        if not _check_operand(o, "where") or o is None:
            raise TypeError(f"np.where: cannot use {type(o).__name__} with AMR fields")
    if not any(_is_field(o) for o in (cond, x, y)):
        raise TypeError("np.where: at least one argument must be an AMR field")
    return _map(arith.where, [cond, x, y], "where")


def _isclose_func(rtol, atol, equal_nan):
    def isclose(x, y):
        return np.isclose(x, y, rtol=rtol, atol=atol, equal_nan=equal_nan)

    return isclose


def _isclose(a, b, rtol=1e-05, atol=1e-08, equal_nan=False):
    from .ops import arith

    for o in (a, b):
        if not _check_operand(o, "isclose") or o is None:
            raise TypeError(f"np.isclose: cannot use {type(o).__name__} with AMR fields")
    func = _isclose_func(rtol, atol, equal_nan)
    return _map(lambda x, y: arith.binary(func, x, y), [a, b], "isclose")


def _bool_reduce(a, how, fname):
    """``np.any``/``np.all`` over every stored node, component and time of a field."""
    # DESIGN-NOTE: any/all look at every stored node (coarse nodes under fine ones included)
    # and at every time of a multi-time field: a boolean answer is not ambiguous.
    func = np.any if how == "any" else np.all
    values = [bool(func(b)) for s in a.snapshots() for ld in s.levels
              for b in ld.buffers().values()]
    return bool(any(values)) if how == "any" else bool(all(values))


def _any_all(how):
    def handler(a, axis=None, **kwargs):
        _require_field(a, how)
        _no_axis(how, axis, kwargs)
        return _bool_reduce(a, how, how)

    return handler


def _allclose(a, b, rtol=1e-05, atol=1e-08, equal_nan=False):
    return _bool_reduce(_isclose(a, b, rtol, atol, equal_nan), "all", "allclose")


_HANDLED = {
    np.mean: _reducer("mean", "mean"),
    np.nanmean: _reducer("mean", "nanmean", nan=True),
    np.min: _reducer("min", "min"),
    np.amin: _reducer("min", "amin"),
    np.nanmin: _reducer("min", "nanmin", nan=True),
    np.max: _reducer("max", "max"),
    np.amax: _reducer("max", "amax"),
    np.nanmax: _reducer("max", "nanmax", nan=True),
    np.gradient: _gradient,
    np.copy: _copy,
    np.clip: _clip,
    np.where: _where,
    np.isclose: _isclose,
    np.allclose: _allclose,
    np.any: _any_all("any"),
    np.all: _any_all("all"),
}
