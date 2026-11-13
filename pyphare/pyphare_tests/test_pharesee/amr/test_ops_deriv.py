#!/usr/bin/env python3
"""Tests of pharesee.amr.ops.deriv: level-aware derivatives on flat AMR data.

The kernels are compared with per-patch reference loops written here (bit-identical),
whose face pass finds neighbours by brute-force global node keys (not through FaceTable).
"""

import unittest

import numpy as np

from pyphare.pharesee.amr.boxes import fortran_strides
from pyphare.pharesee.amr.ops import deriv
from pyphare.pharesee.amr.store import LevelData, Snapshot
from pyphare.pharesee.amr.testing import fill, make_amr


def layouts(periodic_only=False, nonperiodic_only=False):
    """(name, HierGeometry): 1D/2D/3D, 1-3 levels, 2-node patches, singleton shapes."""
    out = []
    for dim in (1, 2, 3):
        if not periodic_only:
            out.append((f"{dim}d-3lvl", make_amr(dim, 3, seed=dim)))
            out.append((f"{dim}d-1lvl-2node", make_amr(dim, 1, seed=30 + dim, two_node=True)))
            out.append((f"{dim}d-singletons", make_amr(dim, 2, seed=20 + dim, singletons=True,
                                                       two_node=True)))
        if not nonperiodic_only:
            out.append((f"{dim}d-2node-periodic", make_amr(dim, 2, seed=10 + dim, two_node=True,
                                                           periodic=True)))
            out.append((f"{dim}d-mixed-periodic",
                        make_amr(dim, 3, seed=40 + dim, singletons=True,
                                 periodic=[True] + [False] * (dim - 1))))
    return out


# --------------------------------------------------------------------------- helpers


def global_index(geom):
    """(dim, nnodes) int64 global node index of every flat node."""
    return np.stack([np.rint(geom.node_coords(a) / geom.dl[a]).astype(np.int64)
                     for a in range(geom.dim)])


def node_key(geom, G):
    """Integer key of global node indices (periodic axes taken modulo N); -1 if outside."""
    N = geom.domain_cells
    key = np.zeros(G.shape[1], dtype=np.int64)
    bad = np.zeros(G.shape[1], dtype=bool)
    mult = 1
    for a in range(geom.dim):
        g = G[a] % N[a] if geom.periodic[a] else G[a]
        bad |= (g < 0) | (g > N[a])
        key += g * mult
        mult *= int(N[a]) + 1
    key[bad] = -1
    return key


def random_level(geom, seed, dtype=np.float64, names=("value",)):
    """Random value per GLOBAL node (equal on shared nodes and periodic images)."""
    key = node_key(geom, global_index(geom))
    size = int(np.prod(geom.domain_cells + 1))
    arrays = {}
    for k, n in enumerate(names):
        table = np.random.default_rng(seed + 7 * geom.level + 101 * k).standard_normal(size)
        arrays[n] = table[key].astype(dtype)
    return LevelData.from_arrays(geom, 0.0, arrays)


def patch_of_node(geom):
    return np.repeat(np.arange(geom.npatch), geom.size)


def reference_faces(geom, axis):
    """Face nodes having a same-level neighbour, found by brute force.

    Returns (face, inner, outer, sign) flat indices. A neighbour exists when some patch
    (or periodic image) holds the edge between the face node and the node one step out.
    """
    G = global_index(geom)
    patch = patch_of_node(geom)
    loc = G[axis] - geom.lower[patch, axis]
    n = geom.shape[patch, axis]
    stride = fortran_strides(geom.shape)[patch, axis]
    keys = node_key(geom, G)
    edge_keys = np.unique(keys[loc < n - 1])  # edge from node to node + e inside a patch
    uk, first = np.unique(keys, return_index=True)  # any holder of a node
    out = []
    for sign, on_face in ((-1, loc == 0), (1, loc == n - 1)):
        face = np.nonzero(on_face)[0]
        Go = G[:, face].copy()
        Go[axis] += sign
        okey = node_key(geom, Go)
        Ge = G[:, face].copy()
        if sign < 0:
            Ge[axis] -= 1
        ekey = node_key(geom, Ge)
        has = (okey >= 0) & (ekey >= 0) & np.isin(ekey, edge_keys)
        face = face[has]
        outer = first[np.searchsorted(uk, okey[has])]
        inner = face - sign * stride[face]
        out.append((face, inner, outer, np.full(len(face), sign)))
    return tuple(np.concatenate(x) for x in zip(*out))


def ref_ddx(geom, f, axis):
    """Per-patch reference of ddx (same formulas, same floating operations)."""
    dt = f.dtype
    two_h, h = dt.type(2.0 * geom.dl[axis]), dt.type(geom.dl[axis])
    d = np.empty_like(f)
    for slot in range(geom.npatch):
        v = np.moveaxis(geom.view(f, slot), axis, 0)
        r = np.moveaxis(geom.view(d, slot), axis, 0)
        r[1:-1] = (v[2:] - v[:-2]) / two_h
        if v.shape[0] >= 3:
            r[0] = (-3 * v[0] + 4 * v[1] - v[2]) / two_h
            r[-1] = (3 * v[-1] - 4 * v[-2] + v[-3]) / two_h
        else:
            r[0] = r[1] = (v[1] - v[0]) / h
    face, inner, outer, sign = reference_faces(geom, axis)
    plus = np.where(sign > 0, outer, inner)
    minus = np.where(sign > 0, inner, outer)
    d[face] = (f[plus] - f[minus]) / two_h
    return d


def ref_d2dx2(geom, f, axis):
    dt = f.dtype
    h2 = dt.type(geom.dl[axis] * geom.dl[axis])
    d = np.empty_like(f)
    for slot in range(geom.npatch):
        v = np.moveaxis(geom.view(f, slot), axis, 0)
        r = np.moveaxis(geom.view(d, slot), axis, 0)
        n = v.shape[0]
        r[1:-1] = (v[2:] - 2 * v[1:-1] + v[:-2]) / h2
        if n >= 4:
            r[0] = (2 * v[0] - 5 * v[1] + 4 * v[2] - v[3]) / h2
            r[-1] = (2 * v[-1] - 5 * v[-2] + 4 * v[-3] - v[-4]) / h2
        elif n == 3:
            r[0] = r[-1] = r[1]
        else:
            r[0] = r[1] = 0
    face, inner, outer, sign = reference_faces(geom, axis)
    plus = np.where(sign > 0, outer, inner)
    minus = np.where(sign > 0, inner, outer)
    d[face] = (f[plus] - 2 * f[face] + f[minus]) / h2
    return d


def lacks_formula(geom, axis, min_nodes):
    """Bool per node: end-layer node of a patch with < min_nodes along axis, no neighbour."""
    patch = patch_of_node(geom)
    G = global_index(geom)
    loc = G[axis] - geom.lower[patch, axis]
    n = geom.shape[patch, axis]
    bad = ((loc == 0) | (loc == n - 1)) & (n < min_nodes)
    bad[geom.faces(axis).face] = False
    return bad


def snapshot(hier, funcs, dtype=np.float64):
    return Snapshot(hier, tuple(fill(hier, funcs, dtype=dtype)))


# --------------------------------------------------------------------------- tests


class TestAgainstReference(unittest.TestCase):
    def test_ddx_bit_identical(self):
        for name, hier in layouts():
            for dtype in (np.float64, np.float32):
                for geom in hier.levels:
                    ld = random_level(geom, 1, dtype)
                    f = ld.buffer("value")
                    for axis in range(geom.dim):
                        got = deriv.ddx(ld, "value", axis)
                        self.assertEqual(got.dtype, dtype)
                        np.testing.assert_array_equal(got, ref_ddx(geom, f, axis),
                                                      err_msg=f"{name} loc{geom.level} a{axis}")

    def test_d2dx2_bit_identical(self):
        for name, hier in layouts():
            for geom in hier.levels:
                ld = random_level(geom, 2, np.float32)
                f = ld.buffer("value")
                for axis in range(geom.dim):
                    np.testing.assert_array_equal(deriv.d2dx2(ld, "value", axis),
                                                  ref_d2dx2(geom, f, axis),
                                                  err_msg=f"{name} loc{geom.level} a{axis}")

    def test_small_blocks_same_result(self):
        from unittest import mock

        from pyphare.pharesee.amr.ops import arith

        hier = make_amr(3, 2, seed=5, two_node=True, periodic=True)
        ld = random_level(hier.levels[1], 3)
        ref = [deriv.ddx(ld, "value", a) for a in range(3)]
        ref2 = [deriv.d2dx2(ld, "value", a) for a in range(3)]
        with mock.patch.object(arith, "BLOCK_NODES", 37):
            for a in range(3):
                np.testing.assert_array_equal(deriv.ddx(ld, "value", a), ref[a])
                np.testing.assert_array_equal(deriv.d2dx2(ld, "value", a), ref2[a])

    def test_shared_nodes_same_derivative(self):
        for name, hier in layouts():
            for geom in hier.levels:
                ld = random_level(geom, 3)
                keys = node_key(geom, global_index(geom))
                for axis in range(geom.dim):
                    d = deriv.ddx(ld, "value", axis)
                    centred = ~lacks_formula(geom, axis, 10**9)  # not one-sided along axis
                    k, v = keys[centred], d[centred]
                    order = np.argsort(k, kind="stable")
                    k, v = k[order], v[order]
                    same = k[1:] == k[:-1]
                    np.testing.assert_array_equal(v[1:][same], v[:-1][same],
                                                  err_msg=f"{name} loc{geom.level} a{axis}")
                    if geom.level == 0 and geom.npatch > 1:
                        self.assertTrue(same.any(), name)


class TestExactness(unittest.TestCase):
    def test_linear_exact_everywhere(self):
        slopes = (1.5, -2.25, 0.75)
        for name, hier in layouts(nonperiodic_only=True):
            dim = hier.dim
            s = snapshot(hier, lambda *X, dim=dim: 3.0 + sum(slopes[a] * X[a] for a in range(dim)))
            for ld in s.levels:
                for axis in range(3):
                    d = deriv.ddx(ld, "value", axis)
                    expect = slopes[axis] if axis < dim else 0.0
                    np.testing.assert_allclose(d, expect, rtol=0, atol=1e-10, err_msg=name)

    def test_quadratic_exact_except_two_node_ends(self):
        for name, hier in layouts(nonperiodic_only=True):
            dim = hier.dim

            def f(*X, dim=dim):
                return sum((a + 1) * X[a] ** 2 for a in range(dim)) + X[0] * X[-1]

            def df(axis, X, dim=dim):
                r = 2 * (axis + 1) * X[axis]
                if dim > 1:
                    r = r + (X[-1] if axis == 0 else X[0] if axis == dim - 1 else 0)
                else:
                    r = r + 2 * X[0]
                return r

            s = snapshot(hier, f)
            for ld in s.levels:
                X = [ld.geom.node_coords(a) for a in range(dim)]
                for axis in range(dim):
                    d = deriv.ddx(ld, "value", axis)
                    ok = ~lacks_formula(ld.geom, axis, 3)
                    scale = np.abs(ld.buffer("value")).max() / ld.geom.dl[axis]
                    err = np.abs(d - df(axis, X))[ok].max() / scale
                    self.assertLess(err, 1e-13, f"{name} loc{ld.level} a{axis}")
                    # 2-node ends without neighbour really are first order (not exact)
                    bad = lacks_formula(ld.geom, axis, 3)
                    if bad.any():
                        self.assertGreater(np.abs(d - df(axis, X))[bad].max() / scale, 1e-6)

    def test_laplacian_exact_on_cubics(self):
        for name, hier in layouts(nonperiodic_only=True):
            dim = hier.dim
            c = (1.0, -2.0, 0.5)

            def f(*X, dim=dim, c=c):
                return sum(c[a] * X[a] ** 3 for a in range(dim)) + X[0] ** 2 * X[-1]

            def lap(*X, dim=dim, c=c):
                r = sum(6 * c[a] * X[a] for a in range(dim))
                return r + (6 * X[0] if dim == 1 else 2 * X[-1])  # 1D: x**2 * x = x**3

            s = snapshot(hier, f)
            out = deriv.laplacian(s)
            for ld, lo in zip(s.levels, out.levels):
                X = [ld.geom.node_coords(a) for a in range(dim)]
                ok = np.ones(ld.geom.nnodes, dtype=bool)
                for axis in range(dim):
                    ok &= ~lacks_formula(ld.geom, axis, 4)
                scale = np.abs(ld.buffer("value")).max() / ld.geom.dl.min() ** 2
                err = np.abs(lo.buffer("value") - lap(*X))[ok].max() / scale
                self.assertLess(err, 1e-12, f"{name} loc{ld.level}")

    def test_periodic_wrap_1d(self):
        hier = make_amr(1, 1, seed=3, periodic=True)
        geom = hier.levels[0]
        ld = random_level(geom, 4)
        f, d = ld.buffer("value"), deriv.ddx(ld, "value", 0)
        gi = global_index(geom)[0]
        N, h = int(geom.domain_cells[0]), geom.dl[0]
        val = {int(g) % N: v for g, v in zip(gi, f)}
        for g in (0, N):  # both images of the boundary node are centred across the period
            i = int(np.nonzero(gi == g)[0][0])
            self.assertEqual(d[i], (val[(g + 1) % N] - val[(g - 1) % N]) / (2 * h))
        self.assertFalse(deriv.one_sided_mask(geom).any())

    def test_float32_and_int_dtypes(self):
        hier = make_amr(2, 2, seed=1)
        s32 = snapshot(hier, lambda x, y: x + y, dtype=np.float32)
        self.assertEqual(deriv.ddx(s32.levels[0], "value", 0).dtype, np.float32)
        self.assertEqual(deriv.grad(s32).levels[1].dtype, np.float32)
        si = snapshot(hier, lambda x, y: np.rint(10 * x), dtype=np.int64)
        self.assertEqual(deriv.ddx(si.levels[0], "value", 0).dtype, np.float64)


class TestOperators(unittest.TestCase):
    def test_grad_components_and_missing_axes(self):
        for dim in (1, 2, 3):
            hier = make_amr(dim, 2, seed=dim)
            s = snapshot(hier, lambda *X, dim=dim: sum((a + 2.0) * X[a] for a in range(dim)))
            g = deriv.grad(s)
            self.assertEqual(g.names, ("x", "y", "z"))
            for ld in g.levels:
                for a, c in enumerate("xyz"):
                    np.testing.assert_allclose(ld.buffer(c), a + 2.0 if a < dim else 0.0,
                                               atol=1e-10)

    def test_curl_linear(self):
        for dim in (1, 2, 3):
            hier = make_amr(dim, 2, seed=dim + 3)
            # A = (y, z, x) restricted to existing axes, plus a 3x term in Az
            comp = {
                "x": lambda *X, dim=dim: X[1] if dim > 1 else 0 * X[0],
                "y": lambda *X, dim=dim: X[2] if dim > 2 else 0 * X[0],
                "z": lambda *X: 3 * X[0],
            }
            s = snapshot(hier, comp)
            r = deriv.curl(s)
            # curl = (dAz/dy - dAy/dz, dAx/dz - dAz/dx, dAy/dx - dAx/dy)
            expect = {"x": -1.0 if dim > 2 else 0.0, "y": -3.0, "z": -1.0 if dim > 1 else 0.0}
            for ld in r.levels:
                for c in "xyz":
                    np.testing.assert_allclose(ld.buffer(c), expect[c], atol=1e-10,
                                               err_msg=f"{dim}d {c}")

    def test_div_curl_roundoff(self):
        for name, hier in layouts(periodic_only=True):
            if not hier.periodic.all():
                continue
            dim = hier.dim
            k = 2 * np.pi / hier.domain

            def A(i, j, dim=dim, k=k):
                def f(*X):
                    Xs = list(X) + [0 * X[0]] * (3 - dim)
                    kk = list(k) + [0.0] * (3 - dim)
                    return np.sin(kk[i] * Xs[i] + 0.3) * np.cos(kk[j] * Xs[j] + 1.1)
                return f

            s = snapshot(hier, {"x": A(1, 2), "y": A(2, 0), "z": A(0, 1)})
            c = deriv.curl(s)
            dv = deriv.div(c)
            for cl, dl in zip(c.levels, dv.levels):
                geom = cl.geom
                ok = ~dilated_one_sided(geom)
                scale = max(np.abs(cl.buffer(n)).max() for n in "xyz") / geom.dl.min()
                err = np.abs(dl.buffer("value"))[ok].max() / scale
                self.assertLess(err, 1e-13, f"{name} loc{geom.level}")
                self.assertTrue(ok.any())

    def test_div_of_averaged_yee_field_is_zero(self):
        """Node-averaged discrete div-free Yee B (from A_z) has centred node div == 0."""
        for name, hier in [(n, h) for n, h in layouts() if h.dim == 2]:
            levels = []
            for geom in hier.levels:
                levels.append(yee_averaged_b(geom, seed=geom.level + 5))
            s = Snapshot(hier, tuple(levels))
            plain = deriv.div(s)
            masked = deriv.div(s, mask_faces=True)
            for ld, p, m in zip(s.levels, plain.levels, masked.levels):
                one = deriv.one_sided_mask(ld.geom)
                scale = np.abs(ld.buffer("x")).max() / ld.geom.dl.min()
                val = p.buffer("value")
                self.assertLess(np.abs(val[~one]).max() / scale, 1e-13, f"{name} loc{ld.level}")
                np.testing.assert_array_equal(np.isnan(m.buffer("value")), one)
                np.testing.assert_array_equal(m.buffer("value")[~one], val[~one])
                if one.any():  # the one-sided seam is visible (not ~0) on random data
                    self.assertGreater(np.abs(val[one]).max() / scale, 1e-3)

    def test_laplacian_vector_per_component(self):
        hier = make_amr(2, 2, seed=8)
        s = snapshot(hier, {"x": lambda x, y: x * x, "y": lambda x, y: y * y * 2,
                            "z": lambda x, y: x * y})
        r = deriv.laplacian(s)
        self.assertEqual(r.names, ("x", "y", "z"))
        for ld in r.levels:
            np.testing.assert_allclose(ld.buffer("x"), 2.0, atol=1e-9)
            np.testing.assert_allclose(ld.buffer("y"), 4.0, atol=1e-9)
            np.testing.assert_allclose(ld.buffer("z"), 0.0, atol=1e-9)

    def test_type_errors(self):
        hier = make_amr(2, 1, seed=1)
        sc = snapshot(hier, lambda x, y: x)
        vec = snapshot(hier, {"x": lambda x, y: x, "y": lambda x, y: y, "z": lambda x, y: x})
        ten = snapshot(hier, {"xx": lambda x, y: x, "yy": lambda x, y: y})
        for func, arg in ((deriv.curl, sc), (deriv.div, sc), (deriv.grad, vec),
                          (deriv.laplacian, ten), (deriv.grad, 3.0)):
            with self.assertRaises(TypeError) as ctx:
                func(arg)
            self.assertIn("expects", str(ctx.exception))
        with self.assertRaises(ValueError):
            deriv.ddx(sc.levels[0], "value", 3)
        with self.assertRaises(KeyError):
            deriv.ddx(sc.levels[0], "nope", 0)


def dilated_one_sided(geom, steps=2):
    """Nodes within ``steps`` nodes (box neighbourhood) of a one-sided node, by global index.

    The grid has N nodes on periodic axes (node N is node 0, rolls wrap) and N + 1 otherwise.
    """
    N = geom.domain_cells
    size = tuple(int(N[a]) if geom.periodic[a] else int(N[a]) + 1 for a in range(geom.dim))
    G = global_index(geom)
    Gk = tuple(G[a] % N[a] if geom.periodic[a] else G[a] for a in range(geom.dim))
    grid = np.zeros(size, dtype=bool)
    one = deriv.one_sided_mask(geom)
    grid[tuple(g[one] for g in Gk)] = True
    for _ in range(steps):
        for a in range(geom.dim):
            grown = grid.copy()
            for s in (-1, 1):
                r = np.roll(grid, s, axis=a)
                if not geom.periodic[a]:
                    edge = [slice(None)] * geom.dim
                    edge[a] = 0 if s == 1 else -1
                    r[tuple(edge)] = False
                grown |= r
            grid = grown
    return grid[Gk]


def yee_averaged_b(geom, seed):
    """2D: node B = 2-point average of the Yee B = curl(A_z e_z), A_z random per node."""
    N = geom.domain_cells
    rng = np.random.default_rng(seed)
    az = rng.standard_normal(tuple(N + 3))  # global nodes -1 .. N + 1

    def A(i, j):
        if geom.periodic[0]:
            i = i % N[0]
        if geom.periodic[1]:
            j = j % N[1]
        return az[i + 1, j + 1]

    G = global_index(geom)
    i, j = G[0], G[1]
    dx, dy = geom.dl
    # Yee: Bx(i, j+1/2) = (Az[i, j+1] - Az[i, j]) / dy ; By(i+1/2, j) = -(Az[i+1, j] - Az[i, j]) / dx
    bx = 0.5 * ((A(i, j + 1) - A(i, j)) / dy + (A(i, j) - A(i, j - 1)) / dy)
    by = -0.5 * ((A(i + 1, j) - A(i, j)) / dx + (A(i, j) - A(i - 1, j)) / dx)
    return LevelData.from_arrays(geom, 0.0, {"x": bx, "y": by, "z": np.zeros_like(bx)})


if __name__ == "__main__":
    unittest.main()
