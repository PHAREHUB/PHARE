"""The five reconstructions of PHARE's MHD solver, ported to numpy for the drawings of numerics/mhd_reconstruction.

Each function gives the value at the face between cell j and cell j+1, seen from the left (L) or from the right (R),
from the cell values u (an array; the caller keeps j-3 .. j+3 inside it). Formulas are copied from the C++ files named
in the comments (src/core/numerics/reconstructions/, pinned commit dc62fdd4e); the C++ names ul, u, ur are kept.
"""
import numpy as np

# ---- limiters: src/core/numerics/slope_limiters/{van_leer,min_mod}.hpp
def van_leer(a, b):                      # van_leer.hpp:10
    return 2.0 * a * b / (a + b) if a * b > 0.0 else 0.0


def minmod(*args):                       # min_mod.hpp:15-25: 0 unless all same sign, else smallest |.|
    if all(x > 0 for x in args) or all(x < 0 for x in args):
        return min(args, key=abs)
    return 0.0


# ---- constant.hpp:18-22
def constant_L(u, j): return u[j]
def constant_R(u, j): return u[j + 1]


# ---- linear.hpp:19-29, :47-67 (slope = limiter(u - ul, ur - u); face value = u +/- slope/2)
def linear_L(u, j, limiter):
    return u[j] + 0.5 * limiter(u[j] - u[j - 1], u[j + 1] - u[j])


def linear_R(u, j, limiter):
    return u[j + 1] - 0.5 * limiter(u[j + 1] - u[j], u[j + 2] - u[j + 1])


# ---- weno3.hpp:46-80 (eps = 1e-6, alpha = d / (beta + eps)^2)
def _weno3_weights(ul, u, ur, d0, d1, eps=1e-6):
    b0, b1 = (u - ul) ** 2, (ur - u) ** 2
    a0, a1 = d0 / (b0 + eps) ** 2, d1 / (b1 + eps) ** 2
    return a0 / (a0 + a1), a1 / (a0 + a1)


def _weno3_L(ul, u, ur):
    w0, w1 = _weno3_weights(ul, u, ur, 1 / 3, 2 / 3)
    return w0 * (-0.5 * ul + 1.5 * u) + w1 * (0.5 * u + 0.5 * ur)


def _weno3_R(ul, u, ur):
    w0, w1 = _weno3_weights(ul, u, ur, 2 / 3, 1 / 3)
    return w0 * (0.5 * u + 0.5 * ul) + w1 * (-0.5 * ur + 1.5 * u)


def weno3_L(u, j): return _weno3_L(u[j - 1], u[j], u[j + 1])
def weno3_R(u, j): return _weno3_R(u[j], u[j + 1], u[j + 2])


# ---- wenoz.hpp:58-107 (eps = 1e-40, tau5 = |beta0 - beta2|, alpha = d (1 + tau5 / (beta + eps)))
def _wenoz_weights(ull, ul, u, ur, urr, d0, d1, d2, eps=1e-40):
    b0 = 13 / 12 * (ull - 2 * ul + u) ** 2 + 1 / 4 * (ull - 4 * ul + 3 * u) ** 2
    b1 = 13 / 12 * (ul - 2 * u + ur) ** 2 + 1 / 4 * (ul - ur) ** 2
    b2 = 13 / 12 * (u - 2 * ur + urr) ** 2 + 1 / 4 * (3 * u - 4 * ur + urr) ** 2
    tau5 = abs(b0 - b2)
    a = [d0 * (1 + tau5 / (b0 + eps)), d1 * (1 + tau5 / (b1 + eps)), d2 * (1 + tau5 / (b2 + eps))]
    s = sum(a)
    return [x / s for x in a]


def _wenoz_L(ull, ul, u, ur, urr):
    w0, w1, w2 = _wenoz_weights(ull, ul, u, ur, urr, 1 / 10, 3 / 5, 3 / 10)
    return (w0 * (1 / 3 * ull - 7 / 6 * ul + 11 / 6 * u) + w1 * (-1 / 6 * ul + 5 / 6 * u + 1 / 3 * ur)
            + w2 * (1 / 3 * u + 5 / 6 * ur - 1 / 6 * urr))


def _wenoz_R(ull, ul, u, ur, urr):
    w0, w1, w2 = _wenoz_weights(ull, ul, u, ur, urr, 3 / 10, 3 / 5, 1 / 10)
    return (w0 * (1 / 3 * u + 5 / 6 * ul - 1 / 6 * ull) + w1 * (-1 / 6 * ur + 5 / 6 * u + 1 / 3 * ul)
            + w2 * (1 / 3 * urr - 7 / 6 * ur + 11 / 6 * u))


def wenoz_L(u, j): return _wenoz_L(u[j - 2], u[j - 1], u[j], u[j + 1], u[j + 2])
def wenoz_R(u, j): return _wenoz_R(u[j - 1], u[j], u[j + 1], u[j + 2], u[j + 3])


# ---- mp5.hpp:69-88 (alpha = 4; accepted unchanged if (f - u)(f - fMP) < 0, else clipped into [fmin, fmax])
def _mp5(v_m2, v_m1, u, v_p1, v_p2):
    alpha = 4.0
    f = (2 * v_m2 - 13 * v_m1 + 47 * u + 27 * v_p1 - 3 * v_p2) / 60
    Dil, Dir = u - v_m1, v_p1 - u
    fMP = u + minmod(Dir, alpha * Dil)
    fUL = u + alpha * Dil
    di, dir_ = Dir - Dil, (v_p2 - v_p1) - Dir
    d4i = minmod(4 * di - dir_, 4 * dir_ - di, di, dir_)
    fMD = (u + v_p1) / 2 - d4i / 2
    dil = Dil - (v_m1 - v_m2)
    d4im = minmod(4 * dil - di, 4 * di - dil, dil, di)
    fLC = u + Dil / 2 + 4 / 3 * d4im
    fmin = max(min(u, v_p1, fMD), min(u, fUL, fLC))
    fmax = min(max(u, v_p1, fMD), max(u, fUL, fLC))
    return f if (f - u) * (f - fMP) < 0.0 else min(max(f, fmin), fmax)


def mp5_L(u, j): return _mp5(u[j - 2], u[j - 1], u[j], u[j + 1], u[j + 2])
def mp5_R(u, j): return _mp5(u[j + 3], u[j + 2], u[j + 1], u[j], u[j - 1])   # mirrored call, mp5.hpp:64-66


METHODS = {
    "Constant": (constant_L, constant_R, 1),
    "Linear + MinMod": (lambda u, j: linear_L(u, j, lambda a, b: minmod(a, b)),
                        lambda u, j: linear_R(u, j, lambda a, b: minmod(a, b)), 2),
    "Linear + VanLeer": (lambda u, j: linear_L(u, j, van_leer), lambda u, j: linear_R(u, j, van_leer), 2),
    "WENO3": (weno3_L, weno3_R, 2),
    "WENOZ": (wenoz_L, wenoz_R, 3),
    "MP5": (mp5_L, mp5_R, 3),
}
