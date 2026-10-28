"""Numbers the physics report quotes. Pure numpy: tested without running PHARE."""

import numpy as np


def growth_rate_summary(gamma, theory, tol):
    err = abs(float(gamma) - theory)
    return {
        "gamma": float(gamma),
        "abs_error": err,
        "within_deck_tolerance": bool(err < tol),
    }


def shifted_profile_error(x, f0, f1, shift, length):
    x, f0, f1 = (np.asarray(a, dtype=float) for a in (x, f0, f1))
    expected = np.interp(np.mod(x - shift, length), x, f0, period=length)
    scale = np.mean(np.abs(f0 - f0.mean())) or 1.0
    return float(np.mean(np.abs(f1 - expected)) / scale)
