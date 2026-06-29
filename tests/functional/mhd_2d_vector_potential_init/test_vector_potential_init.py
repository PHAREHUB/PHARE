#!/usr/bin/env python3
"""
Vector-potential init of the MHD magnetic field (2D): end-to-end smoke test.

A field built as B = curl(A_z z_hat) with the discrete curl is divergence-free in the discrete
(Yee) sense by construction (div . curl = 0). The *discrete* div-B property is proven directly
by the C++ unit test tests/core/numerics/mhd_vector_potential_init; the dumped diagnostics cannot
measure the face-staggered divergence (the MHD total-B reconstruction is cell-centred / node
interpolated). So this functional test only checks that the potential init (a0z/a1z) drives a
real simulation end-to-end and produces a finite, bounded state, and that the legacy component
init still runs.

Run:  phare-run python test_vector_potential_init.py
"""

import os
import subprocess
import sys

import numpy as np
from pyphare.pharesee.run import Run

HERE = os.path.dirname(os.path.abspath(__file__))
V_BOUND = 1e-2  # rest-ish state: |V| must stay small over the short run (sanity, not machine zero)


def _run_case(b_init, ncells="64"):
    env = dict(os.environ)
    env.update({"PHARE_B_INIT": b_init, "PHARE_NCELLS": ncells})
    subprocess.run([sys.executable, "case.py"], cwd=HERE, env=env, check=True)


def _max_abs_v(diag_dir):
    # Read mhd_V at its native centring (no interpolation, which would pull in the NaN-filled
    # boundary ghosts of the dump). Measure the interior (domain) cells only.
    run = Run(os.path.join(HERE, diag_dir))
    m = 0.0
    finite = True
    for t in run.times("mhd_V"):
        h = run._get_hierarchy(float(t), "mhd_V.h5")
        for patch in h.level(0, float(t)).patches:
            for pd in patch.patch_datas.values():
                g = pd.ghosts_nbr
                a = pd.dataset[:]
                sl = tuple(slice(int(gi), -int(gi) if int(gi) else None) for gi in g)
                interior = a[sl]
                finite &= bool(np.all(np.isfinite(interior)))
                m = max(m, float(np.max(np.abs(interior))))
    return m, finite


def main():
    _run_case("potential")
    _run_case("components")

    ok = True
    for b_init in ("potential", "components"):
        mv, finite = _max_abs_v(f"phare_outputs_{b_init}_n64_WENOZ")
        passed = finite and np.isfinite(mv) and mv < V_BOUND
        print(
            f"[{'PASS' if passed else 'FAIL'}] {b_init} init: finite={finite}, "
            f"max|V| = {mv:.3e} (< {V_BOUND:.0e})"
        )
        ok &= passed

    if not ok:
        sys.exit(1)
    print("vector-potential init end-to-end OK")


if __name__ == "__main__":
    main()
