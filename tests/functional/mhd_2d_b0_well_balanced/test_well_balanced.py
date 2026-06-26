#!/usr/bin/env python3
"""
Regression test: the MHD scheme must be well-balanced w.r.t. the background field B0.

Exact solution of every case below is REST. A well-balanced scheme keeps max|V| at
machine zero. A non-well-balanced scheme develops a spurious velocity from the B0
gradient (the bug that bends magnetospheric field lines upstream of the bow shock).

Cases:
  - dipole B0, B1 = 0        : B0-self Maxwell-stress well-balancing
  - dipole B0, B1 = uniform  : cross-term well-balancing, IMF-in-B1 convention

Run:  phare-run python test_well_balanced.py
"""

import os
import subprocess
import sys

import numpy as np
from pyphare.pharesee.run import Run

HERE = os.path.dirname(os.path.abspath(__file__))
TOL = 1e-12  # machine-zero: the rest state must be preserved to round-off


def _run_case(env_overrides):
    env = dict(os.environ)
    env.update(env_overrides)
    subprocess.run([sys.executable, "case.py"], cwd=HERE, env=env, check=True)


def _max_abs_v(diag_dir):
    run = Run(os.path.join(HERE, diag_dir))
    m = 0.0
    for t in run.times("mhd_V"):
        h = run._get_hier_for(float(t), "mhd_V")
        for patch in h.level(0, float(t)).patches:
            for pd in patch.patch_datas.values():
                m = max(m, float(np.max(np.abs(pd.dataset[:]))))
    return m


def check(name, env, diag_dir):
    _run_case(env)
    mv = _max_abs_v(diag_dir)
    ok = mv < TOL
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: max|V| = {mv:.3e} (tol {TOL:.0e})")
    return ok


def main():
    n = "64"  # small + fast for a regression gate
    ok = True
    ok &= check(
        "dipole B0, B1=0",
        {"PHARE_B0_MODE": "dipole", "PHARE_B1_MODE": "zero", "PHARE_NCELLS": n},
        f"phare_outputs_dipole_b1zero_n{n}_WENOZ",
    )
    ok &= check(
        "dipole B0, B1=uniform",
        {"PHARE_B0_MODE": "dipole", "PHARE_B1_MODE": "uniform", "PHARE_NCELLS": n},
        f"phare_outputs_dipole_b1uniform_n{n}_WENOZ",
    )
    if not ok:
        sys.exit(1)
    print("all well-balanced checks passed")


if __name__ == "__main__":
    main()
