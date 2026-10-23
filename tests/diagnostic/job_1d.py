#!/usr/bin/env python3

import os

import pyphare.pharein as ph
from tests.diagnostic import dump_all_diags
from tests.simulator import basicSimulatorArgs, makeBasicModel

out = os.environ["PHARE_DIAGS_TEST_DIR"]  # per test binary, see test-diagnostics.cpp
simInput = {
    "diag_options": {"format": "phareh5", "options": {"dir": out, "mode": "overwrite"}}
}

ph.Simulation(**basicSimulatorArgs(ndim=1, interp=1, **simInput))
model = makeBasicModel()
ph.ElectronModel(closure="isothermal", Te=0.12)
dump_all_diags(model.populations)
