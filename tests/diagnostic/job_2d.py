#!/usr/bin/env python3

import pyphare.pharein as ph
from tests.diagnostic import dump_all_diags
from tests.simulator import basicSimulatorArgs, makeBasicModel

from pyphare import cpp  # concurrent serial vs parallel tests

if ph.PHARE_EXE:
    out = f"phare_outputs/diags_2d/{cpp.mpi_size()}"
    simInput = {
        "diag_options": {
            "format": "phareh5",
            "options": {"dir": out, "mode": "overwrite"},
        }
    }

    ph.Simulation(**basicSimulatorArgs(ndim=2, interp=1, **simInput))
    model = makeBasicModel()
    ph.ElectronModel(closure="isothermal", Te=0.12)
    dump_all_diags(model.populations)
