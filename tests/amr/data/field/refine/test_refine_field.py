import numpy as np

from pyphare.core import box as boxm
from pyphare.core.gridlayout import GridLayout
from pyphare.pharesee.hierarchy import FieldData


# only supports ndim == 1, updates needed for > 1
# assumes refinement_ratio is 2

def refine(field, **kwargs):
    assert isinstance(field, FieldData)
    kwarg_keys = ["data"] # doing if data == None for input data=None, fails if it's an np array
    assert any([key in kwarg_keys for key in list(kwargs.keys())])

    refinement_ratio = 2
    primal_directions = field.primal_directions()
    data   = kwargs.get("data", field.dataset[:])
    ghosts = field.ghosts_nbr

    fine_box    = boxm.refine(field.box, refinement_ratio)
    fine_layout = GridLayout(fine_box, field.origin, field.layout.dl/refinement_ratio, interp_order=field.layout.interp_order)

    assert field.box.ndim == 1

    fine_data = np.zeros(fine_box.shape + primal_directions + (field.ghosts_nbr[0] * 2))

    if primal_directions[0]:
        fine_data[5:-5:2] = data[5:-5] # coarse primal on top of fine
        fine_data[6:-4:2] = 0.5*data[5:-5] + 0.5*data[6:-4]
    else:
        fine_data[5:-6:2] = 0.25*data[4:-6] + 0.75*data[5:-5]
        fine_data[6:-4:2] = 0.25*data[6:-4] + 0.75*data[5:-5]

    return FieldData(fine_layout, field.field_name, data=fine_data)
