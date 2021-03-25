import numpy as np

from pyphare.core import box as boxm
from pyphare.core.gridlayout import GridLayout
from pyphare.pharesee.hierarchy import FieldData


# only supports ndim == 1, updates needed for > 1
# assumes refinement_ratio is 2

def refine(field, **kwargs):
    """
      optional kwargs
        data : for overriding the field dataset so you don't have to make a temporary copy of a field just to have a different dataset
    """
    assert isinstance(field, FieldData)
    kwarg_keys = ["data"]
    assert any([key in kwarg_keys for key in list(kwargs.keys())])

    refinement_ratio = 2
    primal_directions = field.primal_directions()

    data = kwargs.get("data", field.dataset[:])
    if "data" in kwargs:
        assert data.shape == field.dataset.shape

    fine_box    = boxm.refine(field.box, refinement_ratio)
    fine_layout = GridLayout(fine_box, field.origin, field.layout.dl/refinement_ratio, interp_order=field.layout.interp_order)

    assert field.box.ndim == 1
    fine_data = np.zeros(fine_box.shape + primal_directions + (field.ghosts_nbr[0] * 2))

    ghostX = field.ghosts_nbr[0]
    assert ghostX > 0

    cadence = 2

    if primal_directions[0]:
        fine_data[ghostX:-ghostX:cadence] = data[ghostX:-ghostX] # coarse primal on top of fine
        fine_data[ghostX + 1:-(ghostX - 1):cadence] = 0.5*data[ghostX:-ghostX] + 0.5*data[ghostX + 1:-(ghostX - 1)]
    else:
        fine_data[ghostX:-(ghostX + 1):cadence] = 0.25*data[4:-(ghostX + 1)] + 0.75*data[ghostX:-ghostX]
        fine_data[ghostX + 1:-(ghostX - 1):cadence] = 0.25*data[ghostX + 1:-(ghostX - 1)] + 0.75*data[ghostX:-ghostX]

    return FieldData(fine_layout, field.field_name, data=fine_data)
