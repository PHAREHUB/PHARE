import numpy as np

from pyphare.core import box as boxm
from pyphare.core.gridlayout import GridLayout
from pyphare.core.phare_utilities import refinement_ratio
from pyphare.pharesee.hierarchy.patchdata import FieldData


# below is a drawing  representing a 1D coarse FieldData, its ghost box and the
# associated refined  FieldData this module aims at creating.
# o : coarse primal node
# + : fine primal node
# * : fine primal nodes created for simplicity of the refinement algorithm but
#     which are later discarded

#                              Field ghost box
#     ________________________________________________________________
#    |                                                                |
#    v                                                                v
#                                 Field box
#                 _________________________________________
#                |                                         |
#                v                                         v
# 0     1     2     3     4     5     6     7     8     9    10    11    12
#             |                                               |
# o     o     o     o     o     o     o     o     o     o     o     o     o
#             |                                               |
# *  *  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  *  *
#             |                                               |
# 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24
#
#


def fine_layout_from(field):
    return GridLayout(
        boxm.refine(field.box, refinement_ratio),
        field.origin,
        field.layout.dl / refinement_ratio,
        interp_order=field.layout.interp_order,
    )


def cropToFieldData(fine_data, field):
    # reminder: we were working on "fine_data" that covers the whole ghost box of the given coarse field
    # here we drop the parts of data that are beyond the refined patch ghost box.
    if field.box.ndim == 1:
        fine_data = fine_data[field.ghosts_nbr[0] : -field.ghosts_nbr[0]]
    elif field.box.ndim == 2:
        fine_data = fine_data[
            field.ghosts_nbr[0] : -field.ghosts_nbr[0],
            field.ghosts_nbr[1] : -field.ghosts_nbr[1],
        ]
    elif field.box.ndim == 3:
        fine_data = fine_data[
            field.ghosts_nbr[0] : -field.ghosts_nbr[0],
            field.ghosts_nbr[1] : -field.ghosts_nbr[1],
            field.ghosts_nbr[2] : -field.ghosts_nbr[2],
        ]
    fine_layout = fine_layout_from(field)
    return FieldData(fine_layout, field.field_name, data=fine_data)


def refine_electric(field, **kwargs):
    """
    This function refines the electric field from the given field
    See electric_field_refiner.hpp for more details on the refinement algo.
    """
    fine_data = makeDataset(field, **kwargs)
    coarse_data = kwargs.get("data", field.dataset[:])
    if "data" in kwargs:
        assert coarse_data.shape == field.dataset.shape
    primal_directions = field.primal_directions()

    if field.box.ndim == 1:
        # in 1D Ex, Ey and Ey cell edges are on top of coarse cell edges
        # the algo for electric refinement is that fine edges falling on top of coarse
        # ones share their value. So here we just take the coarse value whatever the
        # centering.
        print(field.field_name)
        print(field.box.shape)
        print(field.dataset.shape)
        fine_data[::2] = coarse_data
        fine_data[1::2] = coarse_data

    elif field.box.ndim == 2:
        if primal_directions[0] and not primal_directions[1]:
            # this is Ex
            # in 2D along X Ex fine edges are collocate with coarse edges
            # but along Y only even edges do, ones fall in  between 2 coarse edges
            # and take the average of their value

            # first do even Y indices
            # same for even and odd X fine indices
            fine_data[::2, ::2] = coarse_data[:, :]
            fine_data[1::2, ::2] = coarse_data[:, :]

            # now odd Y indices take the average of surrounding coarse edges
            # again, same for odd and even X fine indices
            fine_data[::2, 1::2] = 0.5 * (coarse_data[:, 1:] + coarse_data[:, :-1])
            fine_data[1::2, 1::2] = 0.5 * (coarse_data[:, 1:] + coarse_data[:, :-1])

    elif field.box.ndim == 3:
        raise NotImplementedError("3D electric field refinement not implemented")

    return cropToFieldData(fine_data, field)


def refine_magnetic(field, **kwargs):
    """
    This function refines the magnetic field from the given field
    See magnetic_field_refiner.hpp for more details on the refinement algo.
    """

    # reminder the "fine_data" we get covers the whole ghost box of the given coarse field
    fine_data = makeDataset(field, **kwargs)

    coarse_data = kwargs.get("data", field.dataset[:])
    if "data" in kwargs:
        assert coarse_data.shape == field.dataset.shape

    primal_directions = field.primal_directions()

    if field.box.ndim == 1:
        if primal_directions[0]:
            #  we're primal and 1D so this is Bx
            # even primal fine nodes fall on top of coarse ones so take the value
            # odd primal fine nodes fall in the middle of coarse ones so take the average
            fine_data[::2] = coarse_data
            fine_data[1::2] = 0.5 * (coarse_data[:-1] + coarse_data[1:])

        else:
            # we're 1D and dual so this is By or Bz
            # by flux conservation, fine dual nodes take the same value as the one
            # of the coarse cell they fall in
            # this means that the coarse values are  copied in both even and odd fine dual nodes
            # note that because the refinement ratio is  2 the number of fine cells is always even
            fine_data[::2] = coarse_data[:]
            fine_data[1::2] = fine_data[::2]

    elif field.box.ndim == 2:
        if primal_directions[0] and not primal_directions[1]:
            # this is Bx. even x faces fall on top of coarse x faces
            # the coarse face value is copied into even and odd fine faces
            fine_data[::2, ::2] = coarse_data[:, :]
            fine_data[::2, 1::2] = coarse_data[:, :]

            # odd x faces are in the middle of coarse x faces
            # and  take the average of their flux
            # as before the coarse face values are just copied into fine odd and even faces along y
            fine_data[1::2, ::2] = 0.5 * (coarse_data[:-1:, :] + coarse_data[1:, :])
            fine_data[1::2, 1::2] = fine_data[1::2, ::2]

        elif not primal_directions[0] and primal_directions[1]:
            # this is now By
            # even y indices fall on top of coarse faces so just take their value
            # copy the coarse value to odd and even x fine indices
            fine_data[::2, ::2] = coarse_data[:, :]
            fine_data[1::2, ::2] = fine_data[::2, ::2]

            # now odd y indices  are fine faces falling in between coarse faces
            # so we need to average the surrounding coarse faces in the y direction
            # as before we do that for even and odd x indices
            fine_data[::2, 1::2] = 0.5 * (coarse_data[:, 1:] + coarse_data[:, :-1])
            fine_data[1::2, 1::2] = fine_data[::2, 1::2]

        elif not all(primal_directions):
            # all directions are dual, this is Bz
            # Bz is easy since all 4 fine faces take the same value as the coarse one
            # and since we are in 2D the fine face is always colocated with the coarse one
            fine_data[::2, ::2] = coarse_data[:, :]
            fine_data[::2, 1::2] = fine_data[::2, ::2]
            fine_data[1::2, ::2] = fine_data[::2, ::2]
            fine_data[1::2, 1::2] = fine_data[::2, ::2]

    elif field.box.ndim == 3:
        raise NotImplementedError("3D magnetic field refinement not implemented")
    else:
        raise RuntimeError("impossible layout for a magnetic field")

    return cropToFieldData(fine_data, field)


def makeDataset(field, **kwargs):
    """
    this function returns the dataset for a refined field given a FieldData

    This is done by refining the coarse field given as argument over its whole
    ghost box and then crop it to obtain the final fine field.

    in the drawing above, this means that the  fine data is computed on the * and +
    nodes (here for primal)

    """
    assert isinstance(field, FieldData)
    kwarg_keys = ["data"]
    assert any([key in kwarg_keys for key in list(kwargs.keys())])

    refinement_ratio = 2
    primal_directions = field.primal_directions()

    # fine box has extra ghosts for padding against coarser such that at the normal number of ghosts are filled
    fine_box = boxm.shrink(
        boxm.refine(field.ghost_box, refinement_ratio), field.ghosts_nbr
    )
    fine_data = np.zeros(fine_box.shape + primal_directions + (field.ghosts_nbr * 2))

    return fine_data


def default_refine(field, **kwargs):
    """
    This function returns a FieldData refined from the given field, and which
    domain overlaps the domain of the given field.

    this refinement is the "default" refinement, not used for the electric or
    magnetic field.

    optional kwargs
      data : for overriding the field dataset so you don't have to make a temporary copy of a field just to have a different dataset
    """

    # reminder the "fine_data" we get covers the whole ghost box of the given coarse field
    fine_data = makeDataset(field, **kwargs)
    ghostX = field.ghosts_nbr[0]
    assert ghostX > 1

    coarse_data = kwargs.get("data", field.dataset[:])
    if "data" in kwargs:
        assert coarse_data.shape == field.dataset.shape

    cadence = 2
    assert cadence == refinement_ratio
    primal_directions = field.primal_directions()

    gX = 2  # level 0 ghost buffer from edge of normal coarse data box
    rgX = 4  # level 1 ghost buffer from edge of extra large fine data box

    if field.box.ndim == 1:
        if primal_directions[0]:
            # coarse primal on top of fine
            fine_data[rgX:-rgX:cadence] = coarse_data[gX:-gX]
            fine_data[rgX + 1 : -(rgX - 1) : cadence] = (
                0.5 * coarse_data[gX:-gX] + 0.5 * coarse_data[gX + 1 : -(gX - 1)]
            )
        else:
            fine_data[rgX : -(rgX + 1) : cadence] = (
                0.25 * coarse_data[gX - 1 : -(gX + 1)] + 0.75 * coarse_data[gX:-gX]
            )
            fine_data[rgX + 1 : -(rgX - 1) : cadence] = (
                0.25 * coarse_data[gX + 1 : -(gX - 1)] + 0.75 * coarse_data[gX:-gX]
            )

    if field.box.ndim == 2:
        assert field.ghosts_nbr[1] > 1
        gY = 2
        rgY = 4
        cad = cadence

        if all(primal_directions):
            fine_data[rgX:-rgX:cad, rgY:-rgY:cad] = coarse_data[gX:-gX, gY:-gY]
            fine_data[rgX + 1 : -(rgX - 1) : cad, rgY:-rgY:cad] = (
                0.5 * coarse_data[gX:-gX, gY:-gY]
                + 0.5 * coarse_data[gX + 1 : -(gX - 1), gY:-gY]
            )
            fine_data[rgX:-rgX:cad, rgY + 1 : -(rgY - 1) : cad] = (
                0.5 * coarse_data[gX:-gX, gY:-gY]
                + 0.5 * coarse_data[gX:-gX, gY + 1 : -(gY - 1)]
            )
            fine_data[rgX + 1 : -(rgX - 1) : cad, rgY + 1 : -(rgY - 1) : cad] = (
                0.25 * coarse_data[gX:-gX, gY:-gY]
                + 0.25 * coarse_data[gX + 1 : -(gX - 1), gY:-gY]
                + 0.25 * coarse_data[gX:-gX, gY + 1 : -(gY - 1)]
                + 0.25 * coarse_data[gX + 1 : -(gX - 1), gY + 1 : -(gY - 1)]
            )

        elif primal_directions[0] and not primal_directions[1]:
            fine_data[rgX:-rgX:cad, rgY:-rgY:cad] = (
                0.25 * coarse_data[gX:-gX, gY - 1 : -(gY + 1)]
                + 0.75 * coarse_data[gX:-gX, gY:-gY]
            )
            fine_data[rgX + 1 : -(rgX - 1) : cad, rgY:-rgY:cad] = 0.75 * (
                0.5 * coarse_data[gX + 1 : -(gX - 1), gY:-gY]
                + 0.5 * coarse_data[gX:-gX, gY:-gY]
            ) + 0.25 * (
                0.5 * coarse_data[gX + 1 : -(gX - 1), gY - 1 : -(gY + 1)]
                + 0.5 * coarse_data[gX:-gX, gY - 1 : -(gY + 1)]
            )
            fine_data[rgX:-rgX:cad, rgY + 1 : -(rgY - 1) : cad] = (
                0.25 * coarse_data[gX:-gX, gY + 1 : -(gY - 1)]
                + 0.75 * coarse_data[gX:-gX, gY:-gY]
            )
            fine_data[rgX + 1 : -(rgX - 1) : cad, rgY + 1 : -(rgY - 1) : cad] = 0.75 * (
                0.5 * coarse_data[gX + 1 : -(gX - 1), gY:-gY]
                + 0.5 * coarse_data[gX:-gX, gY:-gY]
            ) + 0.25 * (
                0.5 * coarse_data[gX + 1 : -(gX - 1), gY + 1 : -(gY - 1)]
                + 0.5 * coarse_data[gX:-gX, gY + 1 : -(gY - 1)]
            )

        elif not primal_directions[0] and primal_directions[1]:
            fine_data[rgX:-rgX:cad, rgY:-rgY:cad] = (
                0.25 * coarse_data[gX - 1 : -(gX + 1), gY:-gY]
                + 0.75 * coarse_data[gX:-gX, gY:-gY]
            )
            fine_data[rgX + 1 : -(rgX - 1) : cad, rgY:-rgY:cad] = (
                0.25 * coarse_data[gX + 1 : -(gX - 1), gY:-gY]
                + 0.75 * coarse_data[gX:-gX, gY:-gY]
            )
            fine_data[rgX:-rgX:cad, rgY + 1 : -(rgY - 1) : cad] = 0.5 * (
                0.25 * coarse_data[gX - 1 : -(gX + 1), gY:-gY]
                + 0.75 * coarse_data[gX:-gX, gY:-gY]
            ) + 0.5 * (
                0.25 * coarse_data[gX - 1 : -(gX + 1), gY + 1 : -(gY - 1)]
                + 0.75 * coarse_data[gX:-gX, gY + 1 : -(gY - 1)]
            )
            fine_data[rgX + 1 : -(rgX - 1) : cad, rgY + 1 : -(rgY - 1) : cad] = 0.5 * (
                0.25 * coarse_data[gX + 1 : -(gX - 1), gY:-gY]
                + 0.75 * coarse_data[gX:-gX, gY:-gY]
            ) + 0.5 * (
                0.25 * coarse_data[gX + 1 : -(gX - 1), gY + 1 : -(gY - 1)]
                + 0.75 * coarse_data[gX:-gX, gY + 1 : -(gY - 1)]
            )

        elif not any(primal_directions):
            fine_data[rgX:-rgX:cad, rgY:-rgY:cad] = 0.75 * (
                0.25 * coarse_data[gX - 1 : -(gX + 1), gY:-gY]
                + 0.75 * coarse_data[gX:-gX, gY:-gY]
            ) + 0.25 * (
                0.25 * coarse_data[gX - 1 : -(gX + 1), gY - 1 : -(gY + 1)]
                + 0.75 * coarse_data[gX:-gX, gY - 1 : -(gY + 1)]
            )
            fine_data[rgX + 1 : -(rgX - 1) : cad, rgY:-rgY:cad] = 0.75 * (
                0.25 * coarse_data[gX + 1 : -(gX - 1), gY:-gY]
                + 0.75 * coarse_data[gX:-gX, gY:-gY]
            ) + 0.25 * (
                0.25 * coarse_data[gX + 1 : -(gX - 1), gY - 1 : -(gY + 1)]
                + 0.75 * coarse_data[gX:-gX, gY - 1 : -(gY + 1)]
            )
            fine_data[rgX:-rgX:cad, rgY + 1 : -(rgY - 1) : cad] = 0.75 * (
                0.25 * coarse_data[gX - 1 : -(gX + 1), gY:-gY]
                + 0.75 * coarse_data[gX:-gX, gY:-gY]
            ) + 0.25 * (
                0.25 * coarse_data[gX - 1 : -(gX + 1), gY + 1 : -(gY - 1)]
                + 0.75 * coarse_data[gX:-gX, gY + 1 : -(gY - 1)]
            )
            fine_data[rgX + 1 : -(rgX - 1) : cad, rgY + 1 : -(rgY - 1) : cad] = 0.75 * (
                0.25 * coarse_data[gX + 1 : -(gX - 1), gY:-gY]
                + 0.75 * coarse_data[gX:-gX, gY:-gY]
            ) + 0.25 * (
                0.25 * coarse_data[gX + 1 : -(gX - 1), gY + 1 : -(gY - 1)]
                + 0.75 * coarse_data[gX:-gX, gY + 1 : -(gY - 1)]
            )

    return cropToFieldData(fine_data, field)


def refine_time_interpolate(
    datahier,
    quantities,
    coarse_ilvl,
    coarsest_time_before,
    coarsest_time_after,
    fine_subcycle_times,
):
    """
    returns {qty : { subcycle_time: [refined_time_interpolated_fields]}}
    """

    from tests.core.numerics.interpolator.interpolator_test import time_interpolate

    def _sort(patches):
        return sorted(patches, key=lambda p: p.origin.all())

    interpolated_fields = {qty: {} for qty in quantities}

    coarse_before_patches = _sort(
        datahier.level(coarse_ilvl, coarsest_time_before).patches
    )
    coarse_after_patches = _sort(
        datahier.level(coarse_ilvl, coarsest_time_after).patches
    )
    assert len(coarse_before_patches) == len(coarse_after_patches)

    for qty in quantities:
        if qty[0] == "B":
            refine_algo = refine_magnetic
        elif qty[0] == "E":
            refine_algo = refine_electric
        else:
            refine_algo = default_refine

        for fine_subcycle_time in fine_subcycle_times:
            interpolated_fields[qty][fine_subcycle_time] = []
            for coarsePatch_idx in range(len(coarse_before_patches)):
                coarse_before_patch = coarse_before_patches[coarsePatch_idx]
                coarse_after_patch = coarse_after_patches[coarsePatch_idx]
                assert coarse_before_patch.box == coarse_after_patch.box
                coarseBefore_pd = coarse_before_patch.patch_datas[qty]
                coarseAfter_pd = coarse_after_patch.patch_datas[qty]
                interpolated_fields[qty][fine_subcycle_time] += [
                    refine_algo(
                        coarseBefore_pd,
                        data=time_interpolate(
                            coarsest_time_before,
                            coarsest_time_after,
                            fine_subcycle_time,
                            coarseBefore_pd.dataset[:],
                            coarseAfter_pd.dataset[:],
                        ),
                    )
                ]

    return interpolated_fields
