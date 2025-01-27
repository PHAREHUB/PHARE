#! /usr/bin/env python
# -*- coding: utf-8 -*-


import h5py
import sys
import numpy as np


def toFlatPrimal(ph_bx, ph_by, ph_bz, npx, npy, npz, gn=2):

    nbrPoints = npx * npy * npz
    b = np.zeros((nbrPoints, 3), dtype="f")

    # pure primal arrays
    bx = np.zeros((npx, npy, npz), dtype=np.float32)
    by = np.zeros((npx, npy, npz), dtype=np.float32)
    bz = np.zeros((npx, npy, npz), dtype=np.float32)

    # converts yee to pure primal
    # by averaging in the dual direction
    # first and last will be fucked
    bx[:, :, 0] = (
        ph_bx[gn:-gn, gn - 1 : -gn + 1][:, 1:] + ph_bx[gn:-gn, gn - 1 : -gn + 1][:, :-1]
    ) * 0.5
    by[:, :, 0] = (
        ph_by[gn - 1 : -gn + 1, gn:-gn][1:, :] + ph_by[gn - 1 : -gn + 1, gn:-gn][:-1, :]
    ) * 0.5
    bz[:, :, 0] = (
        ph_bz[gn - 1 : -gn + 1, gn - 1 : -gn + 1][1:, 1:]
        + ph_bz[gn - 1 : -gn + 1, gn - 1 : -gn + 1][:-1, :-1]
    ) * 0.5

    bx[:, :, 1] = bx[:, :, 0]
    by[:, :, 1] = by[:, :, 0]
    bz[:, :, 1] = bz[:, :, 0]

    b[:, 0] = bx.flatten(order="F")
    b[:, 1] = by.flatten(order="F")
    b[:, 2] = bz.flatten(order="F")

    return b


def boxFromPatch(patch):
    lower = patch.attrs["lower"]
    upper = patch.attrs["upper"]
    return [lower[0], upper[0], lower[1], upper[1], 0, 0]  # 2D


def nbrNodes(box):
    lower = box[0], box[2], box[4]
    upper = box[1], box[3], box[5]
    npx = upper[0] - lower[0] + 1 + 1
    npy = upper[1] - lower[1] + 1 + 1
    npz = upper[2] - lower[2] + 1 + 1
    return npx, npy, npz


def main():

    path = sys.argv[1]
    phare_h5 = h5py.File(f"{path}/EM_B.h5", "r")
    times_str = list(phare_h5["t"].keys())
    times = np.asarray([float(time) for time in times_str])
    times.sort()
    root_spacing = phare_h5.attrs["cell_width"]

    max_nbr_level = 0
    for time in times_str:
        nbrLevels = len(phare_h5["t"][time].keys())
        if max_nbr_level < nbrLevels:
            max_nbr_level = nbrLevels

    vtk = h5py.File(f"{path}/EM_B.vtkhdf", "w")
    root = vtk.create_group("VTKHDF", track_order=True)
    root.attrs["Version"] = (2, 2)
    type = b"OverlappingAMR"
    root.attrs.create("Type", type, dtype=h5py.string_dtype("ascii", len(type)))
    description = b"XYZ"
    root.attrs.create(
        "GridDescription",
        description,
        dtype=h5py.string_dtype("ascii", len(description)),
    )
    origin = [0, 0, 0]
    root.attrs.create("Origin", origin, dtype="f")

    numberOfTimes = 3  # times.size

    steps = root.create_group("Steps")
    steps.attrs.create("NSteps", data=numberOfTimes, dtype="i8")
    steps.create_dataset("Values", data=times[:numberOfTimes])

    print(f"maxlevel:{max_nbr_level}")

    for ilvl in range(max_nbr_level)[:]:

        print(f"Processing level {ilvl}")

        lvl = root.create_group(f"Level{ilvl}")
        lvl_spacing = root_spacing
        lvl_spacing = [dl / 2**ilvl for dl in lvl_spacing] + [0.0]
        lvl.attrs.create("Spacing", lvl_spacing, dtype="f")
        steps_lvl = steps.create_group(f"Level{ilvl}")

        AMRBox = []
        step_nbrBoxes = []
        AMRBoxOffsets = []
        dataOffsets = []

        cellData_g = lvl.create_group("CellData")
        pointData_g = lvl.create_group("PointData")
        fieldData_g = lvl.create_group("FieldData")
        cellDataOffset_g = steps_lvl.create_group("CellDataOffset")
        pointDataOffset_g = steps_lvl.create_group("PointDataOffset")
        FieldDataOffset_g = steps_lvl.create_group("FieldDataOffset")

        first = True

        # these are values for the first time
        # they will be reset for other times
        current_size = 0
        amr_box_offset = 0

        for time in times[:numberOfTimes]:
            nbr_boxes = 0
            print(f"time {time}")
            time_str = f"{time:.10f}"
            phare_lvl_name = f"pl{ilvl}"
            dataOffsets += [current_size]
            print("amr_box_offset", amr_box_offset)
            AMRBoxOffsets += [amr_box_offset]
            # pass this time if this level does not exist
            if phare_lvl_name not in phare_h5["t"][time_str].keys():
                print(f"no level {ilvl} at time {time}")
                continue
            phare_lvl = phare_h5["t"][time_str][phare_lvl_name]

            for patch_id in list(phare_lvl.keys())[:]:
                # print(f"patch {patch_id}")
                patch = phare_lvl[patch_id]

                ph_bx = patch["EM_B_x"][:]
                ph_by = patch["EM_B_y"][:]
                ph_bz = patch["EM_B_z"][:]

                box = boxFromPatch(patch)
                AMRBox.append(box)
                nbr_boxes += 1
                npx, npy, npz = nbrNodes(box)
                b = toFlatPrimal(ph_bx, ph_by, ph_bz, npx, npy, npz)
                # b = np.ones(12 * 14 * 2)  # np.ones((12 * 14, 3), dtype="f")
                # print(f"b shape :{b.shape} npx:{npx} npy:{npy}")
                # print(f"box:{box}")

                if first:
                    # this is the first patch of the first time
                    # for this level, we need to create the dataset
                    # the first dimension is the total # of points
                    # which is unknown hence the None for maxshape
                    # print(f"b shape :{b.shape[0]}")
                    pointData_b = pointData_g.create_dataset(
                        "B", data=b, maxshape=(None, 3)
                    )
                    first = False

                else:
                    # print("else")
                    # print(f"b shape :{b.shape[0]}")
                    # dataset already created with shape (current_size,3)
                    # we add b.shape[0] points (=npx*npy) to the first dim
                    # hence need to resize the dataset.
                    pointData_b.resize(current_size + b.shape[0], axis=0)
                    pointData_b[current_size:, :] = b
                # pass

                current_size += b.shape[0]
            # print(f"current_size:{current_size}")

            # of of the patch loops at that time

            # number of boxes added for ilvl across all times
            amr_box_offset += nbr_boxes
            step_nbrBoxes += [nbr_boxes]

        lvl.create_dataset("AMRBox", data=AMRBox)
        steps_lvl.create_dataset("AMRBoxOffset", data=AMRBoxOffsets)
        steps_lvl.create_dataset("NumberOfAMRBox", data=step_nbrBoxes)
        pointDataOffset_g.create_dataset("B", data=dataOffsets)

    vtk.close()


if __name__ == "__main__":
    main()
