#! /usr/bin/env python
# -*- coding: utf-8 -*-


import h5py
import sys
import numpy as np
import os


def ndim_from(npx, npy, npz):
    if npx > 2 and npy > 2 and npz > 2:
        return 3
    elif npx > 2 and npy > 2 and npz == 2:
        return 2
    elif npx > 2 and npy == 2 and npz == 2:
        return 1
    else:
        raise ValueError(
            f" cannot infer dimension from (npx, npy, npz) = {npx} {npy} {npz}"
        )


def BtoFlatPrimal(ph_bx, ph_by, ph_bz, npx, npy, npz, gn=2):

    nbrPoints = npx * npy * npz
    b = np.zeros((nbrPoints, 3), dtype="f")

    # pure primal arrays
    bx = np.zeros((npx, npy, npz), dtype=np.float32)
    by = np.zeros((npx, npy, npz), dtype=np.float32)
    bz = np.zeros((npx, npy, npz), dtype=np.float32)

    domainP1 = slice(gn - 1, -gn + 1)
    domain = slice(gn, -gn)

    # converts yee to pure primal
    # we average in the dual direction so we need one extract ghost data in that direction

    if ndim_from(npx, npy, npz) == 2:
        bx[:, :, 0] = (
            ph_bx[gn:-gn, gn - 1 : -gn + 1][:, 1:]
            + ph_bx[gn:-gn, gn - 1 : -gn + 1][:, :-1]
        ) * 0.5
        by[:, :, 0] = (
            ph_by[gn - 1 : -gn + 1, gn:-gn][1:, :]
            + ph_by[gn - 1 : -gn + 1, gn:-gn][:-1, :]
        ) * 0.5

        bz[:, :, 0] = 0.5 * (
            0.5 * (ph_bz[domainP1, domain][1:, :] + ph_bz[domainP1, domain][:-1, :])
            + 0.5 * (ph_bz[domain, domainP1][:, 1:] + ph_bz[domain, domainP1][:, :-1])
        )

        bx[:, :, 1] = bx[:, :, 0]
        by[:, :, 1] = by[:, :, 0]
        bz[:, :, 1] = bz[:, :, 0]

    elif ndim_from(npx, npy, npz) == 3:
        # Bx is (primal, dual, dual)
        print("bx", ph_bx.shape)
        print(
            "bx[domain, domainP1, domain][:,1:,:]",
            ph_bx[domain, domainP1, domain][:, 1:, :].shape,
        )
        print(
            "bx[domain, domainP1, domain][:,:-1,:]",
            ph_bx[domain, domainP1, domain][:, :-1, :].shape,
        )
        bx[:, :, :] = 0.5 * (
            0.5
            * (
                ph_bx[domain, domainP1, domain][:, 1:, :]
                + ph_bx[domain, domainP1, domain][:, :-1, :]
            )
            + 0.5
            * (
                ph_bx[domain, domain, domainP1][:, :, 1:]
                + ph_bx[domain, domain, domainP1][:, :, :-1]
            )
        )
        # By is (dual, primal, dual)
        by[:, :, :] = 0.5 * (
            0.5
            * (
                ph_by[domain, domainP1, domain][:, 1:, :]
                + ph_by[domain, domainP1, domain][:, :-1, :]
            )
            + 0.5
            * (
                ph_bz[domain, domain, domainP1][:, :, 1:]
                + ph_bz[domain, domain, domainP1][:, :, :-1]
            )
        )
        # Bz is (dual, dual, primal)
        bz[:, :, :] = 0.5 * (
            0.5
            * (
                ph_bz[domain, domainP1, domain][:, 1:, :]
                + ph_bz[domain, domainP1, domain][:, :-1, :]
            )
            + 0.5
            * (
                ph_bz[domainP1, domain, domain][1:, :, :]
                + ph_bz[domainP1, domain, domain][:-1, :, :]
            )
        )

    b[:, 0] = bx.flatten(order="F")
    b[:, 1] = by.flatten(order="F")
    b[:, 2] = bz.flatten(order="F")

    return b


def EtoFlatPrimal(ph_ex, ph_ey, ph_ez, npx, npy, npz, gn=2):

    nbrPoints = npx * npy * npz
    e = np.zeros((nbrPoints, 3), dtype="f")

    # pure primal arrays
    ex = np.zeros((npx, npy, npz), dtype=np.float32)
    ey = np.zeros((npx, npy, npz), dtype=np.float32)
    ez = np.zeros((npx, npy, npz), dtype=np.float32)

    domainP1 = slice(gn - 1, -gn + 1)
    domain = slice(gn, -gn)

    # converts yee to pure primal
    # we average in the dual direction so we need one extract ghost data in that direction
    if ndim_from(npx, npy, npz) == 2:
        ex[:, :, 0] = (
            ph_ex[domainP1, domain][1:, :] + ph_ex[domainP1, domain][:-1, :]
        ) * 0.5

        ey[:, :, 0] = (
            ph_ey[domain, domainP1][:, 1:] + ph_ey[domain, domainP1][:, :-1]
        ) * 0.5

        # ez already primal in 2D
        ez[:, :, 0] = ph_ez[domain, domain][:, :]

        ex[:, :, 1] = ex[:, :, 0]
        ey[:, :, 1] = ey[:, :, 0]
        ez[:, :, 1] = ez[:, :, 0]

    elif ndim_from(npx, npy, npz) == 3:
        ex[:, :, :] = (
            ph_ex[domainP1, domain, domain][1:, :, :]
            + ph_ex[domainP1, domain, domain][:-1, :, :]
        ) * 0.5

        ey[:, :, :] = (
            ph_ey[domain, domainP1, domain][:, 1:, :]
            + ph_ey[domain, domainP1, domain][:, :-1, :]
        ) * 0.5

        ez[:, :, :] = (
            ph_ez[domain, domain, domainP1][:, :, 1:]
            + ph_ez[domain, domain, domainP1][:, :, :-1]
        ) * 0.5

    e[:, 0] = ex.flatten(order="F")
    e[:, 1] = ey.flatten(order="F")
    e[:, 2] = ez.flatten(order="F")

    return e


def primalScalarToFlatPrimal(ph_scalar, npx, npy, npz, gn=2):

    scalar3d = np.zeros((npx, npy, npz), dtype="f")
    if ndim_from(npx, npy, npz) == 2:
        domain = slice(gn, -gn)
        scalar3d[:, :, 0] = ph_scalar[domain, domain]
        scalar3d[:, :, 1] = ph_scalar[domain, domain]
        return scalar3d.flatten(order="F")
    elif ndim_from(npx, npy, npz) == 3:
        domain = slice(gn, -gn)
        scalar3d[:, :, :] = ph_scalar[domain, domain, domain]
        return scalar3d.flatten(order="F")


def primalVectorToFlatPrimal(ph_vx, ph_vy, ph_vz, npx, npy, npz, gn=2):
    nbrPoints = npx * npy * npz
    v = np.zeros((nbrPoints, 3), dtype="f")

    vx = primalScalarToFlatPrimal(ph_vx, npx, npy, npz, gn)
    vy = primalScalarToFlatPrimal(ph_vy, npx, npy, npz, gn)
    vz = primalScalarToFlatPrimal(ph_vz, npx, npy, npz, gn)

    v[:, 0] = vx.flatten(order="F")
    v[:, 1] = vy.flatten(order="F")
    v[:, 2] = vz.flatten(order="F")

    return v


def boxFromPatch(patch):
    lower = patch.attrs["lower"]
    upper = patch.attrs["upper"]
    if len(lower) == 3:
        return [lower[0], upper[0], lower[1], upper[1], lower[2], upper[2]]
    elif len(lower) == 2:
        return [lower[0], upper[0], lower[1], upper[1], 0, 0]
    elif len(lower) == 1:
        return [lower[0], upper[0], 0, 0, 0, 0]
    else:
        raise ValueError(f"Unknown dimension {len(lower)}")


def nbrNodes(box):
    lower = box[0], box[2], box[4]
    upper = box[1], box[3], box[5]

    # +1 for nbr cells from lo to up, and +1 for primal
    npx = upper[0] - lower[0] + 1 + 1
    npy = upper[1] - lower[1] + 1 + 1
    npz = upper[2] - lower[2] + 1 + 1
    return npx, npy, npz


def primalFlattener(diag_filename):
    if "_B" in diag_filename:
        print("Converting B fields")
        return BtoFlatPrimal
    elif "_E" in diag_filename:
        print("Converting E fields")
        return EtoFlatPrimal
    elif "_bulkVelocity" in diag_filename:
        print("Converting bulk velocity")
        return primalVectorToFlatPrimal
    elif "_density" in diag_filename:
        print("Converting ion density")
        return primalScalarToFlatPrimal
    else:
        raise ValueError(f"Unknown diagnostic {diag_filename}")


def max_nbr_levels_in(phare_h5):
    max_nbr_level = 0
    times_str = list(phare_h5["t"].keys())
    for time in times_str:
        nbrLevels = len(phare_h5["t"][time].keys())
        if max_nbr_level < nbrLevels:
            max_nbr_level = nbrLevels
    return max_nbr_level


def times_in(phare_h5):
    times_str = list(phare_h5["t"].keys())
    times = np.asarray([float(time) for time in times_str])
    times.sort()
    return times


def is_vector_data(patch):
    return len(patch.keys()) == 3


def level_spacing_from(root_spacing, ilvl):
    # hard-coded 2D adds 0 for last dim spacing
    return [dl / 2**ilvl for dl in root_spacing]


def make3d(root_spacing):
    if len(root_spacing) == 1:
        return root_spacing + [0, 0]
    elif len(root_spacing) == 2:
        return root_spacing + [0]
    return root_spacing


def main():

    if len(sys.argv) != 2 or sys.argv[1] in ["-h", "--help"]:
        print(f"Usage: {os.path.basename(sys.argv[0])} <path_to_phare_h5>")
        print("Works for EM fields, bulk velocity and density")
        sys.exit(1)

    path = sys.argv[1]
    phare_h5 = h5py.File(path, "r")
    root_spacing = make3d(phare_h5.attrs["cell_width"])
    times = times_in(phare_h5)
    max_nbr_level = max_nbr_levels_in(phare_h5)
    numberOfTimes = times.size

    phare_fn = os.path.basename(path).split(".")[0]
    data_directory = os.path.dirname(path)
    vtk_fn = f"{data_directory}/{phare_fn}.vtkhdf"

    toFlatPrimal = primalFlattener(phare_fn)

    print("PHARE H5 to VTKHDF conversion")
    print("-----------------------------")
    print(f"Converting {phare_fn} to {vtk_fn} in {data_directory}")
    print(f"Max number of levels : {max_nbr_level}")
    print(f"Number of time steps : {numberOfTimes}")

    vtk = h5py.File(vtk_fn, "w")
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

    steps = root.create_group("Steps")
    steps.attrs.create("NSteps", data=numberOfTimes, dtype="i8")
    steps.create_dataset("Values", data=times[:numberOfTimes])

    for ilvl in range(max_nbr_level)[:]:

        print(f"Processing level {ilvl}")

        lvl = root.create_group(f"Level{ilvl}")
        lvl_spacing = level_spacing_from(root_spacing, ilvl)
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
            AMRBoxOffsets += [amr_box_offset]
            # pass this time if this level does not exist
            if phare_lvl_name not in phare_h5["t"][time_str].keys():
                print(f"no level {ilvl} at time {time}")
                continue
            phare_lvl = phare_h5["t"][time_str][phare_lvl_name]

            for patch_id in list(phare_lvl.keys())[:]:
                patch = phare_lvl[patch_id]
                box = boxFromPatch(patch)
                AMRBox.append(box)
                nbr_boxes += 1
                npx, npy, npz = nbrNodes(box)

                if is_vector_data(patch):
                    x_name, y_name, z_name = list(patch.keys())
                    ph_x = patch[x_name][:]
                    ph_y = patch[y_name][:]
                    ph_z = patch[z_name][:]

                    data = toFlatPrimal(ph_x, ph_y, ph_z, npx, npy, npz)

                    if first:
                        # this is the first patch of the first time
                        # for this level, we need to create the dataset
                        # the first dimension is the total # of points
                        # which is unknown hence the None for maxshape
                        pointData = pointData_g.create_dataset(
                            "data", data=data, maxshape=(None, 3)
                        )
                        first = False

                    else:
                        # dataset already created with shape (current_size,3)
                        # we add b.shape[0] points (=npx*npy) to the first dim
                        # hence need to resize the dataset.
                        pointData.resize(current_size + data.shape[0], axis=0)
                        pointData[current_size:, :] = data
                    # pass

                    current_size += data.shape[0]
                else:
                    assert len(patch.keys()) == 1
                    dataset_name = list(patch.keys())[0]
                    ph_data = patch[dataset_name][:]

                    data = toFlatPrimal(ph_data, npx, npy, npz)

                    if first:
                        # 1D resizable datasets maxshape MUST be (None,) tuple
                        pointData = pointData_g.create_dataset(
                            "data",
                            data=data,
                            maxshape=(None,),
                        )
                        first = False

                    else:
                        pointData.resize(current_size + data.shape[0], axis=0)
                        pointData[current_size:] = data
                    # pass

                    current_size += data.shape[0]

            # of of the patch loops at that time

            # number of boxes added for ilvl across all times
            amr_box_offset += nbr_boxes
            step_nbrBoxes += [nbr_boxes]

        lvl.create_dataset("AMRBox", data=AMRBox)
        steps_lvl.create_dataset("AMRBoxOffset", data=AMRBoxOffsets)
        steps_lvl.create_dataset("NumberOfAMRBox", data=step_nbrBoxes)
        pointDataOffset_g.create_dataset("data", data=dataOffsets)

    vtk.close()


if __name__ == "__main__":
    main()
