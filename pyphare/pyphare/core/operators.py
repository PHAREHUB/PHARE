import numpy as np
from pyphare.pharesee.hierarchy import PatchLevel, Patch, FieldData
from pyphare.pharesee.hierarchy import ScalarField, VectorField , Tensor2Field


def sqrt(hierarchy):

    """
    returns a Hierarchy of the same type of the input (which can be
    ScalarField, VectorField or Tensor2Field) containing the square root
    of each dataset

    """

    patch_levels = hierarchy.patch_levels
    domain_box = hierarchy.domain_box

    num_of_components = len(hierarchy.levels()[0].patches[0].patch_datas.keys())

    if num_of_components == 1:
        names = ['value']
    elif num_of_components == 3:
        names = ['x', 'y', 'z']
    elif num_of_components == 6:
        names = ['xx', 'xy', 'xz', 'yy', 'yz', 'zz']
    elif num_of_components == 9:
        names = ['xx', 'xy', 'xz', 'yx', 'yy', 'yz', 'zx', 'zy', 'zz']


    for time in list(hierarchy.time_hier.keys()):
        new_patch_level = {}

        for ilvl, lvl in hierarchy.levels(time).items():
            new_patches = {}

            for patch in lvl.patches:
                new_pd = {}
                layout = patch.layout

                num_of_components = len(patch.patch_datas.keys())

                for (name, pd_name) in zip(names, patch.patch_datas.keys()):
                    dset = np.sqrt(np.asarray(patch.patch_datas[pd_name].dataset))

                    pd = FieldData(layout, name, dset, centering=patch.patch_datas[pd_name].centerings)
                    new_pd[name] = pd

                if ilvl not in new_patches:
                    new_patches[ilvl] = []

                new_patches[ilvl].append(Patch(new_pd, patch.id))

            new_patch_level[ilvl] = PatchLevel(ilvl, new_patches[ilvl])

    if num_of_components == 1:
        return ScalarField(new_patch_level, domain_box, time=time)
    if num_of_components == 3:
        return VectorField(new_patch_level, domain_box, time=time)
    if num_of_components == 6 or num_of_components == 9:
        return Tensor2Field(new_patch_level, domain_box, time=time)

