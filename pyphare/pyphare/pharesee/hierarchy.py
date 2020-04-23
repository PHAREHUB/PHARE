
from . import box as boxm
import numpy as np

class PatchData:

    def __init__(self, layout, quantity):
        self.quantity = quantity
        self._box     = layout.box
        self.origin   = layout.origin
        self.layout   = layout


class FieldData(PatchData):
    def __init__(self, layout, field_name, data):
        super().__init__(layout, 'field')

        self.layout = layout
        self.field_name = field_name
        self.dx = layout.dl[0]

        centering = layout.centering["X"][field_name]
        self._ghosts_nbr = layout.nbrGhosts(layout.interp_order, centering)
        self.ghost_box = boxm.grow(layout.box, self._ghosts_nbr)

        if centering == "primal":
            self.size = self.ghost_box.size() + 1
        else:
            self.size = self.ghost_box.size()

        self.x = self.origin[0] - self._ghosts_nbr * self.dx + np.arange(self.size) * self.dx
        self.dataset = data



class ParticleData(PatchData):
    def __init__(self, layout, data):
        super().__init__(layout, 'particles')
        self.dataset = data
        if layout.interp_order == 1:
            self._ghosts_nbr = 1
        elif layout.interp_order == 2 or layout.interp_order == 3:
            self._ghosts_nbr = 2
        else:
            raise RuntimeError("invalid interpolation order")
        self.ghost_box = boxm.grow(layout.box, self._ghosts_nbr)


class Patch:
    """
    A patch represents a hyper-rectangular region of space
    """

    def __init__(self, patch_datas):

        # all patch_datas are assumed to have the reference to
        # same layout object so take data from the first
        pdata0 = list(patch_datas.values())[0]
        self.box = pdata0.layout.box
        self.origin = pdata0.layout.origin
        self.dx = pdata0.layout.dl[0]
        self.patch_datas = patch_datas


class PatchLevel:
    """is a collection of patches """

    def __init__(self, lvl_nbr, patches):
        self.level_number = lvl_nbr
        self.patches = patches


class PatchHierarchy:
    """is a collection of patch levels """

    def __init__(self, levels, domain_box, ratio):
        self.patch_levels = levels
        self.domain_box = domain_box
        self.refinement_ratio = ratio

    def refined_domain_box(self, level_number):
        return boxm.refine(self.domain_box, self.refinement_ratio ** level_number)

    def __str__(self):
        s = "Hierarchy: \n"
        for ilvl, lvl in enumerate(self.patch_levels):
            s = s + "Level {}\n".format(ilvl)
            for ip, patch in enumerate(lvl.patches):
                for qty_name, pd in patch.patch_datas.items():
                    pdstr = "    P{ip} {pdname} box is {box} and ghost box is {gbox}"
                    s = s + pdstr.format(ip=ip, pdname=qty_name,
                                         box=patch.box, gbox=pd.ghost_box)
                    s = s + "\n"
        return s


def is_root_lvl(patch_level):
    return patch_level.level_number == 0







