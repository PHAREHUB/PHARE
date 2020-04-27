from ..core import box as boxm
import numpy as np

class PatchData:
    """
    base class for FieldData and ParticleData
    this class just factors common geometrical properties
    """
    def __init__(self, layout, quantity):
        """
        :param layout: a GridLayout representing the domain on which the data
             is defined
        :param quantity: ['field', 'particle']
        """
        self.quantity = quantity
        self.box     = layout.box
        self.origin   = layout.origin
        self.layout   = layout




class FieldData(PatchData):
    """
    Concrete type of PatchData representing a physical quantity
    defined on a grid.
    """
    def __init__(self, layout, field_name, data):
        """
        :param layout: A GridLayout representing the domain on which data is defined
        :param field_name: the name of the field (e.g. "Bx")
        :param data: the dataset from which data can be accessed
        """
        super().__init__(layout, 'field')

        self.layout = layout
        self.field_name = field_name
        self.dx = layout.dl[0]

        centering = layout.centering["X"][field_name]
        self._ghosts_nbr = layout.nbrGhosts(layout.interp_order, centering)
        self.ghost_box = boxm.grow(layout.box, self._ghosts_nbr)

        if centering == "primal":
            self.size = self.ghost_box.size() + 1
            offset = 0
        else:
            self.size = self.ghost_box.size()
            offset = 0.5*self.dx

        self.x = self.origin[0] - self._ghosts_nbr * self.dx + np.arange(self.size) * self.dx + offset
        self.dataset = data



class ParticleData(PatchData):
    """
    Concrete type of PatchData representing particles in a region
    """
    def __init__(self, layout, data):
        """
        :param layout: A GridLayout object representing the domain in which particles are
        :param data: dataset containing particles
        """
        super().__init__(layout, 'particles')
        self.dataset = data
        if layout.interp_order == 1:
            self._ghosts_nbr = 1
        elif layout.interp_order == 2 or layout.interp_order == 3:
            self._ghosts_nbr = 2
        else:
            raise RuntimeError("invalid interpolation order {}".format(layout.interp_order))
        self.ghost_box = boxm.grow(layout.box, self._ghosts_nbr)



class Patch:
    """
    A patch represents a hyper-rectangular region of space
    """

    def __init__(self, patch_datas):
        """
        :param patch_datas: a list of PatchData objects
        these are assumed to "belong" to the Patch so to
        share the same origin, mesh size and box.
        """
        pdata0 = list(patch_datas.values())[0] #0 represents all others
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

    def __init__(self, patch_levels, domain_box, refinement_ratio):
        self.patch_levels = patch_levels
        self.domain_box = domain_box
        self.refinement_ratio = refinement_ratio


    def refined_domain_box(self, level_number):
        """
        returns the domain box refined for a given level number
        """
        assert(level_number>=0)
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


    def plot(self):

        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(10, 3))
        for ilvl, lvl in enumerate(self.patch_levels):
            lvl_offset = ilvl * 0.1
            for patch in lvl.patches:
                dx = patch.dx
                origin = patch.origin
                x0 = patch.box.lower * dx
                x1 = patch.box.upper * dx
                xcells = np.arange(x0, x1 + dx, dx)
                y = lvl_offset + np.zeros_like(xcells)
                ax.plot(xcells, y, marker=".")

        fig.savefig("hierarchy.png")



def is_root_lvl(patch_level):
    return patch_level.level_number == 0







