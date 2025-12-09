import numpy as np

from ...core import phare_utilities as phut
from ...core import box as boxm
from ...core.box import Box


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
        self.box = layout.box
        self.origin = layout.origin
        self.layout = layout

    def __deepcopy__(self, memo):
        no_copy_keys = ["dataset"]  # do not copy these things
        return phut.deep_copy(self, memo, no_copy_keys)


class FieldData(PatchData):
    """
    Concrete type of PatchData representing a physical quantity
    defined on a grid.
    """

    @property
    def x(self):
        withGhosts = self.field_name != "tags"
        if self._x is None:
            self._x = self.layout.yeeCoordsFor(
                self.field_name,
                "x",
                withGhosts=withGhosts,
                centering=self.centerings[0],
            )
        return self._x

    @property
    def y(self):
        withGhosts = self.field_name != "tags"
        if self._y is None:
            self._y = self.layout.yeeCoordsFor(
                self.field_name,
                "y",
                withGhosts=withGhosts,
                centering=self.centerings[1],
            )
        return self._y

    @property
    def z(self):
        withGhosts = self.field_name != "tags"
        if self._z is None:
            self._z = self.layout.yeeCoordsFor(
                self.field_name,
                "z",
                withGhosts=withGhosts,
                centering=self.centerings[2],
            )
        return self._z

    def primal_directions(self):
        return self.size - self.ghost_box.shape

    def __str__(self):
        return "FieldData: (box=({}, {}), key={})".format(
            self.layout.box, self.layout.box.shape, self.field_name
        )

    def __repr__(self):
        return self.__str__()

    def compare(self, that, atol=1e-16):
        return self.field_name == that.field_name and phut.fp_any_all_close(
            self.dataset[:], that.dataset[:], atol=atol
        )

    def __eq__(self, that):
        return self.compare(that)

    def __ne__(self, that):
        return not (self == that)

    def select(self, box):
        """
        return view of internal data based on overlap of input box
           returns a view +1 in size in primal directions
        """
        assert isinstance(box, Box) and box.ndim == self.box.ndim

        gbox = self.ghost_box.copy()
        gbox.upper += self.primal_directions()

        box = box.copy()
        box.upper += self.primal_directions()

        overlap = box * gbox
        if overlap is not None:
            lower = self.layout.AMRToLocal(overlap.lower)
            upper = self.layout.AMRToLocal(overlap.upper)

            if box.ndim == 1:
                return self.dataset[lower[0] : upper[0] + 1]
            if box.ndim == 2:
                return self.dataset[lower[0] : upper[0] + 1, lower[1] : upper[1] + 1]
        return np.array([])

    def __getitem__(self, box_or_slice):
        if isinstance(box_or_slice, slice):
            return self.dataset[box_or_slice]
        return self.select(box_or_slice)

    def __setitem__(self, box_or_slice, val):
        self.__getitem__(box_or_slice)[:] = val

    def __init__(self, layout, field_name, data, **kwargs):
        """
        :param layout: A GridLayout representing the domain on which data is defined
        :param field_name: the name of the field (e.g. "Bx")
        :param data: the dataset from which data can be accessed
        """
        super().__init__(layout, "field")
        self._x = None
        self._y = None
        self._z = None

        self.field_name = field_name
        self.name = field_name
        self.dl = np.asarray(layout.dl)
        self.ndim = layout.box.ndim
        self.ghosts_nbr = np.zeros(self.ndim, dtype=int)

        if field_name in layout.centering["X"]:
            directions = ["X", "Y", "Z"][: layout.box.ndim]  # drop unused directions
            self.centerings = [
                layout.qtyCentering(field_name, direction) for direction in directions
            ]
        elif "centering" in kwargs:
            if isinstance(kwargs["centering"], list):
                self.centerings = kwargs["centering"]
                assert len(self.centerings) == self.ndim
            else:
                if self.ndim != 1:
                    raise ValueError(
                        "FieldData invalid dimenion for centering argument, expected list for dim > 1"
                    )
                self.centerings = [kwargs["centering"]]
        else:
            raise ValueError(
                f"centering not specified and cannot be inferred from field name : {field_name}"
            )

        if self.field_name != "tags":
            for i, centering in enumerate(self.centerings):
                self.ghosts_nbr[i] = layout.nbrGhosts(layout.interp_order, centering)

        self.ghost_box = boxm.grow(layout.box, self.ghosts_nbr)

        self.size = np.copy(self.ghost_box.shape)
        self.offset = np.zeros(self.ndim)

        for i, centering in enumerate(self.centerings):
            if centering == "primal":
                self.size[i] = self.ghost_box.shape[i] + 1
            else:
                self.size[i] = self.ghost_box.shape[i]
                self.offset[i] = 0.5 * self.dl[i]

        self.dataset = data

    def meshgrid(self, select=None):
        def grid():
            if self.ndim == 1:
                return [self.x]
            if self.ndim == 2:
                return np.meshgrid(self.x, self.y, indexing="ij")
            return np.meshgrid(self.x, self.y, self.z, indexing="ij")

        mesh = grid()
        if select is not None:
            return tuple(g[select] for g in mesh)
        return mesh

    def zeros_like(self):
        from copy import deepcopy

        copy = deepcopy(self)
        assert id(copy.dataset) == id(self.dataset)
        copy.dataset = np.zeros(deepcopy(self.dataset[:].shape))
        assert id(copy.dataset) != id(self.dataset)
        return copy


class ParticleData(PatchData):
    """
    Concrete type of PatchData representing particles in a region
    """

    def __init__(self, layout, data, pop_name):
        """
        :param layout: A GridLayout object representing the domain in which particles are
        :param data: dataset containing particles
        """
        super().__init__(layout, "particles")
        self.dataset = data
        self.pop_name = pop_name
        self.name = pop_name
        self.ndim = layout.box.ndim

        self.pop_name = pop_name
        if layout.interp_order == 1:
            self.ghosts_nbr = np.array([1] * layout.box.ndim)
        elif layout.interp_order == 2 or layout.interp_order == 3:
            self.ghosts_nbr = np.array([2] * layout.box.ndim)
        else:
            raise RuntimeError(
                "invalid interpolation order {}".format(layout.interp_order)
            )

        self.ghost_box = boxm.grow(layout.box, self.ghosts_nbr)
        assert (self.box.lower == self.ghost_box.lower + self.ghosts_nbr).all()

    def select(self, box):
        return self.dataset[box]

    def __getitem__(self, box):
        return self.select(box)

    def size(self):
        return self.dataset.size()

    def compare(self, that, *args, **kwargs):
        """args/kwargs may include atol for consistency with field::compare"""
        return self.name == that.name and self.dataset == that.dataset

    def __eq__(self, that):
        return self.compare(that)
