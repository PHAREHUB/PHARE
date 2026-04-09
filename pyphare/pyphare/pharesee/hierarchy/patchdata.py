#
#
#


import numpy as np


from ...core import gridlayout
from ...core import box as boxm
from ...core import phare_utilities as phut


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
        cpy = phut.deep_copy(self, memo, no_copy_keys)
        cpy.dataset = self.dataset
        return cpy


class FieldData(PatchData):
    """
    Concrete type of PatchData representing a physical quantity
    defined on a grid.
    """

    @property
    def x(self):
        if self._x is None:
            self._x = self.yeeCoordsFor(0)
        return self._x

    @property
    def y(self):
        if self._y is None:
            self._y = self.yeeCoordsFor(1)
        return self._y

    @property
    def z(self):
        if self._z is None:
            self._z = self.yeeCoordsFor(2)
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
        try:
            phut.assert_fp_any_all_close(self[:], that[:], atol=atol)
        except AssertionError as e:
            return phut.EqualityCheck(False, str(e))
        return self.field_name == that.field_name

    def __eq__(self, that):
        return self.compare(that)

    def __ne__(self, that):
        return not (self == that)

    def select(self, box_or_slice):
        """
        return view of internal data based on overlap of input box
           returns a view +1 in size in primal directions
        """
        if isinstance(box_or_slice, (slice, list, tuple)):
            return self.dataset[box_or_slice]

        box = box_or_slice
        assert isinstance(box, boxm.Box) and box.ndim == self.box.ndim

        gbox = self.ghost_box.copy()
        gbox.upper += self.primal_directions()

        box = box.copy()
        box.upper += self.primal_directions()

        overlap = box * gbox
        if overlap is not None:
            lower = self.layout.AMRToLocal(overlap.lower)
            upper = self.layout.AMRToLocal(overlap.upper)
            select = tuple(slice(lower[i], upper[i] + 1) for i in range(box.ndim))
            return self.dataset[select]

        return np.array([])

    def __getitem__(self, box_or_slice):
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

        self.dl = np.asarray(layout.dl)
        self.field_name = field_name
        self.name = field_name
        self.ndim = layout.box.ndim

        self.centerings = self._resolve_centering(**kwargs)
        self.ghosts_nbr = self._resolve_ghost_nbr(**kwargs)

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

    def copy_as(self, data=None, **kwargs):
        data = data if data is not None else self.dataset
        name = kwargs.get("name", self.field_name)
        kwargs["centering"] = kwargs.get("centering", self.centerings)
        kwargs["ghosts_nbr"] = kwargs.get("ghosts_nbr", self.ghosts_nbr)
        kwargs["layout"] = kwargs.get("layout", self.layout)
        return FieldData(field_name=name, data=data, **kwargs)

    def yeeCoordsFor(self, idx):
        return self.layout.yeeCoordsFor(
            self.field_name,
            gridlayout.directions[idx],
            withGhosts=any(self.ghosts_nbr) and self.field_name != "tags",
            centering=self.centerings[idx],
            nbrGhosts=self.ghosts_nbr,
        )

    def _resolve_ghost_nbr(self, **kwargs):
        layout = self.layout
        ghosts_nbr = kwargs.get("ghosts_nbr", np.zeros(self.ndim, dtype=int))
        if "ghosts_nbr" not in kwargs:
            if self.field_name != "tags":
                for i, centering in enumerate(self.centerings):
                    ghosts_nbr[i] = layout.nbrGhosts(layout.interp_order, centering)
        return phut.np_array_ify(ghosts_nbr, layout.box.ndim)

    def _resolve_centering(self, **kwargs):
        field_name = self.field_name
        if "centering" in kwargs:
            if isinstance(kwargs["centering"], list):
                assert len(kwargs["centering"]) == self.ndim
                return kwargs["centering"]
            else:
                if self.ndim != 1:
                    raise ValueError(
                        "FieldData invalid dimenion for centering argument, expected list for dim > 1"
                    )
                return phut.listify(kwargs["centering"])

        if field_name in self.layout.centering["X"]:
            directions = ["X", "Y", "Z"][: self.ndim]  # drop unused directions
            return [
                self.layout.qtyCentering(field_name, direction)
                for direction in directions
            ]
        raise ValueError(
            f"centering not specified and cannot be inferred from field name : {field_name}"
        )

    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        return field_data_array_ufunc(self, ufunc, method, *inputs, **kwargs)

    def __array_function__(self, func, types, args, kwargs):
        return field_data_array_function(self, func, types, args, kwargs)


def field_data_array_ufunc(patch_data, ufunc, method, *inputs, **kwargs):
    if method != "__call__":
        raise NotImplementedError

    in_ = [i.dataset if isinstance(i, FieldData) else i for i in inputs]
    out_ = getattr(ufunc, method)(*in_, **kwargs)

    if isinstance(out_, np.ndarray):
        return patch_data.copy_as(
            out_,
            layout=patch_data.layout,
            name=patch_data.field_name,
            centering=patch_data.centerings,
        )

    raise NotImplementedError


def field_data_array_function(patch_data, func, types, args, kwargs):
    in_ = [a.dataset if isinstance(a, FieldData) else a for a in args]
    out_ = func(*in_, **kwargs)

    if phut.is_scalar(out_):
        return out_

    return patch_data.copy_as(
        out_,
        name=patch_data.field_name,
        layout=patch_data.layout,
        centering=patch_data.centerings,
    )


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
