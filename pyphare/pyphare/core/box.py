import numpy as np
from .phare_utilities import np_array_ify, is_scalar, is_nd_array


class Box:
    """represents a box in AMR index cell
    lower, upper : lower and upper cell AMR indexes
    """

    def __init__(self, lower, upper):
        lower, upper = [np_array_ify(arr) for arr in [lower, upper]]
        assert lower.shape == upper.shape
        assert (lower <= upper).all()
        self.lower = lower.astype(int)  # can't slice with floats
        self.upper = upper.astype(int)
        self.ndim = len(self.lower)

    def __mul__(self, box2):
        """
        return the box resulting from the intersection
        of box2 with self. No intersection returns None
        """
        box1 = self

        lower, upper = np.maximum(box1.lower, box2.lower), np.minimum(
            box1.upper, box2.upper
        )
        if (lower <= upper).all():
            return Box(lower, upper)

    @property
    def shape(self):
        """returns the length per dimension"""
        return (self.upper - self.lower) + 1

    def nCells(self):
        """returns the number of cells in the box"""
        return self.shape.prod()

    def __str__(self):
        return "Box({},{})".format(self.lower.tolist(), self.upper.tolist())

    def __repr__(self):
        return self.__str__()

    def __contains__(self, item):
        """true if item is completely within self"""
        if not isinstance(item, Box):
            item = np_array_ify(item)

        if is_nd_array(item):
            assert len(item) == self.ndim
            item = Box(item, item)

        return (item.lower >= self.lower).all() and (item.upper <= self.upper).all()

    def __eq__(self, other):
        return (
            isinstance(other, Box)
            and (self.lower == other.lower).all()
            and (self.upper == other.upper).all()
        )

    def __sub__(self, other):
        assert isinstance(other, Box)
        return remove(self, other)

    def copy(self):
        return Box(self.lower.copy(), self.upper.copy())


class nDBox(Box):
    def __init__(self, dim, lower, upper):
        def _get(self, p):
            return np.asarray([p] * dim)

        super().__init__(_get(dim, lower), _get(dim, upper))


class Box1D(nDBox):
    def __init__(self, lower, upper):
        super().__init__(1, lower, upper)


class Box2D(nDBox):
    def __init__(self, lower, upper):
        super().__init__(2, lower, upper)


class Box3D(nDBox):
    def __init__(self, lower, upper):
        super().__init__(3, lower, upper)


def refine(box, ratio):
    # box [2,10] refined with ratio 2 gives [4,21]
    # upper = 21 because cell 10 is split into 20,21
    # box [2,10] becomes [6,32] with ratio 3. cell 10 becomes cells 30,31,32
    # thus upper is former upper*ratio + ratio - 1
    return Box(box.lower * ratio, box.upper * ratio + ratio - 1)


def coarsen(box, ratio):
    return Box(box.lower / ratio, ((box.upper + 1) / ratio) - 1)


def shift(box, offset):
    return Box(box.lower + offset, box.upper + offset)


def grow(box, size):
    if is_scalar(size) and box.ndim > 1:
        raise ValueError("box.py: grow must use a list for dimension > 1")
    if (np.asarray(size) < 0).any():
        raise ValueError("size must be >=0")
    return Box(box.lower - size, box.upper + size)


def shrink(box, size):
    if is_scalar(size) and box.ndim > 1:
        raise ValueError("box.py: shrink must use a list for dimension > 1")
    if (np.asarray(size) < 0).any():
        raise ValueError("size must be >=0")
    return Box(box.lower + size, box.upper - size)


def remove(box, to_remove):
    """
    removes "remove" from the box
    this operation returns the list of the boxes
    that are not part of the intersection box*remove
    """
    intersection = box * to_remove

    if intersection is None:
        return [box]

    def copy(arr, replace):
        cpy = np.copy(arr)
        for i, v in replace.items():
            cpy[i] = v
        return cpy

    boxes = {}

    if intersection.lower[0] > box.lower[0]:
        boxes["left"] = Box(box.lower, copy(box.upper, {0: intersection.lower[0] - 1}))
    if intersection.upper[0] < box.upper[0]:
        boxes["right"] = Box(copy(box.lower, {0: intersection.upper[0] + 1}), box.upper)

    if box.ndim > 1:
        minx = intersection.lower[0] if "left" in boxes else box.lower[0]
        maxx = intersection.upper[0] if "right" in boxes else box.upper[0]
        if intersection.lower[1] > box.lower[1]:
            boxes["down"] = Box(
                copy(box.lower, {0: minx}),
                copy(box.upper, {0: maxx, 1: intersection.lower[1] - 1}),
            )
        if intersection.upper[1] < box.upper[1]:
            boxes["up"] = Box(
                copy(box.lower, {0: minx, 1: intersection.upper[1] + 1}),
                copy(box.upper, {0: maxx}),
            )

    if box.ndim > 2:
        miny = intersection.lower[1] if "down" in boxes else box.lower[1]
        maxy = intersection.upper[1] if "up" in boxes else box.upper[1]
        if intersection.lower[2] > box.lower[2]:
            boxes["back"] = Box(
                copy(box.lower, {0: minx, 1: miny}),
                copy(intersection.lower - 1, {0: maxx, 1: maxy}),
            )
        if intersection.upper[2] < box.upper[2]:
            boxes["front"] = Box(
                copy(intersection.upper + 1, {0: minx, 1: miny}),
                copy(box.upper, {0: maxx, 1: maxy}),
            )

    return list(boxes.values())


def amr_to_local(box, ref_box):
    return Box(box.lower - ref_box.lower, box.upper - ref_box.lower)


def select(data, box):
    return data[tuple([slice(l, u + 1) for l, u in zip(box.lower, box.upper)])]

class DataSelector:
    """
        can't assign return to function unless []
        usage
        DataSelector(data)[box] = val
    """
    def __init__(self, data):
        self.data = data
    def __getitem__(self, box_or_slice):
        if isinstance(box_or_slice, Box):
            return select(self.data, box_or_slice)
        return self.data[box_or_slice]

    def __setitem__(self, box_or_slice, val):
        self.__getitem__(box_or_slice)[:] = val

