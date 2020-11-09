import numpy as np
from .phare_utilities import np_array_ify, is_scalar, is_nd_array





class Box:
    """represents a box in AMR index cell
    lower, upper : lower and upper cell AMR indexes
    """

    def __init__(self, lower, upper):
        lower, upper = [np_array_ify(arr) for arr in [lower, upper]]
        assert lower.size == upper.size
        assert (lower <= upper).all()
        self.lower = lower.astype(int) # can't slice with floats
        self.upper = upper.astype(int)

    def dim(self):
        return len(self.lower.shape)

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

    def shape(self):
        """returns the length per dimension"""
        return (self.upper - self.lower) + 1

    def size(self):
        """deprecated, use shape"""
        # DOTO remove use shape
        return self.shape()[0]

    def cells(self):
        """returns the number of cells in the box"""
        return self.shape().prod()


    def __str__(self):
        return "[ {lower},{upper} ]".format(lower=self.lower, upper=self.upper)

    def __repr__(self):
        return self.__str__()

    # only 0 or 1 are valid
    def __getitem__(self, idx):
        assert 0 <= idx <= 1
        if idx == 0:
            return self.lower
        return self.upper

    def __contains__(self, item):
        """true if item is completely within self"""
        if not isinstance(item, Box):
            item = np_array_ify(item)

        dims = len(self.lower)

        if is_nd_array(item):
            assert len(item) == dims
            item = Box(item, item)

        return (item.lower >= self.lower).all() and (item.upper <= self.upper).all()


    def __eq__(self, other):
        return (self.lower == other.lower).all() and (self.upper == other.upper).all()



class nDBox(Box):
    def __init__(self, dim, l, u):
        def _get(self, p):
            return [p for i in range(dim)]

        super().__init__(_get(dim, l), _get(dim, u))


class Box1D(nDBox):
    def __init__(self, l, u):
        super().__init__(1, l, u)


class Box2D(nDBox):
    def __init__(self, l, u):
        super().__init__(2, l, u)


class Box3D(nDBox):
    def __init__(self, l, u):
        super().__init__(3, l, u)


def refine(box, ratio):
    # box [2,10] refined with ratio 2 gives [4,21]
    # upper = 21 because cell 10 is split into 20,21
    # box [2,10] becomes [6,32] with ratio 3. cell 10 becomes cells 30,31,32
    # thus upper is former upper*ratio + ratio - 1
    return Box(box.lower * ratio, box.upper * ratio + ratio - 1)


def shift(box, offset):
    return Box(box.lower + offset, box.upper + offset)


def grow(box, size):
    if is_scalar(size):
        assert box.dim() == 1 # possible overkill, here to block accidents
    if (np.asarray(size) < 0).any():
        raise ValueError("size must be >=0")
    return Box(box.lower - size, box.upper + size)


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

    dims = len(box.lower)
    boxes = {}

    if intersection.lower[0] > box.lower[0]:
        boxes["left"] = Box(box.lower, copy(box.upper, {0: intersection.lower[0] - 1}))
    if intersection.upper[0] < box.upper[0]:
        boxes["right"] = Box(copy(box.lower, {0: intersection.upper[0] + 1}), box.upper)

    if dims > 1:
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

    if dims > 2:
        miny = intersection.lower[1] if "down" in boxes else box.lower[1]
        maxy = intersection.upper[1] if "up" in boxes else box.upper[1]
        if intersection.lower[2] > box.lower[2]:
            boxes["back"] = Box(copy(box.lower, {0: minx, 1: miny}),
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
