


class Box:
    """represents a box in AMR index cell
    lower, upper : lower and upper cell AMR indexes
    """

    def __init__(self, lower, upper):
        assert lower <= upper
        self.lower = lower
        self.upper = upper

    def __mul__(self, box2):
        """
        return the box resulting from the intersection
        of box2 with self. No intersection returns None
        """
        box1 = self

        lower, upper = max(box1.lower, box2.lower), min(box1.upper, box2.upper)
        if lower <= upper:
            return Box(lower, upper)





    def size(self):
        """returns the number of cells in the box"""
        # later should return that number in each direction
        return self.upper - self.lower + 1



    def __str__(self):
        return "[ {lower},{upper} ]".format(lower=self.lower, upper=self.upper)


    def __repr__(self):
        return self.__str__()


    def __contains__(self, item):
        if isinstance(item, int):
            cell = item
            return (cell >= self.lower) and (cell <= self.upper)

        if isinstance(item, Box):
            box = item

            return box.lower >= self.lower and box.upper <= self.upper

    def __eq__(self, other):
        return self.lower == other.lower and self.upper == other.upper


def refine(box, ratio):
    # box [2,10] refined with ratio 2 gives [4,21]
    # upper = 21 because cell 10 is split into 20,21
    # box [2,10] becomes [6,32] with ratio 3. cell 10 becomes cells 30,31,32
    # thus upper is former upper*ratio + ratio - 1
    return Box(box.lower * ratio, box.upper * ratio + ratio  -1)


def shift(box, offset):
    return Box(box.lower + offset, box.upper + offset)


def grow(box, size):
    # in multiple dim, size could be a tuple
    # with a number of cell to grow the box in each dir.
    if size < 0:
        raise ValueError("size must be >=0")
    new_box = Box(box.lower, box.upper)
    new_box.lower = box.lower - size
    new_box.upper = box.upper + size
    return new_box


def remove(box, to_remove):
    """
    removes "remove" from the box
    this operation returns the list of the boxes
    that are not part of the intersection box*remove
    """
    intersection = box * to_remove

    if intersection is None:
        return [box, ]


    if to_remove in box and box.lower < to_remove.lower and box.upper > to_remove.upper:
        #   |----------------|    box
        #        |-----|          remove
        return [Box(box.lower, intersection.lower - 1), Box(intersection.upper + 1, box.upper)]

    if box in to_remove:
        #    |--------------|     box
        #  |-------------------|  remove
        return []

    if box.lower in intersection and box.upper > to_remove.upper:
        #       |---------------| box
        #  |-----------|          remove
        return [Box(intersection.upper + 1, box.upper), ]

    if box.upper in to_remove and box.lower < to_remove.lower:
        #  |---------------|      box
        #            |----------| remove
        return [Box(box.lower, intersection.lower - 1), ]




def amr_to_local(box, ref_box):
    return Box(box.lower - ref_box.lower, box.upper - ref_box.lower)
