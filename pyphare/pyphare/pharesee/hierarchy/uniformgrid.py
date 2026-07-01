#
#
#


from . import patchdata


class UniformGrid(patchdata.FieldData):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def plot(self, **kwargs):
        _plot_this(self, **kwargs)


class UniformGrids:
    """Basically TensorField version"""

    def __init__(self, grids: dict):
        self.grids = grids

    def plot(self, qty, **kwargs):
        self.grids[qty].plot(**kwargs)

    def items(self):
        return self.grids.items()

    def keys(self):
        return self.grids.keys()

    def values(self):
        return self.grids.values()

    def __len__(self):
        return len(self.grids)

    def __getitem__(self, key):
        return self.grids[key]

    def __setitem__(self, key, val):
        self.grids[key] = val

    def __contains__(self, key):
        return key in self.grids


def _plot_this(ungrid, **kwargs):
    from .plotting.plot_fields import plot_field_data

    plot_field_data(ungrid, **kwargs)


def _get_grid(obj, qty=None):
    if type(obj) is UniformGrid:
        return obj
    if type(obj) is UniformGrids:
        if qty is None:
            raise ValueError("UniformGrid::peak_coordinates only works on one quantity")
        return obj[qty]
    raise ValueError("UniformGrid::peak_coordinates only supports UniformGrids")


def peak_coordinates(obj, qty=None, **kwargs):
    grid = _get_grid(obj, qty)
    peaks, _ = find_peaks(grid, **kwargs)
    dl = grid.layout.dl
    return peaks * dl[0]  # only works for 1D!


def find_peaks(obj, qty=None, **kwargs):
    grid = _get_grid(obj, qty)

    if grid.ndim > 1:
        raise ValueError("UniformGrid::find_peaks only supports 1d simulations")

    from scipy.signal import find_peaks

    peaks, props = find_peaks(grid[:], **kwargs)
    if "peak_heights" in props:
        return peaks, props["peak_heights"]
    return peaks, props
