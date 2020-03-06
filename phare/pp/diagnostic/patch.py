class PatchLevel:
    REFINEMENT_RATIO = 2
    dim_to_xyz = {0: "x", 1: "y", 2: "z"}
    xyz_to_dim = {"x": 0, "y": 1, "z": 2}

    def __init__(self, diag, idx, cells):
        self.diag, self.idx = diag, idx
        power = pow(2, idx)
        self.cells = [power * x for x in diag.sim.cells]
        self.cell_width_ = [x / power for x in diag.dl]
        self.patches, self.patchDict = [], {}

    def is_root(self):
        return self.idx == 0

    def position_to_index(self, pos, direction):
        return round(pos / self.cell_width(direction))

    def index_to_position(self, idx, direction):
        return self.cell_width(direction) * idx

    def cell_width(self, direction):
        return self.cell_width_[PatchLevel.xyz_to_dim[direction]]


class Patch:
    def __init__(self, patch_level, h5patch, dtype):
        self.patch_level = patch_level
        self.h5patch = h5patch
        self.id = h5patch.name.split("/")[-1][1:]  # samrai patch id e.g. 0x0
        self.dtype = dtype
        self.origin = [float(v) for v in h5patch.attrs["origin"].split(",")]
        self.cells = [int(v) for v in h5patch.attrs["nbrCells"].split(",")]
        self._reset_points()

    def _reset_points(self):
        self.points = [
            (
                self.min_coord(PatchLevel.dim_to_xyz[i]),
                self.max_coord(PatchLevel.dim_to_xyz[i]),
            )
            for i in range(self.patch_level.diag.dim)
        ]
        return self

    # we copy to transform patch origins for periodic overlap calculations
    def copy(self, transform=[0, 0, 0]):
        p = Patch(self.patch_level, self.h5patch, self.dtype)
        p.origin = [f + transform[i] for i, f in enumerate(self.origin)]
        return p._reset_points()

    def min_coord(self, direction):
        return self.origin[PatchLevel.xyz_to_dim[direction]]

    def max_coord(self, direction):
        return self.min_coord(direction) + self.patch_level.index_to_position(
            self.cells[PatchLevel.xyz_to_dim[direction]], direction
        )
