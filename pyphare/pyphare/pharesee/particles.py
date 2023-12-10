import numpy as np
from ..core.phare_utilities import refinement_ratio, print_trace


class Particles:
    """
    this class represent a set of particles
    particles can either be loaded randomly in a given box or have there attribute set from caller
    """

    def __init__(self, **kwargs):
        if "box" in kwargs:
            box = kwargs["box"]

            self.iCells = np.random.randint(
                box.lower, high=box.upper + 1, size=(box.nCells() * 100, box.ndim)
            )
            self.deltas = np.random.rand(box.nCells() * 100, box.ndim)
            self.v = np.random.randn(box.nCells() * 100, 3)
            self.weights = np.zeros(self.deltas.shape[0]) + 0.01
            self.charges = np.zeros_like(self.weights) + 1
            self.dl = np.zeros((self.weights.size, box.ndim)) + 0.1
            ndim = len(box.lower)

        else:
            self.iCells = kwargs["icells"][:]
            self.deltas = kwargs["deltas"][:]
            self.v = kwargs["v"][:]
            self.weights = kwargs["weights"][:]
            self.charges = kwargs["charges"][:]
            self.dl = kwargs["dl"][:]
            ndim = self.iCells.shape[1]

        self._x = None
        self._y = None

        assert self.iCells.ndim == self.deltas.ndim
        assert self.iCells.shape == self.deltas.shape
        assert self.iCells.shape[0] == self.size()
        self.ndim = ndim

    def xyz(self, i=0):
        return self.dl[:, i] * (self.iCells[:, i] + self.deltas[:, i])

    @property
    def x(self):
        if self._x is None:
            self._x = self.xyz()
        return self._x

    @property
    def y(self):
        if self._y is None:
            self._y = self.xyz(1)
        return self._y

    def _reset(self):
        self._x = None
        self._y = None

    def add(self, particles):
        self.iCells = np.concatenate((self.iCells, particles.iCells))
        self.deltas = np.concatenate((self.deltas, particles.deltas))
        self.v = np.concatenate((self.v, particles.v))
        self.charges = np.concatenate((self.charges, particles.charges))
        self.weights = np.concatenate((self.weights, particles.weights))
        self.dl = np.concatenate((self.dl, particles.dl))
        self._reset()

    def shift_icell(self, offset):
        self.iCells += offset
        self._reset()
        return self

    def size(self):
        return len(self.weights)

    def __eq__(self, that):
        if isinstance(that, Particles):
            # fails on OSX for some reason
            set_check = set(self.as_tuples()) == set(that.as_tuples())
            if set_check:
                return True
            try:
                all_assert_sorted(self, that)
                return True
            except AssertionError as ex:
                print(f"particles.py:Particles::eq failed with:", ex)
                print_trace()
                return False

        return False

    def select(self, box, box_type="cell"):
        """
        select particles from the given box
        assumption, box has AMR indexes of the same level as the data that the current instance is created from
        """
        assert len(box.lower) == self.ndim

        if box_type == "cell":
            if self.ndim == 1:
                idx = np.where((self.iCells >= box.lower) & (self.iCells <= box.upper))[
                    0
                ]
            else:

                def isin(p, b):
                    return p in b

                idx = np.where(np.apply_along_axis(isin, 1, self.iCells, box))[0]

        elif box_type == "pos":
            assert self.ndim == 1  # unhandled otherwise
            idx = np.where((self.x > box.lower[0]) & (self.x < box.upper[0]))[0]

        else:
            raise ValueError("unsupported box type ({})".format(box_type))

        return Particles(
            icells=self.iCells[idx],
            deltas=self.deltas[idx],
            v=self.v[idx, :],
            weights=self.weights[idx],
            charges=self.charges[idx],
            dl=self.dl[idx],
        )

    def __getitem__(self, box):
        return self.select(box)

    def erase(self, idx):
        self.iCells = np.delete(self.iCells, idx, axis=0)
        self.deltas = np.delete(self.deltas, idx, axis=0)
        self.v = np.delete(self.v, idx, axis=0)
        self.weights = np.delete(self.weights, idx, axis=0)
        self.charges = np.delete(self.charges, idx, axis=0)
        self.dl = np.delete(self.dl, idx, axis=0)

    def pop(self, idx):
        particles = Particles(
            icells=self.iCells[idx].copy(),
            deltas=self.deltas[idx].copy(),
            v=self.v[idx].copy(),
            weights=self.weights[idx].copy(),
            charges=self.charges[idx].copy(),
            dl=self.dl[idx].copy(),
        )
        self.erase(idx)
        return particles

    def split(self, sim):  # REQUIRES C++ PYBIND PHARE LIB
        from pyphare.cpp import split_pyarrays_fn

        split_pyarrays = split_pyarrays_fn(
            sim.ndim, sim.interp_order, sim.refined_particle_nbr
        )((self.iCells, self.deltas, self.weights, self.charges, self.v))
        return Particles(
            icells=split_pyarrays[0].reshape(
                int(len(split_pyarrays[0]) / self.ndim), self.ndim
            ),
            deltas=split_pyarrays[1].reshape(
                int(len(split_pyarrays[1]) / self.ndim), self.ndim
            ),
            weights=split_pyarrays[2].reshape(int(len(split_pyarrays[2])), 1),
            charges=split_pyarrays[3].reshape(int(len(split_pyarrays[3])), 1),
            v=split_pyarrays[4].reshape(int(len(split_pyarrays[4]) / 3), 3),
            dl=self.dl[0] / refinement_ratio
            + np.zeros((split_pyarrays[2].size, self.ndim)),
        )

    def as_tuples(self):
        return [
            (
                *self.iCells[i],
                *self.deltas[i],
                *self.v[i],
                *self.dl[i],
                *self.weights[i],
                *self.charges[i],
            )
            for i in range(self.size())
        ]


def all_assert_sorted(part1, part2):
    idx1 = _arg_sort(part1)
    idx2 = _arg_sort(part2)

    np.testing.assert_equal(part1.ndim, part2.ndim)
    np.testing.assert_equal(part1.size(), part2.size())

    deltol = (
        1e-6
        if any([part.deltas.dtype == np.float32 for part in [part1, part2]])
        else 1e-12
    )

    np.testing.assert_array_equal(part1.iCells[idx1], part2.iCells[idx2])
    np.testing.assert_allclose(part1.deltas[idx1], part2.deltas[idx2], atol=deltol)

    np.testing.assert_allclose(part1.v[idx1, 0], part2.v[idx2, 0], atol=1e-12)
    np.testing.assert_allclose(part1.v[idx1, 1], part2.v[idx2, 1], atol=1e-12)
    np.testing.assert_allclose(part1.v[idx1, 2], part2.v[idx2, 2], atol=1e-12)


def any_assert(part1, part2):
    np.testing.assert_equal(part1.size(), part2.size())
    assert part1.size() < 10  # slowness

    for part1_idx in range(part1.size()):
        assert part1.deltas[part1_idx] in part2.deltas
        assert part1.iCells[part1_idx] in part2.iCells
        assert part1.v[part1_idx] in part2.v
        assert part1.weights[part1_idx] in part2.weights
        assert part1.charges[part1_idx] in part2.charges
        assert part1.dl[part1_idx] in part2.dl


def aggregate(particles_in):
    assert all([isinstance(particles, Particles) for particles in particles_in])

    from copy import copy

    particles_out = copy(particles_in[0])  # use first, concat rest
    for particles in particles_in[1:]:
        particles_out.add(particles)

    assert particles_out.size() == sum([particles.size() for particles in particles_in])
    return particles_out


def remove(particles, idx):
    """
    returns a Particles object where particles indexed "idx"
    have been removed from "particles"
    """
    if len(idx) == particles.size():
        return None

    icells = np.delete(particles.iCells, idx, axis=0)
    deltas = np.delete(particles.deltas, idx, axis=0)
    vx = np.delete(particles.v[:, 0], idx)
    vy = np.delete(particles.v[:, 1], idx)
    vz = np.delete(particles.v[:, 2], idx)
    v = np.zeros((vx.size, 3))
    v[:, 0] = vx
    v[:, 1] = vy
    v[:, 2] = vz
    weights = np.delete(particles.weights, idx, axis=0)
    dl = np.zeros((len(weights), particles.ndim))
    for i in range(particles.ndim):
        dl[:, i] = np.delete(particles.dl[:, i], idx)

    charges = np.delete(particles.charges, idx, axis=0)

    return Particles(
        icells=icells, deltas=deltas, v=v, weights=weights, charges=charges, dl=dl
    )


def _arg_sort(particles):
    x1 = particles.iCells[:, 0] + particles.deltas[:, 0]

    if particles.ndim == 1:
        return np.argsort(x1)

    if particles.ndim > 1:
        y1 = particles.iCells[:, 1] + particles.deltas[:, 1]
        return np.argsort(np.sqrt((x1**2 + y1**2)) / (x1 / y1))

    if particles.ndim == 3:
        z1 = particles.iCells[:, 2] + particles.deltas[:, 2]
        return np.argsort(np.sqrt((x1**2 + y1**2 + z1**2)) / (x1 / y1 / z1))
