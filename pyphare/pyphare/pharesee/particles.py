
import numpy as np
from ..core.phare_utilities import refinement_ratio, print_trace, is_scalar

class Particles:
    """
    this class represent a set of particles
    particles can either be loaded randomly in a given box or have there attribute set from caller
    """
    def __init__(self, **kwargs):
        if "box" in kwargs:
            box = kwargs["box"]

            self.iCells  = np.random.randint(box.lower, high=box.upper+1, size=(box.nCells()*100, box.ndim))
            self.deltas  = np.random.rand(box.nCells()*100, box.ndim)
            self.v       = np.random.randn(box.nCells()*100, 3)
            self.weights = np.zeros(self.deltas.shape[0]) + 0.01
            self.charges = np.zeros_like(self.weights) + 1
            self.dl      = np.zeros((self.weights.size, box.ndim))+0.1
            ndim = len(box.lower)

        else:
            self.iCells  = kwargs["icells"][:]
            self.deltas  = kwargs["deltas"][:]
            self.v       = kwargs["v"][:]
            self.weights = kwargs["weights"][:]
            self.charges = kwargs["charges"][:]
            self.dl      = kwargs["dl"][:]
            ndim = self.iCells.ndim

        self._x = None
        self._y = None
        assert len(self.weights.shape) == 1
        assert len(self.charges.shape) == 1

        assert self.iCells.ndim == self.deltas.ndim
        assert self.iCells.shape == self.deltas.shape

        # we like to maintain iCells/deltas in reshaped form per ndim
        needs_reshaping = ndim == 1 and int(len(self.iCells) / self.size()) > 1
        if needs_reshaping:
            nICells = len(self.iCells)
            assert nICells % self.size() == 0 # perfect division
            ndim = int(nICells / self.size())
            assert ndim in [1, 2, 3]
            self.iCells = np.asarray(self.iCells).reshape(int(nICells / ndim), ndim)
            self.deltas = np.asarray(self.deltas).reshape(int(nICells / ndim), ndim)

        assert self.iCells.shape[0] == self.size()
        assert self.iCells.ndim == self.deltas.ndim
        assert self.iCells.shape == self.deltas.shape

        self.ndim = ndim



    def xyz(self, i=0):
        if self.ndim == 1:
            return self.dl[:,0]*(self.iCells[:] + self.deltas[:])
        return self.dl[:,i]*(self.iCells[:,i] + self.deltas[:,i])


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
        self.iCells   = np.concatenate((self.iCells, particles.iCells))
        self.deltas   = np.concatenate((self.deltas, particles.deltas))
        self.v        = np.concatenate((self.v, particles.v))
        self.charges  = np.concatenate((self.charges, particles.charges))
        self.weights  = np.concatenate((self.weights, particles.weights))
        self.dl       = np.concatenate((self.dl, particles.dl))
        self._reset()


    def shift_icell(self, offset):
        self.iCells += offset
        self._reset()
        return self

    def size(self):
          return len(self.weights)


    def __eq__(self, that):
        if isinstance(that, Particles):
            try:
                all_assert(self, that)
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

        if box_type=="cell":

            if self.ndim == 1:
                idx = np.where((self.iCells >= box.lower) & (self.iCells <= box.upper))[0]
            else:
                def isin(p, b):
                    return p in b
                idx = np.where(np.apply_along_axis(isin, 1, self.iCells, box))[0]

        elif box_type=="pos":
            assert self.ndim == 1 # unhandled otherwise
            idx = np.where((self.x > box.lower[0]) & (self.x < box.upper[0]))[0]

        else:
            raise ValueError("unsupported box type ({})".format(box_type))

        return Particles(icells=self.iCells[idx],
                         deltas=self.deltas[idx],
                         v = self.v[idx,:],
                         weights=self.weights[idx],
                         charges=self.charges[idx],
                         dl = self.dl[idx])

    def __getitem__(self, box):
        return self.select(box)


    def erase(self, idx_or_list):
        idx_list = [idx_or_list] if is_scalar(idx_or_list) else idx_or_list
        self.iCells  = np.delete(self.iCells, idx_list, axis=0)
        self.deltas  = np.delete(self.deltas, idx_list, axis=0)
        self.v       = np.delete(self.v, idx_list, axis=0)
        self.weights = np.delete(self.weights, idx_list, axis=0)
        self.charges = np.delete(self.charges, idx_list, axis=0)
        self.dl      = np.delete(self.dl, idx_list, axis=0)


    def pop(self, idx_or_list):
        idx_list = [idx_or_list] if is_scalar(idx_or_list) else idx_or_list
        particles = Particles(
            icells  = self.iCells[idx_list].copy(),
            deltas  = self.deltas[idx_list].copy(),
            v       = self.v[idx_list].copy(),
            weights = self.weights[idx_list].copy(),
            charges = self.charges[idx_list].copy(),
            dl      = self.dl[idx_list].copy(),
        )
        self.erase(idx_list)
        return particles



    def split(self, sim): # REQUIRES C++ PYBIND PHARE LIB
        from pyphare.cpp import split_pyarrays_fn

        split_pyarrays = split_pyarrays_fn(sim.ndim, sim.interp_order, sim.refined_particle_nbr)(
          (self.iCells, self.deltas, self.weights, self.charges, self.v)
        )
        return Particles(
          icells=split_pyarrays[0],
          deltas=split_pyarrays[1],
          weights=split_pyarrays[2],
          charges=split_pyarrays[3],
          v=np.asarray(split_pyarrays[4]).reshape(int(len(split_pyarrays[4]) / 3), 3),
          dl = self.dl[0]/refinement_ratio + np.zeros((split_pyarrays[2].size,self.ndim))
        )


def all_assert(part1, part2):
    np.testing.assert_equal(part1.size(), part2.size())

    idx1 = np.argsort(part1.iCells + part1.deltas)
    idx2 = np.argsort(part2.iCells + part2.deltas)

    np.testing.assert_equal(len(idx1), len(idx2))

    np.testing.assert_array_equal(part1.iCells[idx1], part2.iCells[idx2])

    deltol = 1e-6 if any([part.deltas.dtype == np.float32 for part in [part1, part2]] ) else 1e-12
    np.testing.assert_allclose(part1.deltas[idx1], part2.deltas[idx2], atol=deltol)

    np.testing.assert_allclose(part1.v[idx1,0], part2.v[idx2,0], atol=1e-12)
    np.testing.assert_allclose(part1.v[idx1,1], part2.v[idx2,1], atol=1e-12)
    np.testing.assert_allclose(part1.v[idx1,2], part2.v[idx2,2], atol=1e-12)

    np.testing.assert_allclose(part1.dl[idx1], part2.dl[idx2], atol=1e-12)


def any_assert(part1, part2):
    np.testing.assert_equal(part1.size(), part2.size())
    assert part1.size() < 10 # slowness

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

    particles_out = copy(particles_in[0]) # use first, concat rest
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

    icells = np.delete(particles.iCells, idx)
    deltas = np.delete(particles.deltas, idx)
    vx = np.delete(particles.v[:,0], idx)
    vy = np.delete(particles.v[:,1], idx)
    vz = np.delete(particles.v[:,2], idx)
    v = np.zeros((vx.size, 3))
    v[:,0] = vx
    v[:,1] = vy
    v[:,2] = vz
    weights = np.delete(particles.weights, idx)
    dl = np.zeros((len(weights),particles.ndim))
    for i in range(particles.ndim):
        dl[:,i] = np.delete(particles.dl[:,i], idx)

    charges = np.delete(particles.charges, idx)
    return Particles(icells=icells,
                     deltas=deltas,
                     v = v,
                     weights=weights,
                     charges=charges,
                     dl = dl
                    )
