
import numpy as np


class Particles:
    """
    this class represent a set of particles
    particles can either be loaded randomly in a given box or have there attribute set from caller
    """
    def __init__(self, **kwargs):
        if "box" in kwargs:
            box = kwargs["box"]

            self.iCells  = np.random.randint(box.lower, high=box.upper+1, size=box.size()*100)
            self.deltas  = np.random.rand(box.size()*100)
            self.v       = np.random.randn(box.size()*100, 3)
            self.weights = np.zeros_like(self.deltas) + 0.01
            self.charges = np.zeros_like(self.weights) + 1

        else:
            self.iCells  = kwargs["icells"]
            self.deltas  = kwargs["deltas"]
            self.v       = kwargs["v"]
            self.weights = kwargs["weights"]
            self.charges = kwargs["charges"]



    def add(self, particles):
        self.iCells   = np.concatenate((self.iCells, particles.iCells))
        self.deltas   = np.concatenate((self.deltas, particles.deltas))
        self.v        = np.concatenate((self.v, particles.v))
        self.charges  = np.concatenate((self.charges, particles.charges))
        self.weights  = np.concatenate((self.weights, particles.weights))


    def shift_icell(self, offset):
        self.iCells += offset
        return self

    def size(self):
          return len(self.weights)

    def select(self, box):
        """
        select particles from the given box
        assumption, box has AMR indexes of the same level as the data that the current instance is created from
        """
        idx = np.where((self.iCells >= box.lower) & (self.iCells <= box.upper))[0]
        return Particles(icells=self.iCells[idx],
                         deltas=self.deltas[idx],
                         v = self.v[idx,:],
                         weights=self.weights[idx],
                         charges=self.charges[idx])


    def split(self, sim): # REQUIRES C++ PYBIND PHARE LIB
        from pyphare.cpp import split_pyarrays_fn

        split_pyarrays = split_pyarrays_fn(sim.dims, sim.interp_order, sim.refined_particle_nbr)(
          (self.iCells, self.deltas, self.weights, self.charges, self.v)
        )
        return Particles(
          icells=split_pyarrays[0],
          deltas=split_pyarrays[1],
          weights=split_pyarrays[2],
          charges=split_pyarrays[3],
          v=np.asarray(split_pyarrays[4]).reshape(int(len(split_pyarrays[4]) / 3), 3)
        )
