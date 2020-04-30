
import numpy as np

class Particles:
    def __init__(self, **kwargs):
        if "box" in kwargs:
            box = kwargs["box"]

            self.iCells  = np.random.randint(box.lower, high=box.upper+1, size=box.size()*100)
            self.deltas  = np.random.rand(box.size()*100)
            self.v       = np.random.randn(box.size()*100, 3)
            self.weights = np.zeros_like(self.deltas) + 0.01
            self.charges = np.zeros_like(self.weights) + 1

        else:
            self.iCells  = kwargs["icell"]
            self.deltas  = kwargs["deltas"]
            self.v       = kwargs["v"]
            self.weights = kwargs["weights"]
            self.charges = kwargs["charges"]


    def select(self, box):
        idx = np.where((self.iCells >= box.lower) & (self.iCells <= box.upper))[0]
        return Particles(icell=self.iCells[idx],
                         deltas=self.deltas[idx],
                         v=self.v[idx,:],
                         weights = self.weights[idx],
                         charges = self.charges[idx])


    def add(self, particles):

        self.iCells   = np.concatenate((self.iCells, particles.iCells))
        self.deltas   = np.concatenate((self.deltas, particles.deltas))
        self.v        = np.concatenate((self.v, particles.v))
        self.charges  = np.concatenate((self.charges, particles.charges))
        self.weights  = np.concatenate((self.weights, particles.weights))


    def shift_icell(self, offset):
        self.iCells += offset
        return self

