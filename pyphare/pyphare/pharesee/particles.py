
import numpy as np

class Particles():
    """
    this class represent a set of particles
    particles can either be loaded randomly in a given box or have there attribute set from caller
    """
    def __init__(self, **kwargs):
        if "box" in kwargs:
            box = kwargs["box"]
            self.iCells = np.random.randint(box.lower, high=box.upper+1, size=box.size()*100)
            self.deltas = np.random.rand(box.size()*100)
            self.vx = np.random.randn(box.size()*100)
            self.vz = np.random.randn(box.size()*100)
            self.vy = np.random.randn(box.size()*100)
        else:
            self.iCells = kwargs["icell"]
            self.deltas = kwargs["deltas"]
            self.vx = kwargs["vx"]
            self.vy = kwargs["vy"]
            self.vz = kwargs["vz"]

    def select(self, box):
        """
        select particles from the given box
        assumption, box has AMR indexes of the same level as the data that the current instance is created from
        """
        idx = np.where((self.iCells >= box.lower) & (self.iCells <= box.upper))[0]
        return Particles(icell=self.iCells[idx],
                         deltas=self.deltas[idx],
                         vx=self.vx[idx],
                         vy=self.vy[idx],
                         vz=self.vz[idx])


    def add(self, particles):

        self.iCells = np.concatenate((self.iCells, particles.iCells))
        self.deltas = np.concatenate((self.deltas, particles.deltas))
        self.vx = np.concatenate((self.vx, particles.vx))
        self.vy = np.concatenate((self.vy, particles.vy))
        self.vz = np.concatenate((self.vz, particles.vz))


    def shift_icell(self, offset):
        self.iCells += offset
        return self
