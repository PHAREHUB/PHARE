
import numpy as np

class Particles():
    def __init__(self, **kwargs):
        if "box" in kwargs:
            box = kwargs["box"]
            self.iCells = np.arange(box.lower, box.upper + 1)
            self.deltas = np.random.rand(box.size(), 100)
            self.vx = np.random.randn(box.size(), 100)
            self.vz = np.random.randn(box.size(), 100)
            self.vy = np.random.randn(box.size(), 100)
        else:
            self.iCells = kwargs["icell"]
            self.deltas = kwargs["deltas"]
            self.vx = kwargs["vx"]
            self.vy = kwargs["vy"]
            self.vz = kwargs["vz"]


