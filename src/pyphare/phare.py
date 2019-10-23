
import sys, inspect

class Diagnostic:

    def __init__(
          self,
          name : str,
          species : str,
          type : str,
          compute_every = 1,
          write_every = 1,
          start_iteration = 0,
          end_iteration = sys.maxsize):
        _, _, _, kvs = inspect.getargvalues(inspect.currentframe())
        del kvs["self"]
        for k, v in kvs.items():
            object.__setattr__(self, k, v)

    def add_to_simulation(self, pp):
        pp.addDiagnostic(
            self.compute_every,
            self.write_every,
            self.start_iteration,
            self.end_iteration,
            self.name,
            self.species,
            self.type
        )
