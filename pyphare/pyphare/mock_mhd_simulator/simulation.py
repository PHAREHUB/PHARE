from . import global_vars


class Simulation(object):
    def __init__(self, **kwargs):
        if global_vars.sim is not None:
            raise RuntimeError("simulation is already created")

        global_vars.sim = self

        for k, v in kwargs.items():
            object.__setattr__(self, k, v)

        self.model = None

    def set_model(self, model):
        self.model = model

    def clear(self):
        global_vars.sim = None
