# only one simulator is possible because of statics in SAMRAI anyway
#  weakref.ref to the simulator that was last set up, or None
obj = None


def live():
    return obj() if obj is not None else None
