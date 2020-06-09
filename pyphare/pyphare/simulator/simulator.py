

from pybindlibs import cpp
import pyphare.pharein as ph



def create_simulator():

    cpp.reset()
    ph.populateDict()

    import sys

    try:
        hier = cpp.make_hierarchy()
        sim = cpp.make_simulator(hier)
        return [cpp.make_diagnostic_manager(sim, hier), sim, hier]
    except:
        e = sys.exc_info()[0]
        print('Exception caught in "create_simulator": {}'.format(e))
