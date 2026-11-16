from ..fully_kinetic import check_fully_kinetic
from .dict_populator import DictPopulator
from .hybrid import _populate_electromagnetic, _populate_populations


def populateDict(sim):
    """all species, electrons included, are populations of the MaxwellianFluidModel;
    there is no ElectronModel and no Ohm's law"""
    check_fully_kinetic(sim)
    dp = DictPopulator()

    dp.add_bool("simulation/fully_kinetic", True)
    dp.add_double("simulation/c_over_va", sim.c_over_va)
    dp.add_int("simulation/current_filter_passes", sim.current_filter_passes)
    dp.add_string("simulation/algo/ion_updater/pusher/name", sim.particle_pusher)

    _populate_populations(dp, sim)
    _populate_electromagnetic(dp, sim)
