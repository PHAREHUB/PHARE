#
#
#

from pyphare.core import phare_utilities as phut
from pyphare.pharesee.hierarchy import hierarchy_compute as hc

from pyphare.pharesee import hierarchy as hierm

#  import patch, patchdata
# from pyphare.pharesee.hierarchy import ScalarField, VectorField, TensorField

# from pyphare.pharesee.hierarchy.hierarchy_utils import compute_hier_from

from . import utils as uzi


class LazyFieldData(hierm.patchdata.FieldData):
    def __init__(self, pd, default_mutators, select_mutators):
        dic = dict(ghosts_nbr=pd.ghosts_nbr, centering=pd.centerings)
        super().__init__(pd.layout, pd.name, pd.dataset, **dic)

        self.default_mutators = phut.listify(default_mutators) or []
        self.select_mutators = phut.listify(select_mutators) or []

    def update(self, pd):
        self.dataset = pd.dataset
        self.ghosts_nbr = pd.ghosts_nbr
        self.centerings = pd.centerings

    def __getitem__(self, get):
        if self.default_mutators:
            get_all = get == slice(None) or get is None
            if get_all:
                for mut in self.default_mutators:
                    self.update(mut(self))
            for mut in self.select_mutators:
                self.update(mut(self))
            self.default_mutators = []
            self.select_mutators = []
        return self.select(get)


class LazyPatch(hierm.patch.Patch):
    def __init__(self, patch, calculators):
        super().__init__(patch.patch_datas, patch.id, patch.layout, patch.attrs)
        self.pds = {}
        self.calculators = calculators or {}
        if not type(self.calculators) is dict:
            raise RuntimeError(
                "Lazy Patch calculators is expected to be a dict[str, callable]"
            )

    def __getitem__(self, key):
        if key in self.patch_datas:
            return self.patch_datas[key]
        print("self.calculators", self.calculators, key)
        if key not in self.pds:
            assert key in self.calculators
            if key in self.calculators:
                self.pds.update(self.calculators[key](self))
        return self.pds[key]


def lazy_patch_hierarchy(
    hier, default_mutators=None, select_mutators=None, calculators=None
):
    for time, levels in hier.time_hier.items():
        for ilvl, lvl in levels.items():
            lvl.patches = [LazyPatch(patch, calculators) for patch in lvl.patches]
            for patch in lvl.patches:
                patch.patch_datas = {
                    name: LazyFieldData(pd, default_mutators, select_mutators)
                    for name, pd in patch
                }

    return hier


class RunFunc:
    def __init__(self, λ, *args, **kwargs):
        self.λ = λ
        self.args = args
        self.kwargs = {**kwargs}

    def __call__(self, *args, **kwargs):
        return self.λ(*args, *self.args, **kwargs, **self.kwargs)


class RunMan:
    def __init__(self, run, *args, **kwargs):
        self.run = run
        self.patch_mutators = []
        self.level_mutators = []
        for key, val in kwargs.items():
            setattr(self, key, val)

    def RawVi(self, time, **kwargs):
        return self.run._get_hier_for(time, "ions_bulkVelocity")

    def GetVi(self, time, merged=False, interp="nearest", **kwargs):
        return lazy_patch_hierarchy(self.RawVi(time, **kwargs), hc.drop_ghosts)

    def RawN(self, time, pop_name, **kwargs):
        return self.run._get_hier_for(time, f"ions_pop_{pop_name}_density", **kwargs)

    def GetN(self, time, pop_name, merged=False, interp="nearest", **kwargs):
        return hierm.ScalarField(
            lazy_patch_hierarchy(self.RawN(time, pop_name, **kwargs), hc.drop_ghosts)
        )

    def RawFlux(self, time, pop_name, **kwargs):
        return self.run._get_hier_for(time, f"ions_pop_{pop_name}_flux", **kwargs)

    def GetFlux(self, time, pop_name, merged=False, interp="nearest", **kwargs):
        return hierm.VectorField(
            lazy_patch_hierarchy(self.RawFlux(time, pop_name, **kwargs), hc.drop_ghosts)
        )

    def GetPressure(self, time, pop_name, merged=False, interp="nearest", **kwargs):
        hier = self.run._get_hierarchy(
            time, f"ions_pop_{pop_name}_momentum_tensor.h5", **kwargs
        )

        hier = self.RawFlux(time, pop_name, hier=hier)
        hier = self.RawN(time, pop_name, hier=hier)
        pressure_kwargs = dict(
            mass=self.run.GetMass(pop_name, **kwargs), popname=pop_name
        )
        pxx = RunFunc(uzi._compute_pop_pressure_xx, **pressure_kwargs)
        return hierm.TensorField(  # CAN'T USE RENAME!
            lazy_patch_hierarchy(
                hier,
                hc.drop_ghosts,
                calculators={f"{pop_name}_Pxx": pxx},
            ),
        )
