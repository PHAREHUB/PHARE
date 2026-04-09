#
#
#

from pyphare.core import phare_utilities as phut
from pyphare.pharesee.hierarchy import hierarchy_compute as hc

from pyphare.pharesee import hierarchy as hierm

from . import utils as uzi


class LazyFieldData(hierm.patchdata.FieldData):
    def __init__(self, pd, default_mutators, select_mutators):
        dic = dict(ghosts_nbr=pd.ghosts_nbr, centering=pd.centerings)
        super().__init__(pd.layout, pd.name, pd.dataset, **dic)

        self.default_mutators = phut.listify(default_mutators or [])
        self.select_mutators = phut.listify(select_mutators or [])

    def update(self, pd):
        self.dataset = pd.dataset
        self.ghosts_nbr = pd.ghosts_nbr
        self.centerings = pd.centerings
        return self

    def __getitem__(self, get):
        default_mutators, select_mutators = self.default_mutators, self.select_mutators
        self.default_mutators = []  # only compute once!
        self.select_mutators = []

        for mut in default_mutators:
            self.update(mut(self))
        for mut in select_mutators:
            self.update(mut(self.select(get)))

        return self.select(get)

    def copy_as(self, data=None, **kwargs):
        # if self.default_mutators or self.select_mutators:
        return LazyFieldData(
            super().copy_as(data=data, **kwargs),
            self.default_mutators,
            self.select_mutators,
        )
        # return super().copy_as(data=data, **kwargs)


class LazyPatch(hierm.patch.Patch):
    def __init__(self, *args, **kwargs):
        self.calculators = kwargs.pop("calculators", {}) or {}
        super().__init__(*args, **kwargs)
        self.lazy_patch_datas = {}

    def __getitem__(self, key):
        if key in self.patch_datas:
            return self.patch_datas[key]
        if key not in self.lazy_patch_datas:
            if key in self.calculators:
                self.lazy_patch_datas.update(self.calculators[key](self))
            else:
                for k in self.patch_datas.keys():
                    if k.endswith(key):
                        return self.patch_datas[k]
        return self.lazy_patch_datas[key]


def lazy_patch_from_base(patch, calculators):
    return LazyPatch(
        patch.patch_datas, patch.id, patch.layout, patch.attrs, calculators=calculators
    )


def resolve_selected_level_patches(time, hier, level, ilvl):
    if not hier.selection_box[time]:
        return level.patches

    overlapped = []
    boxes = phut.listify(hier.selection_box[time] or [])
    for patch in level:
        for box in boxes:
            level_box = hier.refined_box(ilvl, box)
            if patch.box * level_box is not None:
                overlapped.append(patch)
                break

    return overlapped


def lazy_patch_hierarchy(
    hier, default_mutators=None, select_mutators=None, calculators=None
):
    for time, levels in hier.time_hier.items():
        for ilvl, lvl in levels.items():
            lvl.patches = [
                lazy_patch_from_base(patch, calculators) for patch in lvl.patches
            ]
            for patch in resolve_selected_level_patches(time, hier, lvl, ilvl):
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
        for key, val in kwargs.items():
            setattr(self, key, val)

    def RawB(self, time, **kwargs):
        return self.run._get_hier_for(time, "EM_B", **kwargs)

    def GetB(self, time, **kwargs):
        return hierm.VectorField.FROM(
            lazy_patch_hierarchy(self.RawB(time, **kwargs), uzi._compute_to_primal)
        )

    def RawVi(self, time, **kwargs):
        return self.run._get_hier_for(time, "ions_bulkVelocity", **kwargs)

    def GetVi(self, time, **kwargs):
        return hierm.VectorField.FROM(
            lazy_patch_hierarchy(self.RawVi(time, **kwargs), hc.drop_ghosts)
        )

    def RawN(self, time, pop_name, **kwargs):
        return self.run._get_hier_for(time, f"ions_pop_{pop_name}_density", **kwargs)

    def GetN(self, time, pop_name, **kwargs):
        return hierm.ScalarField.FROM(
            lazy_patch_hierarchy(self.RawN(time, pop_name, **kwargs), hc.drop_ghosts)
        )

    def RawFlux(self, time, pop_name, **kwargs):
        return self.run._get_hier_for(time, f"ions_pop_{pop_name}_flux", **kwargs)

    def GetFlux(self, time, pop_name, **kwargs):
        return hierm.VectorField.FROM(
            lazy_patch_hierarchy(self.RawFlux(time, pop_name, **kwargs), hc.drop_ghosts)
        )

    def GetPressure(self, time, pop_name, **kwargs):
        hier = self.run._get_hierarchy(
            time, f"ions_pop_{pop_name}_momentum_tensor.h5", **kwargs
        )

        hier = self.RawFlux(time, pop_name, hier=hier)
        hier = self.RawN(time, pop_name, hier=hier)
        pressure_kwargs = dict(
            mass=self.run.GetMass(pop_name, **kwargs), popname=pop_name
        )
        pxx = RunFunc(uzi._compute_pop_pressure_xx, **pressure_kwargs)
        return hierm.TensorField.FROM(
            lazy_patch_hierarchy(
                hier,
                hc.drop_ghosts,
                calculators={f"{pop_name}_Pxx": pxx},
            ),
        )
