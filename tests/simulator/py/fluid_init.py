#!/usr/bin/env python3
#
# formatted with black

from tests.simulator.py import InitValueValidation

import phare.pharein as ph, numpy as np
from phare.pp.diagnostics import Diagnostics


class FluidInitValidation(InitValueValidation):

    pop_diag_file = lambda fn_name: "ions_pop_ions_${POP}_" + fn_name

    def test_1d(self):
        from phare.pp.diagnostics import _Fluid
        from tests.diagnostic import dump_all_diags

        diag_out_dir = "phare_outputs/fluid_init"
        dic = InitValueValidation.diag_options(diag_out_dir)
        dic.update(
            {
                "diags_fn": lambda model: dump_all_diags(model.populations),
                "populations": {
                    "beta": {"nbr_part_per_cell": 10000, "init": {"seed": 133337}},
                    "gamma": {"nbr_part_per_cell": 100000, "init": {"seed": 1333337}},
                },
            }
        )

        diags = self._simulate_diagnostics(dim=1, interp=1, input=dic)[_Fluid.__name__]
        self._checkFluidFromDensityFN(diags)
        self._checkFluidFromFluxFN(diags)
        self._checkIons(diags)

    def _sum_pop_patches(self, pop_patches, dataset):
        return np.sum(
            [
                self._aggregate_physical_patch_value(patches, True)[dataset]
                for pop, patches in pop_patches.items()
            ],
            axis=0,
        )

    def _checkIons(self, fluids):
        tolerance = 1e-6
        sim, model_dict = ph.globals.sim, ph.globals.sim.model.model_dict

        pop_densities, pop_fluxs = [
            self._get_contig_popdata_4_level0_ions_type(
                fluids, FluidInitValidation.pop_diag_file(fn_name=name)
            )
            for name in ["density", "flux"]
        ]

        sum_pop_density = self._sum_pop_patches(pop_densities, "density")
        ions_density = self._get_contig_data_4_level0_ions_type(fluids, "ions_density")[
            "density"
        ]
        np.testing.assert_allclose(
            sum_pop_density, ions_density, rtol=tolerance,
        )

        sum_pop_flux = {
            key: self._sum_pop_patches(pop_fluxs, key) for key in ["x", "y", "z"]
        }
        ions_bulkVelocity = self._get_contig_data_4_level0_ions_type(
            fluids, "ions_bulkVelocity"
        )
        for key, val in ions_bulkVelocity.items():
            np.testing.assert_allclose(
                sum_pop_flux[key] / sum_pop_density, val, rtol=tolerance
            )

    def _checkFluidFromFluxFN(self, fluids, fn_name="flux"):
        sim, model_dict = ph.globals.sim, ph.globals.sim.model.model_dict
        pops = self._get_contig_popdata_4_level0_ions_type(
            fluids, FluidInitValidation.pop_diag_file(fn_name)
        )

        def x0(pop, key, x):
            return model_dict[pop]["v" + key](x) * model_dict[pop]["density"](x)

        datasets, std_devs = self._get_std_devs_for(sim, pops, model_dict, x0)
        dev_per_ppc = self._sort_std_devs_per_ppc(pops, std_devs, model_dict)
        noises = self._get_noise_from(dev_per_ppc)
        nppc = [x[1] for x in dev_per_ppc]
        ratios = {key: noise[0] * np.sqrt(nppc[0]) for key, noise in noises.items()}
        self.assertTrue(
            all([diff < 0.09 for diff in self._get_diff(ratios, nppc, noises).values()])
        )
        for key, noise in noises.items():
            self._save_std_fig(
                nppc, noise, ratios[key], fn_name + key, test_threashold=0.007
            )
        self._save_dbg_fig(datasets, model_dict, fn_name)

    def _checkFluidFromDensityFN(self, fluids, fn_name="density"):
        sim, model_dict = ph.globals.sim, ph.globals.sim.model.model_dict
        pops = self._get_contig_popdata_4_level0_ions_type(
            fluids, FluidInitValidation.pop_diag_file(fn_name)
        )

        def x0(pop, key, x):
            return model_dict[pop][key](x)

        datasets, std_devs = self._get_std_devs_for(sim, pops, model_dict, x0)
        dev_per_ppc = self._sort_std_devs_per_ppc(pops, std_devs, model_dict)
        noises = self._get_noise_from(dev_per_ppc)
        nppc = [x[1] for x in dev_per_ppc]
        ratios = {key: noise[0] * np.sqrt(nppc[0]) for key, noise in noises.items()}
        self.assertTrue(self._get_diff(ratios, nppc, noises)[fn_name] < 0.09)
        self._save_std_fig(nppc, noises[fn_name], ratios[fn_name], fn_name)
        self._save_dbg_fig(datasets, model_dict, fn_name)


    def _get_contig_data_4_level0_ions_type(self, fluids, Type):
        diags = Diagnostics.get_type(fluids, Type)
        # can't be more than one, shouldn't be 0
        assert len(diags) == 1
        return self._aggregate_physical_patch_value(
            sorted(
                diags[0].levels[0].patches,
                key=lambda x: x.min_coord("x"),
            ),
            True,
        )

    def _get_contig_popdata_4_level0_ions_type(self, diags, Type):
        pop_patches = {}
        for diag in Diagnostics.get_type(diags, Type):
            for patch in diag.levels[0].patches:
                # pull "alpha" from hd5f path eg /t#/pl0/p0#0/ions/pop/ions_alpha/${fn_name}
                pop = patch.dtype.h5.name.split("/")[-2][5:]
                if pop not in pop_patches:
                    pop_patches[pop] = []
                pop_patches[pop].append(patch)

        return {
            pop: sorted(patches, key=lambda x: x.min_coord("x"))
            for pop, patches in pop_patches.items()
        }

    # this might make mores sense to live on PatchLevel and work per dim
    def _aggregate_physical_patch_value(self, patches, shared_patch_border=False):
        physical_datasets = {}
        for p in patches:
            plus = 0
            if shared_patch_border:
                # first value of second patch dataset is last of first
                plus = 0 if p == patches[-1] else 1
            for dataset in p.dtype.keys():
                nGhosts = p.dtype.nGhosts(dataset)
                if dataset not in physical_datasets:
                    physical_datasets[dataset] = np.asarray([])
                hdf5_data = p.dtype.get()[dataset]
                physical_datasets[dataset] = np.concatenate(
                    (
                        physical_datasets[dataset],
                        hdf5_data[nGhosts : int((nGhosts + plus) * -1)],
                    )
                )
        return physical_datasets

    def _get_std_devs_for(self, sim, pops, model_dict, lmbda):
        datasets, std_devs = {}, {}
        for pop, patches in pops.items():
            std_devs[pop] = {}
            cell_width = patches[0].patch_level.cell_width("x")
            datasets[pop] = self._aggregate_physical_patch_value(patches, True)
            for key in patches[0].dtype.keys():
                x = np.arange(sim.cells[0] + 1) * cell_width + sim.origin[0]
                std_devs[pop][key] = np.std(lmbda(pop, key, x) - datasets[pop][key])

        return datasets, std_devs

    def _sort_std_devs_per_ppc(self, pops, std_devs, model_dict):
        return sorted(
            [
                (std_devs[pop], model_dict[pop]["nbrParticlesPerCell"])
                for pop in list(pops.keys())
            ],
            key=lambda x: x[1],
        )

    def _get_noise_from(self, dev_per_ppc):
        noise = {key: [] for key in dev_per_ppc[0][0].keys()}
        for dev in dev_per_ppc:
            for k, v in dev[0].items():
                noise[k].append(v)
        return noise

    def _get_diff(self, ratio, nppc, noises):
        return {
            key: np.std(ratio[key] * 1.0 / np.sqrt(nppc) - noise)
            for key, noise in noises.items()
        }

    def _save_dbg_fig(self, pop_datasets, model_dict, fn_name):
        import matplotlib, matplotlib.pyplot as plt

        matplotlib.use("Agg")  # for systems without GUI

        colours = {1e2: "red", 1e3: "orange", 1e4: "blue", 1e5: "green"}

        dataset_keys = list(pop_datasets.values())[0].keys()

        for dataset_key in dataset_keys:
            for pop, datasets in pop_datasets.items():
                ppc = model_dict[pop]["nbrParticlesPerCell"]
                dataset = datasets[dataset_key]
                plt.plot(
                    np.asarray(dataset), label=str(ppc) + " ppc", color=colours[ppc]
                )

            plt.legend()
            plt.title("Population Flux from FN: " + fn_name + dataset_key)
            plt.savefig(fn_name + dataset_key + "_data.png", bbox_inches="tight")
            plt.clf()  # clear figure

    def _save_std_fig(self, nppc, noise, ratio, fig_name, test_threashold=0.001):
        import matplotlib, matplotlib.pyplot as plt

        matplotlib.use("Agg")  # for systems without GUI

        plt.plot(
            nppc,
            np.abs(ratio * 1.0 / np.sqrt(nppc) - noise),
            marker="o",
            label="theory - measured (noise)",
        )
        plt.axhline(test_threashold, label="test threshold", color="red")
        plt.axhline(np.std(ratio * 1.0 / np.sqrt(nppc) - noise), label="std")
        plt.legend()
        plt.title("Population density noise: theory vs measure")
        plt.xlabel("# Particles Per Cell")
        plt.savefig(fig_name + ".png", bbox_inches="tight")
        plt.clf()  # clear figure


if __name__ == "__main__":
    import unittest

    unittest.main()
