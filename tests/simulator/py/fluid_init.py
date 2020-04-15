#!/usr/bin/env python3
#
# formatted with black

from tests.simulator.py import InitValueValidation

import phare.pharein as ph, numpy as np
from phare.pp.diagnostics import Diagnostics, selectDiag, selectDiags
from phare.pp.diagnostic.patch import aggregate_level0_patch_domain


class FluidInitValidation(InitValueValidation):
    min_interp = 1
    max_interp = 1

    def pop_diag_file(fn_name):
        return "ions_pop_ions_${POP}_" + fn_name

    def test_1d(self):
        from phare.pp.diagnostics import _FluidPatchData
        from tests.diagnostic import dump_all_diags

        dim = 1

        min_interp = FluidInitValidation.min_interp
        max_interp = FluidInitValidation.max_interp + 1
        for interp in range(min_interp, max_interp):
            diag_out_dir = "phare_outputs/fluid_init" + str(dim) + "_" + str(interp)
            simInputs = InitValueValidation.diag_options(diag_out_dir)
            simInputs.update(
                {
                    "diags_fn": lambda model: dump_all_diags(model.populations),
                    "populations": {
                        "beta": {"nbr_part_per_cell": 10000, "init": {"seed": 133337}},
                        "gamma": {
                            "nbr_part_per_cell": 100000,
                            "init": {"seed": 1333337},
                        },
                    },
                }
            )

            diags = self.runAndDump(dim=1, interp=interp, input=simInputs)[
                _FluidPatchData.__name__
            ]
            self.checkPopulationDensityEqualsUserFunction(diags, dim, interp)
            self.checkPopulationFluxEqualsUserFunction(diags, dim, interp)
            self.checkIonDensityAndVelocityConsistentWithPopulations(diags, dim, interp)

    def checkIonDensityAndVelocityConsistentWithPopulations(
        self, fluid_diags, dim, interp
    ):
        tolerance = 1e-6
        sim = ph.globals.sim
        model_dict = ph.globals.sim.model.model_dict

        pop_densities, pop_fluxes = [
            _pop_name_to_patch_level0(
                fluid_diags, FluidInitValidation.pop_diag_file(fn_name=name)
            )
            for name in ["density", "flux"]
        ]

        total_density = _sum_pop_patches(pop_densities, "density")
        ions_density_diag = selectDiag(fluid_diags, "ions_density")
        ions_density = aggregate_level0_patch_domain(
            ions_density_diag.levels[0], shared_patch_border=True
        )["density"]

        np.testing.assert_allclose(total_density, ions_density, atol=tolerance)

        total_flux = {key: _sum_pop_patches(pop_fluxes, key) for key in ["x", "y", "z"]}
        ions_bulkVelocity_diag = selectDiag(fluid_diags, "ions_bulkVelocity")
        ions_bulkVelocity = aggregate_level0_patch_domain(
            ions_bulkVelocity_diag.levels[0], shared_patch_border=True
        )

        for vel, flux in ions_bulkVelocity.items():
            np.testing.assert_allclose(
                total_flux[vel] / total_density, flux, atol=tolerance
            )

    def checkPopulationFluxEqualsUserFunction(self, fluid_diags, dim, interp):
        threshold = 0.007
        fn_name = "flux"
        sim = ph.globals.sim
        model_dict = ph.globals.sim.model.model_dict

        pop_fluxes = _pop_name_to_patch_level0(
            fluid_diags, FluidInitValidation.pop_diag_file(fn_name)
        )

        def user_flux(pop_name, component_name, x):
            return model_dict[pop_name]["v" + component_name](x) * model_dict[pop_name][
                "density"
            ](x)

        datasets, std_devs = _get_std_devs_for(sim, pop_fluxes, model_dict, user_flux)
        dev_per_ppc = _sort_std_devs_per_ppc(pop_fluxes, std_devs, model_dict)
        noises = _get_noise_from(dev_per_ppc)
        nppc = [x[1] for x in dev_per_ppc]
        ratios = {key: noise[0] * np.sqrt(nppc[0]) for key, noise in noises.items()}
        self.assertTrue(
            all([diff < 0.09 for diff in _get_diff(ratios, nppc, noises).values()])
        )
        for key, noise in noises.items():
            _save_std_fig(
                nppc,
                noise,
                ratios[key],
                fn_name + key,
                dim,
                interp,
                test_threshold=threshold,
            )
        _save_dbg_fig(datasets, model_dict, fn_name, dim, interp)

    def checkPopulationDensityEqualsUserFunction(self, fluid_diags, dim, interp):
        threshold = 0.09
        fn_name = "density"
        sim = ph.globals.sim
        model_dict = ph.globals.sim.model.model_dict

        pop_densities = _pop_name_to_patch_level0(
            fluid_diags, FluidInitValidation.pop_diag_file(fn_name)
        )

        def user_density(pop_name, component_name, x):
            return model_dict[pop_name][component_name](x)

        datasets, std_devs = _get_std_devs_for(
            sim, pop_densities, model_dict, user_density
        )
        dev_per_ppc = _sort_std_devs_per_ppc(pop_densities, std_devs, model_dict)
        noises = _get_noise_from(dev_per_ppc)
        nppc = [x[1] for x in dev_per_ppc]
        ratios = {key: noise[0] * np.sqrt(nppc[0]) for key, noise in noises.items()}
        self.assertTrue(_get_diff(ratios, nppc, noises)[fn_name] < threshold)
        _save_std_fig(nppc, noises[fn_name], ratios[fn_name], fn_name, dim, interp)
        _save_dbg_fig(datasets, model_dict, "", dim, interp)  # fn_name == dataset key


def _pop_name_to_patch_level0(diags, Type):
    pop_patch_level = {}
    for diag in selectDiags(diags, Type):
        patches = diag.levels[0].patches_list()
        assert len(patches)
        # pull "alpha" from hd5f path eg /t#/pl0/p0#0/ions/pop/ions_alpha/${fn_name}
        pop_name = patches[0].patch_data.h5item.name.split("/")[-2][5:]
        pop_patch_level[pop_name] = diag.levels[0]

    return pop_patch_level


def _get_std_devs_for(sim, pop_patch_levels, model_dict, lmbda):
    datasets, std_devs = {}, {}
    for pop_name, patch_level in pop_patch_levels.items():
        patches = patch_level.patches_list()
        std_devs[pop_name] = {}
        cell_width = patch_level.cell_width("x")
        datasets[pop_name] = aggregate_level0_patch_domain(
            patch_level, shared_patch_border=True
        )
        for data_name in patch_level.data_names():
            x = np.arange(sim.cells[0] + 1) * cell_width + sim.origin[0]
            std_devs[pop_name][data_name] = np.std(
                lmbda(pop_name, data_name, x) - datasets[pop_name][data_name]
            )

    return datasets, std_devs


def _sort_std_devs_per_ppc(pop_patch_levels, std_devs, model_dict):
    return sorted(
        [
            (std_devs[pop_name], model_dict[pop_name]["nbrParticlesPerCell"])
            for pop_name in list(pop_patch_levels.keys())
        ],
        key=lambda tple: tple[1],
    )


def _get_noise_from(devs_per_ppc):
    """
      devs_per_ppc = list[(data_name: std_dev, n_ppc)] # sorted per ppc
    """
    noises = {}
    for data_name_dev_dict, _ in devs_per_ppc:
        for data_name, std_dev in data_name_dev_dict.items():
            if data_name not in noises:
                noises[data_name] = []
            noises[data_name].append(std_dev)
    return noises  # = data_name:list[std_dev] # list sorted by ppc ascending


def _get_diff(ratio, nppc, noises):
    return {
        data_name: np.std(ratio[data_name] * 1.0 / np.sqrt(nppc) - noise)
        for data_name, noise in noises.items()
    }


def _sum_pop_patches(pop_patch_levels, dataset_name):
    return np.sum(
        [
            aggregate_level0_patch_domain(patch_level, shared_patch_border=True)[
                dataset_name
            ]
            for patch_level in list(pop_patch_levels.values())
        ],
        axis=0,
    )


def _save_dbg_fig(pop_datasets, model_dict, fn_name, dim, interp):
    import matplotlib, matplotlib.pyplot as plt

    matplotlib.use("Agg")  # for systems without GUI

    colours = {1e2: "red", 1e3: "orange", 1e4: "blue", 1e5: "green"}

    for dataset_key in list(pop_datasets.values())[0].keys():
        for pop, datasets in pop_datasets.items():
            ppc = model_dict[pop]["nbrParticlesPerCell"]
            dataset = datasets[dataset_key]
            plt.plot(np.asarray(dataset), label=str(ppc) + " ppc", color=colours[ppc])

        plt.legend()
        plt.title("Population Flux from FN: " + fn_name + dataset_key)
        out_png = (
            fn_name + dataset_key + "_" + str(dim) + "_" + str(interp) + "_dbg.png"
        )
        plt.savefig(out_png, bbox_inches="tight")
        plt.clf()  # clear figure


def _save_std_fig(nppc, noise, ratio, fig_name, dim, interp, test_threshold=0.001):
    import matplotlib, matplotlib.pyplot as plt

    matplotlib.use("Agg")  # for systems without GUI

    plt.plot(
        nppc,
        np.abs(ratio * 1.0 / np.sqrt(nppc) - noise),
        marker="o",
        label="theory - measured (noise)",
    )
    plt.axhline(test_threshold, label="test threshold", color="red")
    plt.axhline(np.std(ratio * 1.0 / np.sqrt(nppc) - noise), label="std")
    plt.legend()
    plt.title("Population density noise: theory vs measure")
    plt.xlabel("# Particles Per Cell")
    out_png = fig_name + "_" + str(dim) + "_" + str(interp) + ".png"
    plt.savefig(out_png, bbox_inches="tight")
    plt.clf()  # clear figure


if __name__ == "__main__":
    import unittest

    unittest.main()
