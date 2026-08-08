#ifndef PHARE_DIAGNOSTIC_COMPUTERS_DIAGNOSTIC_COMPUTERS_HPP_
#define PHARE_DIAGNOSTIC_COMPUTERS_DIAGNOSTIC_COMPUTERS_HPP_

#include "core/data/tensorfield/tensorfield.hpp"
#include "core/numerics/interpolator/interpolator.hpp"
#include "core/numerics/primite_conservative_converter/to_primitive_converter.hpp"

#include <algorithm>
#include <functional>

namespace PHARE::diagnostic
{

// PUBLIC API - writer-agnostic: only ever needs the model view and the writer's current
// level range/timestamp, none of which are specific to any particular writer type (h5, vtk)
void compute_momentum_tensor(auto& modelView, std::size_t minLvl, std::size_t maxLvl,
                             double timestamp);
void compute_pop_momentum_tensor(auto& modelView, auto& pop, std::size_t minLvl,
                                 std::size_t maxLvl, double timestamp);
void compute_kinetic_energy_flux_vector(auto& modelView, std::size_t minLvl, std::size_t maxLvl,
                                        double timestamp);
void compute_pop_kinetic_energy_flux_vector(auto& modelView, auto& pop, std::size_t minLvl,
                                            std::size_t maxLvl, double timestamp);
void compute_mhd_velocity(auto& modelView, std::size_t minLvl, std::size_t maxLvl,
                          double timestamp);
void compute_mhd_pressure(auto& modelView, double gamma, std::size_t minLvl, std::size_t maxLvl,
                          double timestamp);



// PRIVATE IMPL


void _fill_all_pop_schedules(auto& modelView, auto fn, std::size_t minLvl, std::size_t maxLvl,
                             double timestamp)
{
    auto const fill_schedules = [&](auto& lvl) {
        for (std::size_t i = 0; i < modelView.getIons().size(); ++i)
            fn(modelView, lvl, timestamp, i);
    };
    modelView.onLevels(fill_schedules, minLvl, maxLvl);
}

void _fill_one_pop_schedule(auto& modelView, auto fn, auto& pop, std::size_t minLvl,
                            std::size_t maxLvl, double timestamp)
{
    auto& ions                = modelView.getIons();
    auto const i              = ions.pop_index(pop.name());
    auto const fill_schedules = [&](auto& lvl) { fn(modelView, lvl, timestamp, i); };
    modelView.onLevels(fill_schedules, minLvl, maxLvl);
}

// sums resources[0..n) (per population) into resources[n] (the ions-level aggregate),
// component-wise, where n == ions.size()
template<std::size_t n_components>
void _sum_pop_resources_into_aggregate(auto& modelView, auto&& getResource, std::size_t n_pops)
{
    auto& total = getResource(n_pops);
    total.zero();
    for (std::size_t i = 0; i < n_pops; ++i)
    {
        auto& part = getResource(i);
        for (std::uint8_t c = 0; c < n_components; ++c)
            std::transform(part[c].begin(), part[c].end(), total[c].begin(), total[c].begin(),
                           std::plus<>{});
    }
}


template<typename ModelView_t>
struct KineticEnergyFluxVectorComputer
{
    using Model_t      = ModelView_t::Model_t;
    using level_t      = Model_t::amr_types::level_t;
    using GridLayout_t = Model_t::gridlayout_type;
    auto constexpr static fillFn
        = &ModelView_t::template fillPopKineticEnergyFluxVector<level_t, double, std::size_t>;
    auto constexpr static dimension    = GridLayout_t::dimension;
    auto constexpr static interp_order = GridLayout_t::options.interp_order;
    auto constexpr static n_components = core::detail::tensor_field_dim_from_rank<1>();

    void fillAllSchedules(std::size_t minLvl, std::size_t maxLvl, double timestamp)
    {
        _fill_all_pop_schedules(modelView, std::mem_fn(fillFn), minLvl, maxLvl, timestamp);
    }

    void fillPopSchedules(auto& pop, std::size_t minLvl, std::size_t maxLvl, double timestamp)
    {
        _fill_one_pop_schedule(modelView, std::mem_fn(fillFn), pop, minLvl, maxLvl, timestamp);
    }

    auto interpolate(auto& pop, std::size_t popidx, auto& layout)
    {
        auto& kineticEnergyFlux = modelView.tmpVecField(popidx);
        kineticEnergyFlux.zero();
        interpolator(pop, kineticEnergyFlux, pop.domainParticles(), layout);
        interpolator(pop, kineticEnergyFlux, pop.levelGhostParticlesOld(), layout);
    };

    ModelView_t& modelView;
    core::KineticEnergyFluxVectorInterpolator<dimension, interp_order> interpolator{};
};
template<typename ModelView_t>
KineticEnergyFluxVectorComputer(ModelView_t&) -> KineticEnergyFluxVectorComputer<ModelView_t>;

void compute_kinetic_energy_flux_vector(auto& modelView, std::size_t minLvl, std::size_t maxLvl,
                                        double timestamp)
{
    KineticEnergyFluxVectorComputer computer{modelView};
    auto& ions             = modelView.getIons();
    auto const n_pops      = ions.size();
    auto const interpolate = [&](auto& layout, auto&&...) {
        std::size_t i = 0;
        for (auto& pop : ions)
            computer.interpolate(pop, i++, layout);
    };
    modelView.visitHierarchy(interpolate, minLvl, maxLvl);
    computer.fillAllSchedules(minLvl, maxLvl, timestamp);
    modelView.visitHierarchy(
        [&](auto&&...) {
            _sum_pop_resources_into_aggregate<decltype(computer)::n_components>(
                modelView, [&](std::size_t idx) -> auto& { return modelView.tmpVecField(idx); },
                n_pops);
        },
        minLvl, maxLvl);
}


void compute_pop_kinetic_energy_flux_vector(auto& modelView, auto& pop, std::size_t minLvl,
                                            std::size_t maxLvl, double timestamp)
{
    KineticEnergyFluxVectorComputer computer{modelView};
    auto const i           = modelView.getIons().pop_index(pop.name());
    auto const interpolate = [&](auto& layout, auto&&...) { computer.interpolate(pop, i, layout); };
    modelView.visitHierarchy(interpolate, minLvl, maxLvl);
    computer.fillPopSchedules(pop, minLvl, maxLvl, timestamp);
}




template<typename ModelView_t>
struct MomentumTensorComputer
{
    // compute the momentum tensor for each population that requires it
    // compute for all ions but that requires the computation of all pop

    // dumps occur after the last substep but before the next first substep
    // at this time, levelGhostPartsNew is emptied and not yet filled
    // and the former levelGhostPartsNew has been moved to levelGhostPartsOld

    using Model_t      = ModelView_t::Model_t;
    using level_t      = Model_t::amr_types::level_t;
    using GridLayout_t = Model_t::gridlayout_type;
    auto constexpr static fillFn
        = &ModelView_t::template fillPopMomTensor<level_t, double, std::size_t>;
    auto constexpr static dimension    = GridLayout_t::dimension;
    auto constexpr static interp_order = GridLayout_t::options.interp_order;
    auto constexpr static n_components = core::detail::tensor_field_dim_from_rank<2>();

    void fillAllSchedules(std::size_t minLvl, std::size_t maxLvl, double timestamp)
    {
        _fill_all_pop_schedules(modelView, std::mem_fn(fillFn), minLvl, maxLvl, timestamp);
    }

    void fillPopSchedules(auto& pop, std::size_t minLvl, std::size_t maxLvl, double timestamp)
    {
        _fill_one_pop_schedule(modelView, std::mem_fn(fillFn), pop, minLvl, maxLvl, timestamp);
    }

    auto interpolate(auto& pop, std::size_t popidx, auto& layout)
    {
        auto& pop_momentum_tensor = modelView.tmpTensorField(popidx);
        pop_momentum_tensor.zero();
        interpolator(pop.domainParticles(), pop_momentum_tensor, layout, pop.mass());
        interpolator(pop.levelGhostParticlesOld(), pop_momentum_tensor, layout, pop.mass());
    };

    ModelView_t& modelView;
    core::MomentumTensorInterpolator<dimension, interp_order> interpolator{};
};
template<typename ModelView_t>
MomentumTensorComputer(ModelView_t&) -> MomentumTensorComputer<ModelView_t>;


void compute_momentum_tensor(auto& modelView, std::size_t minLvl, std::size_t maxLvl,
                             double timestamp)
{
    MomentumTensorComputer computer{modelView};
    auto& ions             = modelView.getIons();
    auto const n_pops      = ions.size();
    auto const interpolate = [&](auto& layout, auto&&...) {
        std::size_t i = 0;
        for (auto& pop : ions)
            computer.interpolate(pop, i++, layout);
    };
    modelView.visitHierarchy(interpolate, minLvl, maxLvl);
    computer.fillAllSchedules(minLvl, maxLvl, timestamp);
    modelView.visitHierarchy(
        [&](auto&&...) {
            _sum_pop_resources_into_aggregate<decltype(computer)::n_components>(
                modelView, [&](std::size_t idx) -> auto& { return modelView.tmpTensorField(idx); },
                n_pops);
        },
        minLvl, maxLvl);
}


void compute_pop_momentum_tensor(auto& modelView, auto& pop, std::size_t minLvl,
                                 std::size_t maxLvl, double timestamp)
{
    MomentumTensorComputer computer{modelView};
    auto const i           = modelView.getIons().pop_index(pop.name());
    auto const interpolate = [&](auto& layout, auto&&...) { computer.interpolate(pop, i, layout); };
    modelView.visitHierarchy(interpolate, minLvl, maxLvl);
    computer.fillPopSchedules(pop, minLvl, maxLvl, timestamp);
}


// state.V is only kept up to date at t=0 (the MHD solver evolves rho/rhoV/B/Etot, not V/P
// directly) - recompute it from the conservative variables on demand, in place, before a diag
// dump reads it via MHDAccessors.
void compute_mhd_velocity(auto& modelView, std::size_t minLvl, std::size_t maxLvl,
                          double /*timestamp*/)
{
    using GridLayout = typename std::decay_t<decltype(modelView)>::GridLayout;
    auto& state       = modelView.model().state;

    modelView.visitHierarchy(
        [&](GridLayout& layout, std::string const&, std::size_t) {
            core::ToPrimitiveConverter<GridLayout> toPrim{layout};
            toPrim.rhoVToVOnGhostBox(state.rho, state.rhoV, state.V);
        },
        minLvl, maxLvl);
}


// same rationale as compute_mhd_velocity, for the pressure EOS - gamma comes from the
// diagnostic's own fileAttributes (see DiagnosticProperties), not the model, so it is passed in
// rather than read off the state
void compute_mhd_pressure(auto& modelView, double gamma, std::size_t minLvl, std::size_t maxLvl,
                          double /*timestamp*/)
{
    using GridLayout = typename std::decay_t<decltype(modelView)>::GridLayout;
    auto& state       = modelView.model().state;

    modelView.visitHierarchy(
        [&](GridLayout& layout, std::string const&, std::size_t) {
            core::ToPrimitiveConverter<GridLayout> toPrim{layout};
            toPrim.eosEtotToPOnGhostBox(gamma, state.rho, state.rhoV, state.B, state.Etot,
                                        state.P);
        },
        minLvl, maxLvl);
}


} // namespace PHARE::diagnostic

#endif /* PHARE_DIAGNOSTIC_COMPUTERS_DIAGNOSTIC_COMPUTERS_HPP_ */
