#ifndef PHARE_DIAGNOSTIC_COMPUTERS_DIAGNOSTIC_COMPUTERS_HPP_
#define PHARE_DIAGNOSTIC_COMPUTERS_DIAGNOSTIC_COMPUTERS_HPP_

#include "core/numerics/interpolator/interpolator.hpp"

#include <functional>

namespace PHARE::diagnostic
{

// PUBLIC API
void compute_momentum_tensor(auto& h5Writer);
void compute_pop_momentum_tensor(auto& h5Writer, auto& pop);
void compute_pop_kinetic_energy_flux_vector(auto& h5Writer, auto& pop);



// PRIVATE IMPL


void _fill_all_pop_schedules(auto& h5Writer, auto fn)
{
    auto& modelView           = h5Writer.modelView();
    auto const fill_schedules = [&](auto& lvl) {
        for (std::size_t i = 0; i < modelView.getIons().size(); ++i)
            fn(modelView, lvl, h5Writer.timestamp(), i);
    };
    modelView.onLevels(fill_schedules, h5Writer.minLevel, h5Writer.maxLevel);
}

void _fill_one_pop_schedule(auto& h5Writer, auto fn, auto& pop)
{
    auto& modelView           = h5Writer.modelView();
    auto& ions                = modelView.getIons();
    auto const i              = ions.pop_index(pop.name());
    auto const fill_schedules = [&](auto& lvl) { fn(modelView, lvl, h5Writer.timestamp(), i); };
    modelView.onLevels(fill_schedules, h5Writer.minLevel, h5Writer.maxLevel);
}


template<typename FluidWriter>
struct KineticEnergyFluxVectorComputer
{
    using ModelView_t  = FluidWriter::ModelView_t;
    using Model_t      = ModelView_t::Model_t;
    using level_t      = Model_t::amr_types::level_t;
    using GridLayout_t = Model_t::gridlayout_type;
    auto constexpr static fillFn
        = &ModelView_t::template fillPopKineticEnergyFluxVector<level_t, double, std::size_t>;
    auto constexpr static dimension    = GridLayout_t::dimension;
    auto constexpr static interp_order = GridLayout_t::interp_order;

    void fillPopSchedules(auto& h5Writer, auto& pop)
    {
        _fill_one_pop_schedule(h5Writer, std::mem_fn(fillFn), pop);
    }

    auto interpolate(auto& pop, auto& layout)
    {
        pop.kineticEnergyFlux().zero();
        interpolator(pop, pop.domainParticles(), layout);
        interpolator(pop, pop.levelGhostParticlesOld(), layout);
    };

    FluidWriter& h5Writer;
    core::HeatEnergyFluxVectorInterpolator<dimension, interp_order> interpolator{};
};
template<typename FluidWriter>
KineticEnergyFluxVectorComputer(FluidWriter&) -> KineticEnergyFluxVectorComputer<FluidWriter>;

void compute_pop_kinetic_energy_flux_vector(auto& h5Writer, auto& pop)
{
    KineticEnergyFluxVectorComputer computer{h5Writer};
    auto const interpolate = [&](auto& layout, auto&&...) { computer.interpolate(pop, layout); };
    h5Writer.modelView().visitHierarchy(interpolate, h5Writer.minLevel, h5Writer.maxLevel);
    computer.fillPopSchedules(h5Writer, pop);
}




template<typename FluidWriter>
struct MomentumTensorComputer
{
    // compute the momentum tensor for each population that requires it
    // compute for all ions but that requires the computation of all pop

    // dumps occur after the last substep but before the next first substep
    // at this time, levelGhostPartsNew is emptied and not yet filled
    // and the former levelGhostPartsNew has been moved to levelGhostPartsOld

    using ModelView_t  = FluidWriter::ModelView_t;
    using Model_t      = ModelView_t::Model_t;
    using level_t      = Model_t::amr_types::level_t;
    using GridLayout_t = Model_t::gridlayout_type;
    auto constexpr static fillFn
        = &ModelView_t::template fillPopMomTensor<level_t, double, std::size_t>;
    auto constexpr static dimension    = GridLayout_t::dimension;
    auto constexpr static interp_order = GridLayout_t::interp_order;

    void fillAllSchedules(auto& h5Writer)
    {
        _fill_all_pop_schedules(h5Writer, std::mem_fn(fillFn));
    }

    void fillPopSchedules(auto& h5Writer, auto& pop)
    {
        _fill_one_pop_schedule(h5Writer, std::mem_fn(fillFn), pop);
    }

    auto interpolate(auto& pop, auto& layout)
    {
        auto& pop_momentum_tensor = pop.momentumTensor();
        pop_momentum_tensor.zero();
        interpolator(pop.domainParticles(), pop_momentum_tensor, layout, pop.mass());
        interpolator(pop.levelGhostParticlesOld(), pop_momentum_tensor, layout, pop.mass());
    };

    FluidWriter& h5Writer;
    core::MomentumTensorInterpolator<dimension, interp_order> interpolator{};
};
template<typename FluidWriter>
MomentumTensorComputer(FluidWriter&) -> MomentumTensorComputer<FluidWriter>;


void compute_momentum_tensor(auto& h5Writer)
{
    MomentumTensorComputer computer{h5Writer};
    auto& modelView        = h5Writer.modelView();
    auto const minLvl      = h5Writer.minLevel;
    auto const maxLvl      = h5Writer.maxLevel;
    auto& ions             = modelView.getIons();
    auto const interpolate = [&](auto& layout, auto&&...) {
        for (auto& pop : ions)
            computer.interpolate(pop, layout);
    };
    modelView.visitHierarchy(interpolate, minLvl, maxLvl);
    computer.fillAllSchedules(h5Writer);
    modelView.visitHierarchy([&](auto&&...) { ions.computeFullMomentumTensor(); }, minLvl, maxLvl);
}


void compute_pop_momentum_tensor(auto& h5Writer, auto& pop)
{
    MomentumTensorComputer computer{h5Writer};
    auto const interpolate = [&](auto& layout, auto&&...) { computer.interpolate(pop, layout); };
    h5Writer.modelView().visitHierarchy(interpolate, h5Writer.minLevel, h5Writer.maxLevel);
    computer.fillPopSchedules(h5Writer, pop);
}


} // namespace PHARE::diagnostic

#endif /* PHARE_DIAGNOSTIC_COMPUTERS_DIAGNOSTIC_COMPUTERS_HPP_ */
