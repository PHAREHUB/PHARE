#ifndef PHARE_DIAGNOSTIC_COMPUTERS_DIAGNOSTIC_COMPUTERS_HPP_
#define PHARE_DIAGNOSTIC_COMPUTERS_DIAGNOSTIC_COMPUTERS_HPP_

#include "core/numerics/interpolator/interpolator.hpp"
#include "core/numerics/primite_conservative_converter/to_primitive_converter.hpp"

#include <cstddef>
#include <string>

namespace PHARE::diagnostic
{

// PUBLIC API - writer-agnostic: only ever needs the model view and the writer's current
// level range/timestamp, none of which are specific to any particular writer type (h5, vtk)
void compute_momentum_tensor(auto& modelView, std::size_t minLvl, std::size_t maxLvl,
                             double timestamp);
void compute_pop_momentum_tensor(auto& modelView, std::size_t popidx, std::size_t minLvl,
                                 std::size_t maxLvl, double timestamp);
void compute_mhd_velocity(auto& modelView, std::size_t minLvl, std::size_t maxLvl,
                          double timestamp);
void compute_mhd_pressure(auto& modelView, double gamma, std::size_t minLvl, std::size_t maxLvl,
                          double timestamp);



// PRIVATE IMPL


template<typename ModelView_t>
struct MomentumTensorComputer
{
    // compute the momentum tensor for each population that requires it
    // compute for all ions but that requires the computation of all pop

    // dumps occur after the last substep but before the next first substep
    // at this time, levelGhostPartsNew is emptied and not yet filled
    // and the former levelGhostPartsNew has been moved to levelGhostPartsOld

    using GridLayout_t                 = ModelView_t::GridLayout;
    auto constexpr static dimension    = GridLayout_t::dimension;
    auto constexpr static interp_order = GridLayout_t::options.interp_order;

    void fillSchedules(std::size_t popidx, std::size_t minLvl, std::size_t maxLvl, double timestamp)
    {
        modelView.onLevels([&](auto& lvl) { modelView.fillPopMomTensor(lvl, timestamp, popidx); },
                           minLvl, maxLvl);
    }

    void interpolate_domain(auto& pop, auto& layout)
    {
        auto& pop_momentum_tensor = pop.momentumTensor();
        pop_momentum_tensor.zero();
        interpolator(pop.domainParticles(), pop_momentum_tensor, layout, pop.mass());
    }

    // after schedules are filled, so ghost contributions are not duplicated
    void interpolate_ghost(auto& pop, auto& layout)
    {
        auto& pop_momentum_tensor = pop.momentumTensor();
        interpolator(pop.levelGhostParticlesOld(), pop_momentum_tensor, layout, pop.mass());
    }

    ModelView_t& modelView;
    core::MomentumTensorInterpolator<dimension, interp_order> interpolator{};
};


template<typename ModelView_t>
MomentumTensorComputer(ModelView_t&) -> MomentumTensorComputer<ModelView_t>;


void compute_momentum_tensor(auto& modelView, std::size_t minLvl, std::size_t maxLvl,
                             double timestamp)
{
    MomentumTensorComputer computer{modelView};
    auto& ions = modelView.getIons();

    modelView.visitHierarchy(
        [&](auto& layout, auto&&...) {
            for (auto& pop : ions)
                computer.interpolate_domain(pop, layout);
        },
        minLvl, maxLvl);

    for (std::size_t i = 0; i < ions.size(); ++i)
        computer.fillSchedules(i, minLvl, maxLvl, timestamp);

    modelView.visitHierarchy(
        [&](auto& layout, auto&&...) {
            for (auto& pop : ions)
                computer.interpolate_ghost(pop, layout);
            ions.computeFullMomentumTensor();
        },
        minLvl, maxLvl);
}


void compute_pop_momentum_tensor(auto& modelView, std::size_t popidx, std::size_t minLvl,
                                 std::size_t maxLvl, double timestamp)
{
    MomentumTensorComputer computer{modelView};
    auto& pop = modelView.getIons()[popidx];

    modelView.visitHierarchy(
        [&](auto& layout, auto&&...) { computer.interpolate_domain(pop, layout); }, minLvl, maxLvl);
    computer.fillSchedules(popidx, minLvl, maxLvl, timestamp);
    modelView.visitHierarchy(
        [&](auto& layout, auto&&...) { computer.interpolate_ghost(pop, layout); }, minLvl, maxLvl);
}


// V is only kept up to date at t=0 (the MHD solver evolves rho/rhoV/B/Etot, not V/P
// directly) - recompute it from the conservative variables on demand before a diag dump
void compute_mhd_velocity(auto& modelView, std::size_t minLvl, std::size_t maxLvl,
                          double /*timestamp*/)
{
    using GridLayout = std::decay_t<decltype(modelView)>::GridLayout;

    modelView.visitHierarchy(
        [&](GridLayout& layout, std::string const&, std::size_t) {
            core::ToPrimitiveConverter<GridLayout> toPrim{layout};
            toPrim.rhoVToVOnGhostBox(modelView.getRho(), modelView.getRhoV(), modelView.getV());
        },
        minLvl, maxLvl);
}


// same rationale as compute_mhd_velocity, for the pressure EOS - gamma comes from the
// diagnostic's own fileAttributes (see DiagnosticProperties), not the model, so it is passed in
void compute_mhd_pressure(auto& modelView, double gamma, std::size_t minLvl, std::size_t maxLvl,
                          double /*timestamp*/)
{
    using GridLayout = std::decay_t<decltype(modelView)>::GridLayout;

    modelView.visitHierarchy(
        [&](GridLayout& layout, std::string const&, std::size_t) {
            core::ToPrimitiveConverter<GridLayout> toPrim{layout};
            toPrim.eosEtotToPOnGhostBox(gamma, modelView.getRho(), modelView.getRhoV(),
                                        modelView.getB(), modelView.getEtot(), modelView.getP());
        },
        minLvl, maxLvl);
}


} // namespace PHARE::diagnostic

#endif /* PHARE_DIAGNOSTIC_COMPUTERS_DIAGNOSTIC_COMPUTERS_HPP_ */
