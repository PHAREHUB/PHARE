#ifndef PHARE_SIMULATOR_SIMULATOR_MESSENGER_BORDERS_HPP
#define PHARE_SIMULATOR_SIMULATOR_MESSENGER_BORDERS_HPP

// the hybrid messenger moment border fills, split from simulator_messengers.hpp for compiler
// memory, each emit_* function is explicitly instantiated in its own translation unit

#include "phare_solver.hpp"
#include "amr/messengers/hybrid_hybrid_messenger_strategy_borders.hpp"


namespace PHARE
{

template<auto opts>
void emit_hybrid_moment_borders()
{
    if constexpr (has_hybrid_v<opts>)
    {
        using Hybrid   = solver::PHARE_Types<opts>::Hybrid;
        using Model_t  = Hybrid::Model_t;
        using Params_t = Hybrid::RefinementParams_t;
        using Strategy = amr::HybridHybridMessengerStrategy<Model_t, Params_t>;

        amr::HybridMessengerStrategyBorders<Strategy>::emit_moments();
    }
}


template<auto opts>
void emit_hybrid_ion_borders()
{
    if constexpr (has_hybrid_v<opts>)
    {
        using Hybrid   = solver::PHARE_Types<opts>::Hybrid;
        using Model_t  = Hybrid::Model_t;
        using Params_t = Hybrid::RefinementParams_t;
        using Strategy = amr::HybridHybridMessengerStrategy<Model_t, Params_t>;

        amr::HybridMessengerStrategyBorders<Strategy>::emit_ions();
    }
}

} // namespace PHARE

#endif /*PHARE_SIMULATOR_SIMULATOR_MESSENGER_BORDERS_HPP*/
