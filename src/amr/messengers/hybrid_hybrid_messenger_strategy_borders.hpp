#ifndef PHARE_AMR_MESSENGERS_HYBRID_HYBRID_MESSENGER_STRATEGY_BORDERS_HPP
#define PHARE_AMR_MESSENGERS_HYBRID_HYBRID_MESSENGER_STRATEGY_BORDERS_HPP

// HybridHybridMessengerStrategy moment border fills (fillFluxBorders, fillDensityBorders,
// fillIonBorders), kept out of the class so one translation unit does not instantiate the
// whole strategy at once
// include this wherever the strategy is instantiated, or emit them as simulator_messengers does

#include "amr/messengers/hybrid_hybrid_messenger_strategy.hpp"

#include <array>

namespace PHARE::amr
{

template<typename Strategy>
struct HybridMessengerStrategyBorders
{
    // forces the border fills to be emitted in the calling translation unit
    static void emit_moments()
    {
        [[gnu::used]] static constexpr std::array fills{&Strategy::fluxBorders_,
                                                        &Strategy::densityBorders_};
    }

    static void emit_ions() { [[gnu::used]] static constexpr auto fill = &Strategy::ionBorders_; }
};


template<typename HybridModel, typename RefinementParams>
void HybridHybridMessengerStrategy<HybridModel, RefinementParams>::fluxBorders_(
    IonsT& ions, SAMRAI::hier::PatchLevel& level, double const fillTime)
{
    PHARE_LOG_SCOPE(2, "HybridHybridMessengerStrategy::fillFluxBorders");

    auto constexpr N = core::detail::tensor_field_dim_from_rank<1>();

    for (std::size_t i = 0; i < ions.size(); ++i)
    {
        {
            PHARE_LOG_SCOPE(1, "HybridHybridMessengerStrategy::fillFluxBorders::reduce_into");
            patch_exec_(
                level,
                [&](auto& ions_p, auto& scratch_p) {
                    for (std::uint8_t c = 0; c < N; ++c)
                        core::reduce_into(scratch_p[c], ions_p[i].flux()[c]);
                },
                ions, scratch_vecfield());
        }

        {
            PHARE_LOG_SCOPE(1, "HybridHybridMessengerStrategy::fillFluxBorders::schedule");
            popFluxBorderSumRefiners_[i].fill(level.getLevelNumber(), fillTime);
        }

        {
            PHARE_LOG_SCOPE(1, "HybridHybridMessengerStrategy::fillFluxBorders::copy_fields");
            level_exec_(
                level, [&](auto& ions_p, auto&) -> auto& { return ions_p[i].flux(); },
                [&](auto& ions_p, auto& scratch_p, auto const&... tile_idx) {
                    for (std::uint8_t c = 0; c < N; ++c)
                        core::copy_fields(ions_p[i].flux()[c], scratch_p[c], tile_idx...);
                },
                ions, scratch_vecfield());
        }
    }
}


template<typename HybridModel, typename RefinementParams>
void HybridHybridMessengerStrategy<HybridModel, RefinementParams>::densityBorders_(
    IonsT& ions, SAMRAI::hier::PatchLevel& level, double const fillTime)
{
    PHARE_LOG_SCOPE(2, "HybridHybridMessengerStrategy::fillDensityBorders");

    static_assert(core::is_field_v<std::remove_cvref_t<decltype(scratch_field())>>);

    // one particle density per population, the charge density is computed from it
    assert(popDensityBorderSumRefiners_.size() == ions.size());

    auto const border_sum = [&](auto const& density, auto& refiner) {
        PHARE_LOG_SCOPE(1, "HybridHybridMessengerStrategy::fillDensityBorders::border_sum");
        patch_exec_(
            level,
            [&](auto& ions_p, auto& scratch_p) { core::reduce_into(scratch_p, density(ions_p)); },
            ions, scratch_field());

        {
            PHARE_LOG_SCOPE(1, "HybridHybridMessengerStrategy::fillDensityBorders::schedule");
            refiner.fill(level.getLevelNumber(), fillTime);
        }

        level_exec_(
            level, [&](auto& ions_p, auto&) -> auto& { return density(ions_p); },
            [&](auto& ions_p, auto& scratch_p, auto const&... tile_idx) {
                core::copy_fields(density(ions_p), scratch_p, tile_idx...);
            },
            ions, scratch_field());
    };

    for (std::size_t i = 0; i < ions.size(); ++i)
        border_sum([&](auto& ions_p) -> auto& { return ions_p[i].particleDensity(); },
                   popDensityBorderSumRefiners_[i]);
}


template<typename HybridModel, typename RefinementParams>
void HybridHybridMessengerStrategy<HybridModel, RefinementParams>::ionBorders_(
    IonsT& ions, level_t& level, double const fillTime)
{
    PHARE_LOG_SCOPE(2, "HybridHybridMessengerStrategy::fillIonBorders");

    assert(ionFluxBorderMaxRefiners_.size() == 1);    // charge flux
    assert(ionDensityBorderMaxRefiners_.size() == 1); // charge density

    auto constexpr N = core::detail::tensor_field_dim_from_rank<1>();

    auto const border_max = [&](auto& scratch, auto const& moment, auto& refiner) {
        PHARE_LOG_SCOPE(1, "HybridHybridMessengerStrategy::fillIonBorders::border_max");
        patch_exec_(
            level,
            [&](auto& ions_p, auto& scratch_p) { core::copy_fields(scratch_p, moment(ions_p)); },
            ions, scratch);

        {
            PHARE_LOG_SCOPE(1, "HybridHybridMessengerStrategy::fillIonBorders::schedule");
            refiner.fill(level.getLevelNumber(), fillTime);
        }

        level_exec_(
            level, [&](auto& ions_p, auto&) -> auto& { return moment(ions_p); },
            [&](auto& ions_p, auto& scratch_p, auto const&... tile_idx) {
                core::copy_fields(moment(ions_p), scratch_p, tile_idx...);
            },
            ions, scratch);
    };

    border_max(
        scratch_field(), [](auto& ions_p) -> auto& { return ions_p.chargeDensity(); },
        ionDensityBorderMaxRefiners_[0]);

    patch_exec_(
        level,
        [&](auto& ions_p, auto& scratch_p) {
            for (std::uint8_t c = 0; c < N; ++c)
                core::copy_fields(scratch_p[c], ions_p.chargeFlux()[c]);
        },
        ions, scratch_vecfield());
    ionFluxBorderMaxRefiners_[0].fill(level.getLevelNumber(), fillTime);
    level_exec_(
        level, [](auto& ions_p, auto&) -> auto& { return ions_p.chargeFlux(); },
        [&](auto& ions_p, auto& scratch_p, auto const&... tile_idx) {
            for (std::uint8_t c = 0; c < N; ++c)
                core::copy_fields(ions_p.chargeFlux()[c], scratch_p[c], tile_idx...);
        },
        ions, scratch_vecfield());
}

} // namespace PHARE::amr

#endif /*PHARE_AMR_MESSENGERS_HYBRID_HYBRID_MESSENGER_STRATEGY_BORDERS_HPP*/
