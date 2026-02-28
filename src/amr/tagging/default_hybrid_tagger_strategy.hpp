#ifndef DEFAULT_HYBRID_TAGGER_STRATEGY_H
#define DEFAULT_HYBRID_TAGGER_STRATEGY_H

#include "hybrid_tagger_strategy.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/data/vecfield/vecfield_component.hpp"
#include "core/data/ndarray/ndarray_vector.hpp"
#include <cstddef>
#include "initializer/data_provider.hpp"

namespace PHARE::amr
{
template<typename HybridModel>
class DefaultHybridTaggerStrategy : public HybridTaggerStrategy<HybridModel>
{
    using gridlayout_type           = typename HybridModel::gridlayout_type;
    static auto constexpr dimension = HybridModel::dimension;


public:
    DefaultHybridTaggerStrategy(initializer::PHAREDict const& dict)
        : threshold_{cppdict::get_value(dict, "threshold", 0.1)}
    {
    }
    void tag(HybridModel& model, gridlayout_type const& layout, int* tags) const override;

private:
    double threshold_ = 0.1;
};

template<typename HybridModel>
void DefaultHybridTaggerStrategy<HybridModel>::tag(HybridModel& model,
                                                   gridlayout_type const& layout, int* tags) const
{
    auto& Bx = model.state.electromag.B.getComponent(PHARE::core::Component::X);
    auto& By = model.state.electromag.B.getComponent(PHARE::core::Component::Y);
    auto& Bz = model.state.electromag.B.getComponent(PHARE::core::Component::Z);

    auto& N = model.state.ions.chargeDensity();

    // we loop on cell indexes for all qties regardless of their centering
    auto const& [start_x, _]
        = layout.physicalStartToEnd(PHARE::core::QtyCentering::dual, PHARE::core::Direction::X);

    // override end_x because the tag buffer does not have ghost cells
    // and physicalEnd will account for ghost cells
    auto const& end_x = layout.nbrCells()[0] - 1;

    // SAMRAI tags int* buffer is FORTRAN ordering so we set false to the view
    bool constexpr c_ordering = false;
    auto tagsv = core::NdArrayView<dimension, int, c_ordering>(tags, layout.nbrCells());

    if constexpr (dimension == 1 and false)
    {
        // at interporder 1 we choose not to tag the last patch cell since
        // the 5 points stencil may go beyond the last ghost node.
        // for interp order 2 and 3 this is ok
        for (auto iCell = 0u, ix = start_x; iCell <= end_x; ++ix, ++iCell)
        {
            auto crit_by_x = (By(ix + 2) - By(ix)) / (1 + By(ix + 1) - By(ix));
            auto crit_bz_x = (Bz(ix + 2) - Bz(ix)) / (1 + Bz(ix + 1) - Bz(ix));
            auto criter    = std::max(crit_by_x, crit_bz_x);

            if (criter > threshold_)
            {
                tagsv(iCell) = 1;
            }
            else
                tagsv(iCell) = 0;
        }
    }
    if constexpr (dimension == 1)
    {
        // at interporder 1 we choose not to tag the last patch cell since
        // the 5 points stencil may go beyond the last ghost node.
        // for interp order 2 and 3 this is ok
        auto constexpr doLastCell = gridlayout_type::nbrGhosts() > 2;
        std::size_t oneOrZero     = doLastCell ? 1 : 0;
        for (auto iCell = 0u, ix = start_x; iCell < end_x + oneOrZero; ++ix, ++iCell)
        {
            auto Byavg     = 0.2 * (By(ix - 2) + By(ix - 1) + By(ix) + By(ix + 1) + By(ix + 2));
            auto Bzavg     = 0.2 * (Bz(ix - 2) + Bz(ix - 1) + Bz(ix) + Bz(ix + 1) + Bz(ix + 2));
            auto Byavgp1   = 0.2 * (By(ix - 1) + By(ix) + By(ix + 1) + By(ix + 2) + By(ix + 3));
            auto Bzavgp1   = 0.2 * (Bz(ix - 1) + Bz(ix) + Bz(ix + 1) + Bz(ix + 2) + Bz(ix + 3));
            auto criter_by = std::abs(Byavgp1 - Byavg) / (1 + std::abs(Byavg));
            auto criter_bz = std::abs(Bzavgp1 - Bzavg) / (1 + std::abs(Bzavg));
            auto criter_b  = std::sqrt(criter_by * criter_by + criter_bz * criter_bz);
            auto criter    = criter_b;

            if (criter > threshold_)
            {
                tagsv(iCell) = 1;
            }
            else
                tagsv(iCell) = 0;
        }
    }
    if constexpr (dimension == 2)
    {
        auto const& [start_y, __]
            = layout.physicalStartToEnd(PHARE::core::QtyCentering::dual, PHARE::core::Direction::Y);

        auto const& end_y = layout.nbrCells()[1] - 1;

        for (auto iTag_x = 0u, ix = start_x; iTag_x <= end_x; ++ix, ++iTag_x)
        {
            for (auto iTag_y = 0u, iy = start_y; iTag_y <= end_y; ++iy, ++iTag_y)
            {
                auto field_diff = [&](auto const& F) //
                {
                    auto const delta_2x = std::abs(F(ix + 2, iy) - F(ix, iy));
                    auto const delta_2y = std::abs(F(ix, iy + 2) - F(ix, iy));
                    auto const delta_x  = std::abs(F(ix + 1, iy) - F(ix, iy));
                    auto const delta_y  = std::abs(F(ix, iy + 1) - F(ix, iy));

                    auto const criter_x = delta_2x / (1 + delta_x);
                    auto const criter_y = delta_2y / (1 + delta_y);

                    return std::make_tuple(criter_x, criter_y);
                };

                auto const& [Bx_x, Bx_y] = field_diff(Bx);
                auto const& [By_x, By_y] = field_diff(By);
                auto const& [Bz_x, Bz_y] = field_diff(Bz);
                auto crit                = std::max({Bx_x, Bx_y, By_x, By_y, Bz_x, Bz_y});

                if (crit > threshold_)
                {
                    tagsv(iTag_x, iTag_y) = 1;
                }
                else
                {
                    tagsv(iTag_x, iTag_y) = 0;
                }
            }
        }
    }
    if constexpr (dimension == 3)
    {
        auto const& [start_y, __]
            = layout.physicalStartToEnd(PHARE::core::QtyCentering::dual, PHARE::core::Direction::Y);
        auto const& [start_z, ___]
            = layout.physicalStartToEnd(PHARE::core::QtyCentering::dual, PHARE::core::Direction::Z);

        auto const& end_y = layout.nbrCells()[1] - 1;
        auto const& end_z = layout.nbrCells()[2] - 1;

        for (auto iTag_x = 0u, ix = start_x; iTag_x <= end_x; ++ix, ++iTag_x)
        {
            for (auto iTag_y = 0u, iy = start_y; iTag_y <= end_y; ++iy, ++iTag_y)
            {
                for (auto iTag_z = 0u, iz = start_z; iTag_z <= end_z; ++iz, ++iTag_z)
                {
                    auto const field_diff = [&](auto const& F) {
                        auto const delta_2x = std::abs(F(ix + 2, iy, iz) - F(ix, iy, iz));
                        auto const delta_2y = std::abs(F(ix, iy + 2, iz) - F(ix, iy, iz));
                        auto const delta_2z = std::abs(F(ix, iy, iz + 2) - F(ix, iy, iz));
                        auto const delta_x  = std::abs(F(ix + 1, iy, iz) - F(ix, iy, iz));
                        auto const delta_y  = std::abs(F(ix, iy + 1, iz) - F(ix, iy, iz));
                        auto const delta_z  = std::abs(F(ix, iy, iz + 1) - F(ix, iy, iz));

                        auto const criter_x = delta_2x / (1 + delta_x);
                        auto const criter_y = delta_2y / (1 + delta_y);
                        auto const criter_z = delta_2z / (1 + delta_z);

                        return std::make_tuple(criter_x, criter_y, criter_z);
                    };

                    auto const& [Bx_x, Bx_y, Bx_z] = field_diff(Bx);
                    auto const& [By_x, By_y, By_z] = field_diff(By);
                    auto const& [Bz_x, Bz_y, Bz_z] = field_diff(Bz);
                    auto const crit
                        = std::max({Bx_x, Bx_y, Bx_z, By_x, By_y, By_z, Bz_x, Bz_y, Bz_z});

                    if (crit > threshold_)
                    {
                        tagsv(iTag_x, iTag_y, iTag_z) = 1;
                    }
                    else
                    {
                        tagsv(iTag_x, iTag_y, iTag_z) = 0;
                    }
                }
            }
        }
    }
}
} // namespace PHARE::amr

#endif // DEFAULT_HYBRID_TAGGER_STRATEGY_H
