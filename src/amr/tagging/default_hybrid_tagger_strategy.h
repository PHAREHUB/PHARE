#ifndef DEFAULT_HYBRID_TAGGER_STRATEGY_H
#define DEFAULT_HYBRID_TAGGER_STRATEGY_H

#include "hybrid_tagger_strategy.h"
#include "core/data/grid/gridlayoutdefs.h"
#include "core/data/vecfield/vecfield_component.h"
#include "core/data/ndarray/ndarray_vector.h"
#include <cstddef>


namespace PHARE::amr
{
template<typename HybridModel>
class DefaultHybridTaggerStrategy : public HybridTaggerStrategy<HybridModel>
{
    using gridlayout_type           = typename HybridModel::gridlayout_type;
    static auto constexpr dimension = HybridModel::dimension;

public:
    void tag(HybridModel& model, gridlayout_type const& layout, int* tags) const override;
};

template<typename HybridModel>
void DefaultHybridTaggerStrategy<HybridModel>::tag(HybridModel& model,
                                                   gridlayout_type const& layout, int* tags) const
{
    auto& Bx = model.state.electromag.B.getComponent(PHARE::core::Component::X);
    auto& By = model.state.electromag.B.getComponent(PHARE::core::Component::Y);
    auto& Bz = model.state.electromag.B.getComponent(PHARE::core::Component::Z);

    auto& N = model.state.ions.density();

    double threshold = 0.1;

    // we loop on cell indexes for all qties regardless of their centering
    auto const& [start_x, start_y]
        = layout.physicalStartToEnd(PHARE::core::QtyCentering::dual, PHARE::core::Direction::X);

    auto endCell_x = layout.nbrCells()[0] - 1;

    // SAMRAI tags int* buffer is FORTRAN ordering so we set false to the view
    auto tagsv
        = core::NdArrayView<dimension, int, int*, /*c_ordering=*/false>(tags, layout.nbrCells());

    if constexpr (dimension == 1)
    {
        // at interporder 1 we choose not to tag the last patch cell since
        // the 5 points stencil may go beyond the last ghost node.
        // for interp order 2 and 3 this is ok
        auto constexpr doLastCell = gridlayout_type::nbrGhosts() > 2;
        std::size_t oneOrZero     = doLastCell ? 1 : 0;
        for (auto iCell = 0u, ix = start_x; iCell < endCell_x + oneOrZero; ++ix, ++iCell)
        {
            auto Byavg     = 0.2 * (By(ix - 2) + By(ix - 1) + By(ix) + By(ix + 1) + By(ix + 2));
            auto Bzavg     = 0.2 * (Bz(ix - 2) + Bz(ix - 1) + Bz(ix) + Bz(ix + 1) + Bz(ix + 2));
            auto Navg      = 0.2 * (N(ix - 2) + N(ix - 1) + N(ix) + N(ix + 1) + N(ix + 2));
            auto Byavgp1   = 0.2 * (By(ix - 1) + By(ix) + By(ix + 1) + By(ix + 2) + By(ix + 3));
            auto Bzavgp1   = 0.2 * (Bz(ix - 1) + Bz(ix) + Bz(ix + 1) + Bz(ix + 2) + Bz(ix + 3));
            auto Navgp1    = 0.2 * (N(ix - 1) + N(ix) + N(ix + 1) + N(ix + 2) + N(ix + 3));
            auto criter_by = std::abs(Byavgp1 - Byavg) / (1 + std::abs(Byavg));
            auto criter_bz = std::abs(Bzavgp1 - Bzavg) / (1 + std::abs(Bzavg));

            auto criter_b = std::sqrt(criter_by * criter_by + criter_bz * criter_bz);
            auto criter_n = std::abs(Navgp1 - Navg) / (1 + std::abs(Navg));

            auto criter = std::max(criter_b, criter_n);

            if (criter > threshold)
            {
                tagsv(iCell) = 1;
            }
            else
                tagsv(iCell) = 0;
        }
    }
    if constexpr (dimension == 2)
    {
        auto const& [start_y, end_y]
            = layout.physicalStartToEnd(PHARE::core::QtyCentering::dual, PHARE::core::Direction::Y);

        auto endCell_y = layout.nbrCells()[1] - 1;

        // at interporder 1 we choose not to tag the last patch cell in both directions since
        // the 5 points stencill may go beyond the last ghost node.
        // for interp order 2 and 3 this is ok
        auto constexpr doLastCell = gridlayout_type::nbrGhosts() > 2;
        std::size_t oneOrZero     = doLastCell ? 1 : 0;
        for (auto iTag_x = 0u, ix = start_x; iTag_x < endCell_x + oneOrZero; ++ix, ++iTag_x)
        {
            for (auto iTag_y = 0u, iy = start_y; iTag_y < endCell_y + oneOrZero; ++iy, ++iTag_y)
            {
                auto field_avg = [&](auto const& F) {
                    return std::make_tuple(0.2
                                               * (F(ix - 2, iy) + F(ix - 1, iy) + F(ix, iy)
                                                  + F(ix + 1, iy) + F(ix + 2, iy)),
                                           0.2
                                               * (F(ix, iy - 2) + F(ix, iy - 1) + F(ix, iy)
                                                  + F(ix, iy + 1) + F(ix, iy + 2)),

                                           0.2
                                               * (F(ix - 1, iy) + F(ix, iy) + F(ix + 1, iy)
                                                  + F(ix + 2, iy) + F(ix + 3, iy)),
                                           0.2
                                               * (F(ix, iy - 1) + F(ix, iy) + F(ix, iy + 1)
                                                  + F(ix, iy + 2) + F(ix, iy + 3)));
                };

                auto const& [Bxavg_x, Bxavg_y, Bxavgp1_x, Bxavgp1_y] = field_avg(Bx);
                auto const& [Byavg_x, Byavg_y, Byavgp1_x, Byavgp1_y] = field_avg(By);
                auto const& [Bzavg_x, Bzavg_y, Bzavgp1_x, Bzavgp1_y] = field_avg(Bz);
                auto const& [Navg_x, Navg_y, Navgp1_x, Navgp1_y]     = field_avg(N);


                auto criter_bx_x = std::abs(Bxavgp1_x - Bxavg_x) / (1 + std::abs(Bxavg_x));
                auto criter_by_x = std::abs(Byavgp1_x - Byavg_x) / (1 + std::abs(Byavg_x));
                auto criter_bz_x = std::abs(Bzavgp1_x - Bzavg_x) / (1 + std::abs(Bzavg_x));

                auto criter_bx_y = std::abs(Bxavgp1_y - Bxavg_y) / (1 + std::abs(Bxavg_y));
                auto criter_by_y = std::abs(Byavgp1_y - Byavg_y) / (1 + std::abs(Byavg_y));
                auto criter_bz_y = std::abs(Bzavgp1_y - Bzavg_y) / (1 + std::abs(Bzavg_y));

                auto criter_bx = std::max(criter_bx_x, criter_bx_y);
                auto criter_by = std::max(criter_by_x, criter_by_y);
                auto criter_bz = std::max(criter_bz_x, criter_bz_y);

                auto criter_b = std::sqrt(criter_bx * criter_bx + criter_by * criter_by
                                          + criter_bz * criter_bz);


                auto criter_n_x = std::abs(Navgp1_x - Navg_x) / (1 + std::abs(Navg_x));
                auto criter_n_y = std::abs(Navgp1_y - Navg_y) / (1 + std::abs(Navg_y));
                auto criter_n   = std::max(criter_n_x, criter_n_y);

                // for now keep 0¨*criter_n to have only magnetic criteria
                auto criter = std::max(criter_b, 0 * criter_n);


                if (criter > threshold)
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
}
} // namespace PHARE::amr

#endif // LOHNER_HYBRID_TAGGER_STRATEGY_H
