#ifndef DEFAULT_HYBRID_TAGGER_STRATEGY_H
#define DEFAULT_HYBRID_TAGGER_STRATEGY_H

#include "hybrid_tagger_strategy.h"
#include "core/data/grid/gridlayoutdefs.h"
#include "core/data/vecfield/vecfield_component.h"


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

    // we loop on cell indexes for all qties regardless of their centering
    auto start
        = layout.physicalStartIndex(PHARE::core::QtyCentering::dual, PHARE::core::Direction::X);
    auto end = layout.physicalEndIndex(PHARE::core::QtyCentering::dual, PHARE::core::Direction::X);

    auto endCell     = layout.nbrCells()[0] - 1;
    double threshold = 0.05;
    int idx          = 4;
    double epsilon   = 0.1;


    if constexpr (dimension == 1)
    {
        for (auto iCell = 0u, ix = start; iCell <= endCell; ++ix, ++iCell)
        {
            auto Byavgm2 = 0.2 * (By(ix - 4) + By(ix - 3) + By(ix - 2) + By(ix - 1) + By(ix));
            auto Byavgm1 = 0.2 * (By(ix - 3) + By(ix - 2) + By(ix - 1) + By(ix) + By(ix + 1));
            auto Byavg   = 0.2 * (By(ix - 2) + By(ix - 1) + By(ix) + By(ix + 1) + By(ix + 2));
            auto Byavgp1 = 0.2 * (By(ix - 1) + By(ix) + By(ix + 1) + By(ix + 2) + By(ix + 3));
            auto Byavgp2 = 0.2 * (By(ix) + By(ix + 1) + By(ix + 2) + By(ix + 3) + By(ix + 4));

            auto Bzavgm2 = 0.2 * (Bz(ix - 4) + Bz(ix - 3) + Bz(ix - 2) + Bz(ix - 1) + Bz(ix));
            auto Bzavgm1 = 0.2 * (Bz(ix - 3) + Bz(ix - 2) + Bz(ix - 1) + Bz(ix) + Bz(ix + 1));
            auto Bzavg   = 0.2 * (Bz(ix - 2) + Bz(ix - 1) + Bz(ix) + Bz(ix + 1) + Bz(ix + 2));
            auto Bzavgp1 = 0.2 * (Bz(ix - 1) + Bz(ix) + Bz(ix + 1) + Bz(ix + 2) + Bz(ix + 3));
            auto Bzavgp2 = 0.2 * (Bz(ix) + Bz(ix + 1) + Bz(ix + 2) + Bz(ix + 3) + Bz(ix + 4));

            auto Navgp1 = 0.2 * (N(ix - 1) + N(ix) + N(ix + 1) + N(ix + 2) + N(ix + 3));
            auto Navg   = 0.2 * (N(ix - 2) + N(ix - 1) + N(ix) + N(ix + 1) + N(ix + 2));

            auto criter_by = std::abs(Byavgp1 - Byavg) / (1 + std::abs(Byavg));
            auto criter_bz = std::abs(Bzavgp1 - Bzavg) / (1 + std::abs(Bzavg));
            auto criter_b  = std::sqrt(criter_by * criter_by + criter_bz * criter_bz);
            auto criter_n  = std::abs(Navgp1 - Navg) / (1 + std::abs(Navg));

            auto criter = std::max(criter_b, criter_n);

            if (criter > threshold)
            {
                tags[iCell] = 1;
            }
            else
                tags[iCell] = 0;
        }
    }
}
} // namespace PHARE::amr

#endif // LOHNER_HYBRID_TAGGER_STRATEGY_H
