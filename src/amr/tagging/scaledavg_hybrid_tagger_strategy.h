#ifndef SCALEDAVG_HYBRID_TAGGER_STRATEGY_H
#define SCALEDAVG_HYBRID_TAGGER_STRATEGY_H

#include "hybrid_tagger_strategy.h"
#include "core/data/grid/gridlayoutdefs.h"
#include "core/data/vecfield/vecfield_component.h"


namespace PHARE::amr
{
template<typename HybridModel>
class ScaledAvgHybridTaggerStrategy : public HybridTaggerStrategy<HybridModel>
{
    using gridlayout_type           = typename HybridModel::gridlayout_type;
    static auto constexpr dimension = HybridModel::dimension;

public:
    void tag(HybridModel& model, gridlayout_type const& layout, int* tags) const override;
};

template<typename HybridModel>
void ScaledAvgHybridTaggerStrategy<HybridModel>::tag(HybridModel& model,
                                                     gridlayout_type const& layout, int* tags) const
{
    auto& Bx = model.state.electromag.B.getComponent(PHARE::core::Component::X);
    auto& By = model.state.electromag.B.getComponent(PHARE::core::Component::Y);
    auto& Bz = model.state.electromag.B.getComponent(PHARE::core::Component::Z);

    auto start
        = layout.physicalStartIndex(PHARE::core::QtyCentering::dual, PHARE::core::Direction::X);
    auto end = layout.physicalEndIndex(PHARE::core::QtyCentering::dual, PHARE::core::Direction::X);

    auto endCell     = layout.nbrCells()[0] - 1;
    double threshold = 0.1;


    if constexpr (dimension == 1)
    {
        for (auto iCell = 0u, ix = start; iCell <= endCell; ++ix, ++iCell)
        {
            auto Bxavg = (Bx(ix - 1) + Bx(ix) + Bx(ix + 1)) / 3.;
            auto Byavg = (By(ix - 1) + By(ix) + By(ix + 1)) / 3.;
            auto Bzavg = (Bz(ix - 1) + Bz(ix) + Bz(ix + 1)) / 3.;

            auto diffx = std::abs(Bxavg - Bx(ix));
            auto diffy = std::abs(Byavg - By(ix));
            auto diffz = std::abs(Bzavg - Bz(ix));

            auto max = std::max({diffx, diffy, diffz});
            if (max > threshold)
            {
                tags[iCell] = 1;
            }
            else
                tags[iCell] = 0;
        }
    }
}
} // namespace PHARE::amr

#endif // SCALEDAVG_HYBRID_TAGGER_STRATEGY_H
