#ifndef LOHNER_HYBRID_TAGGER_STRATEGY_H
#define LOHNER_HYBRID_TAGGER_STRATEGY_H

#include "hybrid_tagger_strategy.h"
#include "core/data/grid/gridlayoutdefs.h"
#include "core/data/vecfield/vecfield_component.h"


namespace PHARE::amr
{
template<typename HybridModel>
class LohnerHybridTaggerStrategy : public HybridTaggerStrategy<HybridModel>
{
    using gridlayout_type           = typename HybridModel::gridlayout_type;
    static auto constexpr dimension = HybridModel::dimension;

public:
    void tag(HybridModel& model, gridlayout_type const& layout, int* tags) const override;
};

template<typename HybridModel>
void LohnerHybridTaggerStrategy<HybridModel>::tag(HybridModel& model, gridlayout_type const& layout,
                                                  int* tags) const
{
    auto& Bx = model.state.electromag.B.getComponent(PHARE::core::Component::X);
    auto& By = model.state.electromag.B.getComponent(PHARE::core::Component::Y);
    auto& Bz = model.state.electromag.B.getComponent(PHARE::core::Component::Z);

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

            /*auto derp     = std::abs(Byavgp2 - Byavg);
            auto derm     = std::abs(Byavg - Byavgm2);
            auto der2     = std::abs(Byavgp2 - 2 * Byavg + Byavgm2);
            auto der2p    = std::abs(Byavgp2 + 2 * Byavg + Byavgm2);
            auto lohner_y = der2 / (1 + derp + derm + epsilon * der2p);
            auto lohner   = lohner_y;
            std::cout << lohner << "\n";
*/
            auto lohner = std::abs((Byavgp2 - Byavg) - (Byavgp1 - Byavg)) / (1 + std::abs(Byavg));

            /*
                        auto derp       = std::abs(By(ix + idx) - By(ix));
                        derm            = std::abs(By(ix) - By(ix - idx));
                        der2            = std::abs(By(ix + idx) - 2 * By(ix) + By(ix - idx));
                        der2p           = std::abs(By(ix + idx) + 2 * By(ix) + By(ix - idx));
                        auto Byavg      = 0.2 * (By(ix - 2) + By(ix - 1) + By(ix) + By(ix + 1) +
               By(ix + 2)); denom           = 1 + derp + derm + epsilon * der2p; double lohner_y =
               (Byavg > 0.01) ? der2 / denom : 0; lohner_y        = der2 / denom;

                        derp            = std::abs(Bz(ix + idx) - Bz(ix));
                        derm            = std::abs(Bz(ix) - Bz(ix - idx));
                        der2            = std::abs(Bz(ix + idx) - 2 * Bz(ix) + Bz(ix - idx));
                        der2p           = std::abs(Bz(ix + idx) + 2 * Bz(ix) + Bz(ix - idx));
                        auto Bzavg      = 0.2 * (Bz(ix - 2) + Bz(ix - 1) + Bz(ix) + Bz(ix + 1) +
               Bz(ix + 2)); denom           = 1 + derp + derm + epsilon * der2p; double lohner_z =
               (Bzavg > 0.01) ? der2 / denom : 0; lohner_z        = der2 / denom;
            auto lohner
                = std::sqrt(lohner_x * lohner_x + lohner_y * lohner_y + lohner_z * lohner_z);

            */
            if (lohner > threshold)
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
