#ifndef PHARE_AMR_MHD_LEVEL_INITIALIZER_HPP
#define PHARE_AMR_MHD_LEVEL_INITIALIZER_HPP

#include "amr/level_initializer/level_initializer.hpp"
#include "amr/messengers/messenger.hpp"
#include "amr/messengers/mhd_messenger.hpp"
#include "amr/physical_models/physical_model.hpp"
#include "core/numerics/interpolator/interpolator.hpp"
#include "initializer/data_provider.hpp"

namespace PHARE::solver
{
template<typename MHDModel>
class MHDLevelInitializer : public LevelInitializer<typename MHDModel::amr_types>
{
    using amr_types                    = typename MHDModel::amr_types;
    using hierarchy_t                  = typename amr_types::hierarchy_t;
    using level_t                      = typename amr_types::level_t;
    using patch_t                      = typename amr_types::patch_t;
    using IPhysicalModelT              = IPhysicalModel<amr_types>;
    using IMessengerT                  = amr::IMessenger<IPhysicalModelT>;
    using MHDMessenger                 = amr::MHDMessenger<MHDModel>;
    using GridLayoutT                  = typename MHDModel::gridlayout_type;
    static constexpr auto dimension    = GridLayoutT::dimension;
    static constexpr auto interp_order = GridLayoutT::interp_order;

    inline bool isRootLevel(int levelNumber) const { return levelNumber == 0; }

public:
    MHDLevelInitializer() = default;

    void initialize(std::shared_ptr<hierarchy_t> const& hierarchy, int levelNumber,
                    std::shared_ptr<level_t> const& oldLevel, IPhysicalModelT& model,
                    amr::IMessenger<IPhysicalModelT>& messenger, double initDataTime,
                    bool isRegridding) override
    {
        core::Interpolator<dimension, interp_order> interpolate_;
        auto& mhdModel = static_cast<MHDModel&>(model);
        auto& level    = amr_types::getLevel(*hierarchy, levelNumber);

        if (isRegridding)
        {
            PHARE_LOG_LINE_STR("regriding level " + std::to_string(levelNumber));
            PHARE_LOG_START(3, "mhdLevelInitializer::initialize : regriding block");
            messenger.regrid(hierarchy, levelNumber, oldLevel, model, initDataTime);
            PHARE_LOG_STOP(3, "mhdLevelInitializer::initialize : regriding block");
        }
        else
        {
            if (isRootLevel(levelNumber))
            {
                PHARE_LOG_START(3, "mhdLevelInitializer::initialize : root level init");
                model.initialize(level);
                messenger.fillRootGhosts(model, level, initDataTime);
                PHARE_LOG_STOP(3, "mhdLevelInitializer::initialize : root level init");
            }
            else
            {
                PHARE_LOG_START(3, "mhdLevelInitializer::initialize : initlevel");
                messenger.initLevel(model, level, initDataTime);
                PHARE_LOG_STOP(3, "mhdLevelInitializer::initialize : initlevel");
            }
        }
    }
};

} // namespace PHARE::solver


#endif
