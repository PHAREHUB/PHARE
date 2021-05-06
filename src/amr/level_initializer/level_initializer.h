#ifndef PHARE_LEVEL_INITIALIZER_H
#define PHARE_LEVEL_INITIALIZER_H

#include "amr/messengers/messenger.h"
#include "amr/physical_models/physical_model.h"

namespace PHARE
{
namespace solver
{
    template<typename AMRTypes>
    class LevelInitializer
    {
        using level_t         = typename AMRTypes::level_t;
        using patch_t         = typename AMRTypes::patch_t;
        using hierarchy_t     = typename AMRTypes::hierarchy_t;
        using IPhysicalModelT = IPhysicalModel<AMRTypes>;
        using IMessengerT     = amr::IMessenger<IPhysicalModelT>;

    public:
        virtual void initialize(std::shared_ptr<hierarchy_t> const& hierarchy, int levelNumber,
                                std::shared_ptr<level_t> const& oldLevel, IPhysicalModelT& model,
                                amr::IMessenger<IPhysicalModelT>& messenger, double initDataTime,
                                bool isRegridding)
            = 0;


        virtual ~LevelInitializer() {}
    };
} // namespace solver
} // namespace PHARE
#endif
