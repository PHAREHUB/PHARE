#ifndef PHARE_LEVEL_INITIALIZER_HPP
#define PHARE_LEVEL_INITIALIZER_HPP

#include "amr/messengers/messenger.hpp"
#include "amr/physical_models/physical_model.hpp"

namespace PHARE
{
namespace solver
{
    template<typename AMRTypes>
    class LevelInitializer
    {
        using level_ptr       = AMRTypes::level_ptr;
        using hierarchy_ptr   = AMRTypes::hierarchy_ptr;
        using IPhysicalModelT = IPhysicalModel<AMRTypes>;

    public:
        virtual void initialize(hierarchy_ptr const& hierarchy, int levelNumber,
                                level_ptr const& oldLevel, IPhysicalModelT& model,
                                amr::IMessenger<IPhysicalModelT>& messenger, double initDataTime,
                                bool isRegridding)
            = 0;


        virtual ~LevelInitializer() {}
    };
} // namespace solver
} // namespace PHARE
#endif
