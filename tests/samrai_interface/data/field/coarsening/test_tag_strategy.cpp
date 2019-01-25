#include "test_tag_strategy.h"


TagStrategy::TagStrategy(std::map<std::string, int> const &dataToAllocate)
    : dataToAllocate_{dataToAllocate}
{
}

void TagStrategy::initializeLevelData(
    const std::shared_ptr<SAMRAI::hier::PatchHierarchy> &hierarchy, const int levelNumber,
    const double, const bool, const bool, const std::shared_ptr<SAMRAI::hier::PatchLevel> &,
    const bool allocateData)
{
    if (allocateData)
    {
        auto level = hierarchy->getPatchLevel(levelNumber);
        for (auto &patch : *level)
        {
            for (auto const &dataPair : dataToAllocate_)
            {
                auto const &dataId = dataPair.second;
                patch->allocatePatchData(dataId);
            }
        }
    }
}

void TagStrategy::resetHierarchyConfiguration(const std::shared_ptr<SAMRAI::hier::PatchHierarchy> &,
                                              const int, const int)
{
    // do nothing
}
