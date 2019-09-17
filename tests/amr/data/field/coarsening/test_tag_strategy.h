#ifndef PHARE_TEST_TAG_STRATEGY_H
#define PHARE_TEST_TAG_STRATEGY_H


#include <SAMRAI/mesh/StandardTagAndInitStrategy.h>


#include <map>
#include <string>


class TagStrategy : public SAMRAI::mesh::StandardTagAndInitStrategy
{
public:
    explicit TagStrategy(std::map<std::string, int> const& dataToAllocate);

    virtual ~TagStrategy() = default;

    void initializeLevelData(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                             int const levelNumber, double const initDataTime,
                             bool const canBeRefined, bool const initialTime,
                             std::shared_ptr<SAMRAI::hier::PatchLevel> const& oldLevel
                             = std::shared_ptr<SAMRAI::hier::PatchLevel>(),
                             bool const allocateData = true) override;

    void resetHierarchyConfiguration(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                                     int const coarsestLevel, int const finestLevel) override;

private:
    std::map<std::string, int> dataToAllocate_;
};

#endif
