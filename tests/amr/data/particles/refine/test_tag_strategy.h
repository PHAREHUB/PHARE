#ifndef PHARE_TEST_TAG_STRATEGY_H
#define PHARE_TEST_TAG_STRATEGY_H


#include <SAMRAI/geom/CartesianPatchGeometry.h>
#include <SAMRAI/hier/Patch.h>
#include <SAMRAI/hier/RefineOperator.h>
#include <SAMRAI/mesh/StandardTagAndInitStrategy.h>
#include <SAMRAI/xfer/PatchLevelBorderFillPattern.h>
#include <SAMRAI/xfer/PatchLevelEnhancedFillPattern.h>
#include <SAMRAI/xfer/PatchLevelInteriorFillPattern.h>
#include <SAMRAI/xfer/RefineAlgorithm.h>

#include "core/data/grid/gridlayout.h"
#include "core/data/grid/gridlayoutdefs.h"
#include "amr/data/particles/particles_data.h"
#include "amr/data/particles/refine/particles_data_split.h"
#include "core/utilities/constants.h"
#include "core/utilities/point/point.h"


#include <map>
#include <string>
#include <type_traits>

using namespace PHARE::core;
using namespace PHARE::amr;

template<std::size_t dimension, typename Box>
std::array<int, 3> boxBoundsLower(Box const& box)
{
    return {box.lower(dirX), dimension > 1 ? box.lower(dirY) : 0,
            dimension > 2 ? box.lower(dirZ) : 0};
}
template<std::size_t dimension, typename Box>
std::array<int, 3> boxBoundsUpper(Box const& box)
{
    return {box.upper(dirX), dimension > 1 ? box.upper(dirY) : 0,
            dimension > 2 ? box.upper(dirZ) : 0};
}


template<std::size_t dimension>
std::vector<Particle<dimension>> loadCell(int iCellX, int iCellY, int iCellZ)
{
    std::array<int, 3> _3diCell = {iCellX, iCellY, iCellZ};

    float middle = 0.5;
    float delta  = 0.30f;

    Particle<dimension> particle;
    std::vector<Particle<dimension>> particles;

    particle.weight = 1.;
    particle.charge = 1.;
    particle.v      = {{1.0, 0.0, 0.0}};
    particle.delta.fill(middle);

    for (size_t i = 0; i < dimension; i++)
        particle.iCell[i] = _3diCell[i];

    particle.delta[dirX] = middle - delta;
    particles.push_back(particle);

    particle.delta[dirX] = middle + delta;
    particles.push_back(particle);

    particle.delta[dirX] = middle;
    particles.push_back(particle);

    particle.delta[dirX] = middle - delta / 2;
    particles.push_back(particle);

    particle.delta[dirX] = middle + delta / 2;
    particles.push_back(particle);

    particle.delta[dirX] = middle - delta / 3;
    particles.push_back(particle);

    particle.delta[dirX] = middle + delta / 3;
    particles.push_back(particle);

    return particles;
}



template<std::size_t dimension>
class TagStrategy : public SAMRAI::mesh::StandardTagAndInitStrategy
{
public:
    TagStrategy(std::map<std::string, int> const& dataToAllocate,
                std::shared_ptr<SAMRAI::hier::RefineOperator>& refineOperator,
                ParticlesDataSplitType splitType)
        : dataToAllocate_{dataToAllocate}
        , refineOp_{refineOperator}
        , splitType_{splitType}
    {
        for (auto const& nameToIds : dataToAllocate_)
        {
            algorithm_.registerRefine(nameToIds.second, nameToIds.second, nameToIds.second,
                                      refineOp_);
        }
    }

    virtual ~TagStrategy() = default;

    void initializeLevelData(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                             int const levelNumber, double const initDataTime, bool const,
                             bool const,
                             std::shared_ptr<SAMRAI::hier::PatchLevel> const& = std::shared_ptr<
                                 SAMRAI::hier::PatchLevel>(),
                             bool const allocateData = true) override
    {
        auto level = hierarchy->getPatchLevel(levelNumber);
        if (allocateData)
        {
            for (auto& patch : *level)
            {
                for (auto const& [dataName, dataId] : dataToAllocate_)
                {
                    // auto const& dataId = dataPair.second;
                    patch->allocatePatchData(dataId, initDataTime);
                }
            }
        }
        else
        {
            for (auto& patch : *level)
            {
                for (auto const& [dataName, dataId] : dataToAllocate_)
                {
                    patch->setTime(dataId, initDataTime);
                }
            }
        }




        if (levelNumber == 0)
        {
            for (auto& patch : *level)
            {
                for (auto const& [name, dataId] : dataToAllocate_)
                {
                    auto particlesData = std::dynamic_pointer_cast<ParticlesData<dimension>>(
                        patch->getPatchData(dataId));

                    auto& interior = particlesData->domainParticles;

                    auto particlesBox = particlesData->getBox();

                    auto lower = boxBoundsLower<dimension>(particlesBox);
                    auto upper = boxBoundsUpper<dimension>(particlesBox);

                    for (auto iCellZ = lower[dirZ]; iCellZ <= upper[dirZ]; ++iCellZ)
                    {
                        for (auto iCellY = lower[dirY]; iCellY <= upper[dirY]; ++iCellY)
                        {
                            for (auto iCellX = lower[dirX]; iCellX <= upper[dirX]; ++iCellX)
                            {
                                auto particles = loadCell<dimension>(iCellX, iCellY, iCellZ);
                                interior.insert(std::end(interior), std::begin(particles),
                                                std::end(particles));
                            }
                        }
                    }
                }
            }
        }
        else
        {
            // create schedule
            if (splitType_ == ParticlesDataSplitType::coarseBoundary)
            {
                // warning : the refine operator is of type 'coarseBoundary'
                // therefore it can only work with a border fill pattern.
                // using another fill pattern (e.g. interior) will result in
                // the operator putting refined particles in the wrong particle Array
                // in the destination patch data (ex : interior particles in the coarse to fine
                // particle array)
                auto refineScheduleBorder = algorithm_.createSchedule(
                    std::make_shared<SAMRAI::xfer::PatchLevelBorderFillPattern>(),
                    hierarchy->getPatchLevel(levelNumber), nullptr, levelNumber - 1, hierarchy);

                refineScheduleBorder->fillData(initDataTime);
            }
            else if (splitType_ == ParticlesDataSplitType::interior)
            {
                auto refineScheduleInterior = algorithm_.createSchedule(
                    std::make_shared<SAMRAI::xfer::PatchLevelInteriorFillPattern>(),
                    hierarchy->getPatchLevel(levelNumber), nullptr, levelNumber - 1, hierarchy);

                refineScheduleInterior->fillData(initDataTime);
            }
        }

        // whatever the level is, we need to fill ghosts
        // auto ghostFiller = algorithm_.createSchedule(hierarchy->getPatchLevel(levelNumber));
        // ghostFiller->fillData(initDataTime);
    }



    void resetHierarchyConfiguration(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const&,
                                     int const, int const) override
    {
        // do nothing
    }



private:
    std::map<std::string, int> dataToAllocate_;
    std::shared_ptr<SAMRAI::hier::RefineOperator> refineOp_;
    SAMRAI::xfer::RefineAlgorithm algorithm_;

    ParticlesDataSplitType const splitType_;
};

#endif
