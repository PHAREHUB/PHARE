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

#include "data/grid/gridlayout.h"
#include "data/grid/gridlayoutdefs.h"
#include "data/particles/particles_data.h"
#include "data/particles/refine/particles_data_split.h"
#include "utilities/constants.h"
#include "utilities/point/point.h"


#include <map>
#include <string>
#include <type_traits>

using namespace PHARE::core;
using namespace PHARE::amr_interface;

template<std::size_t dimension>
std::vector<Particle<dimension>> loadCell(int iCell)
{
    float middle = 0.5;
    float delta  = 0.30f;

    Particle<dimension> tmpParticle;
    std::vector<Particle<dimension>> particles;

    tmpParticle.weight = 1.;
    tmpParticle.charge = 1.;
    tmpParticle.v      = {{1.0, 0.0, 0.0}};

    tmpParticle.iCell[dirX] = iCell;

    tmpParticle.delta[dirX] = middle - delta;
    particles.push_back(tmpParticle);

    tmpParticle.delta[dirX] = middle + delta;
    particles.push_back(tmpParticle);

    tmpParticle.delta[dirX] = middle;
    particles.push_back(tmpParticle);

    tmpParticle.delta[dirX] = middle - delta / 2;
    particles.push_back(tmpParticle);

    tmpParticle.delta[dirX] = middle + delta / 2;
    particles.push_back(tmpParticle);

    tmpParticle.delta[dirX] = middle - delta / 3;
    particles.push_back(tmpParticle);

    tmpParticle.delta[dirX] = middle + delta / 3;
    particles.push_back(tmpParticle);

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

                    auto particlesBox      = particlesData->getBox();
                    auto particlesGhostBox = particlesData->getGhostBox();
                    // particlesBox           = AMRToLocal(
                    //    static_cast<std::add_const_t<decltype(particlesBox)>>(particlesBox),
                    //    particlesGhostBox);

                    // here we are 1D
                    for (int iCellPos = particlesBox.lower(dirX);
                         iCellPos <= particlesBox.upper(dirX); ++iCellPos)
                    {
                        auto particles = loadCell<dimension>(iCellPos);
                        interior.insert(std::end(interior), std::begin(particles),
                                        std::end(particles));
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
