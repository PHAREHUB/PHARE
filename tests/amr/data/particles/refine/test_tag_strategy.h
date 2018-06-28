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
#include "utilities/constants.h"
#include "utilities/point/point.h"


#include <map>
#include <string>
#include <type_traits>

using namespace PHARE;

template<std::size_t dimension>
class TagStrategy : public SAMRAI::mesh::StandardTagAndInitStrategy
{
public:
    TagStrategy(std::map<std::string, int> const& dataToAllocate,
                std::shared_ptr<SAMRAI::hier::RefineOperator>& refineOperator,
                bool refineOnlyBorder)
        : dataToAllocate_{dataToAllocate}
        , refineOp_{refineOperator}
        , refineOnlyBorder_{refineOnlyBorder}
    {
        for (auto const& nameToIds : dataToAllocate_)
        {
            algorithm_.registerRefine(nameToIds.second, nameToIds.second, nameToIds.second,
                                      refineOp_);
        }
    }

    virtual ~TagStrategy() = default;

    void initializeLevelData(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                             int const levelNumber, double const, bool const, bool const,
                             std::shared_ptr<SAMRAI::hier::PatchLevel> const& = std::shared_ptr<
                                 SAMRAI::hier::PatchLevel>(),
                             bool const allocateData = true) override
    {
        if (allocateData)
        {
            auto level = hierarchy->getPatchLevel(levelNumber);
            for (auto& patch : *level)
            {
                for (auto const& dataPair : dataToAllocate_)
                {
                    auto const& dataId = dataPair.second;
                    patch->allocatePatchData(dataId);
                }
            }
        }




        if (levelNumber == 0)
        {
            auto level = hierarchy->getPatchLevel(levelNumber);
            for (auto& patch : *level)
            {
                for (auto const& variablesId : dataToAllocate_)
                {
                    auto const& dataId = variablesId.second;
                    auto particlesData = std::dynamic_pointer_cast<ParticlesData<dimension>>(
                        patch->getPatchData(dataId));

                    auto& interior = particlesData->domainParticles;

                    auto particlesGhostBox = particlesData->getGhostBox();
                    particlesGhostBox
                        = AMRToLocal(static_cast<std::add_const_t<decltype(particlesGhostBox)>>(
                                         particlesGhostBox),
                                     particlesGhostBox);



                    // here we are 1D
                    for (int iCellPos = particlesGhostBox.lower(dirX);
                         iCellPos <= particlesGhostBox.upper(dirX); ++iCellPos)
                    {
                        float middle = 0.5;
                        float delta  = 0.25;

                        Particle<dimension> particle;

                        particle.weight = 1.;
                        particle.charge = 1.;
                        particle.v      = {{1.0, 0.0, 0.0}};

                        particle.iCell[dirX] = iCellPos;

                        particle.delta[dirX] = middle - delta;

                        interior.push_back(particle);

                        particle.delta[dirX] = middle + delta;

                        interior.push_back(particle);
                    }
                }
            }
        }
        else
        {
            // create schedule
            if (refineOnlyBorder_)
            {
                auto refineScheduleBorder = algorithm_.createSchedule(
                    std::make_shared<SAMRAI::xfer::PatchLevelBorderFillPattern>(),
                    hierarchy->getPatchLevel(levelNumber), nullptr, levelNumber - 1, hierarchy);

                refineScheduleBorder->fillData(0.);
            }
            else
            {
                auto refineScheduleInterior = algorithm_.createSchedule(
                    std::make_shared<SAMRAI::xfer::PatchLevelInteriorFillPattern>(),
                    hierarchy->getPatchLevel(levelNumber), nullptr, levelNumber - 1, hierarchy);

                refineScheduleInterior->fillData(0.);
            }
        }
    }

    void resetHierarchyConfiguration(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const&,
                                     int const, int const) override
    {
        // do nothing
    }

    static double affineFill(Point<double, 1> position)
    {
        double a = 0.5;
        double b = 2.0;
        return a * position[dirX] + b;
    }



private:
    std::map<std::string, int> dataToAllocate_;
    std::shared_ptr<SAMRAI::hier::RefineOperator> refineOp_;
    SAMRAI::xfer::RefineAlgorithm algorithm_;

    bool const refineOnlyBorder_;
};

#endif
