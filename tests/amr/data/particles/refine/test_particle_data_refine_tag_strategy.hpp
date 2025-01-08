#ifndef PHARE_TEST_TAG_STRATEGY_HPP
#define PHARE_TEST_TAG_STRATEGY_HPP

#include "core/def/phare_mpi.hpp"

#include <SAMRAI/geom/CartesianPatchGeometry.h>
#include <SAMRAI/hier/Patch.h>
#include <SAMRAI/hier/RefineOperator.h>
#include <SAMRAI/mesh/StandardTagAndInitStrategy.h>
#include <SAMRAI/xfer/PatchLevelBorderFillPattern.h>
#include <SAMRAI/xfer/PatchLevelEnhancedFillPattern.h>
#include <SAMRAI/xfer/PatchLevelInteriorFillPattern.h>
#include <SAMRAI/xfer/RefineAlgorithm.h>

#include "core/data/grid/gridlayout.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "amr/data/particles/particles_data.hpp"
#include "amr/data/particles/refine/particles_data_split.hpp"
#include "core/utilities/constants.hpp"
#include "core/utilities/point/point.hpp"


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
void loadCell(ParticleArray<dimension>& particles, int iCellX, int iCellY, int iCellZ)
{
    std::array<int, 3> const _3diCell = {iCellX, iCellY, iCellZ};

    double const middle = 0.5;
    double const delta  = 0.30;

    Particle<dimension> particle;

    particle.weight = 1.;
    particle.charge = 1.;
    particle.v      = {{1.0, 0.0, 0.0}};
    particle.delta.fill(middle);

    particle.iCell = sized_array<dimension>(_3diCell);

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
                    auto particlesData
                        = std::dynamic_pointer_cast<ParticlesData<ParticleArray<dimension>>>(
                            patch->getPatchData(dataId));

                    auto& interior = particlesData->domainParticles;

                    auto const particlesBox = particlesData->getBox();

                    auto const lower = boxBoundsLower<dimension>(particlesBox);
                    auto const upper = boxBoundsUpper<dimension>(particlesBox);

                    for (auto iCellX = lower[dirX]; iCellX <= upper[dirX]; ++iCellX)
                        for (auto iCellY = lower[dirY]; iCellY <= upper[dirY]; ++iCellY)
                            for (auto iCellZ = lower[dirZ]; iCellZ <= upper[dirZ]; ++iCellZ)
                                loadCell<dimension>(interior, iCellX, iCellY, iCellZ);
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

    auto domainParticlesForLevel(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                                 int levelNumber)
    {
        std::vector<ParticleArray<dimension>*> particle_arrays;
        auto level = hierarchy->getPatchLevel(levelNumber);
        for (auto& patch : *level)
            for (auto const& [name, dataId] : dataToAllocate_)
                particle_arrays.emplace_back(
                    &std::dynamic_pointer_cast<ParticlesData<ParticleArray<dimension>>>(
                         patch->getPatchData(dataId))
                         ->domainParticles);
        return particle_arrays;
    }

private:
    std::map<std::string, int> dataToAllocate_;
    std::shared_ptr<SAMRAI::hier::RefineOperator> refineOp_;
    SAMRAI::xfer::RefineAlgorithm algorithm_;

    ParticlesDataSplitType const splitType_;
};

#endif
