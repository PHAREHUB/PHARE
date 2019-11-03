#include "data/particles/particle.h"
#include "data/particles/particles_data.h"
#include "data/particles/refine/split.h"
#include "test_basic_hierarchy.h"
#include "test_tag_strategy.h"


#include "gmock/gmock.h"
#include "gtest/gtest.h"


#include <algorithm>
#include <limits>
#include <vector>


using testing::DoubleNear;
using testing::Eq;

template<std::size_t dimension_, int ratio_, std::size_t interpOrder_, int refineParticlesNbr_,
         ParticlesDataSplitType SplitType_>
struct ParticlesDataSplitTestDescriptors
{
    std::size_t constexpr static dimension            = dimension_;
    int constexpr static ratio                        = ratio_;
    std::size_t constexpr static interpOrder          = interpOrder_;
    int constexpr static refineParticlesNbr           = refineParticlesNbr_;
    static constexpr ParticlesDataSplitType SplitType = SplitType_;
};




template<typename Type>
class aSimpleBasicHierarchyWithTwoLevels : public ::testing::Test
{
public:
    static std::size_t constexpr dimension            = Type::dimension;
    static std::size_t constexpr interpOrder          = Type::interpOrder;
    static int constexpr refineParticlesNbr           = Type::refineParticlesNbr;
    static int constexpr ratio                        = Type::ratio;
    static constexpr ParticlesDataSplitType SplitType = Type::SplitType;


    aSimpleBasicHierarchyWithTwoLevels()
        : basicHierarchy{ratio}
        , hierarchy{basicHierarchy.getHierarchy()}
        , level0{hierarchy.getPatchLevel(0)}
        , level1{hierarchy.getPatchLevel(1)}
        , dataID{basicHierarchy.getVariables().find("proton1")->second}

    {
    }



    BasicHierarchy<dimension, interpOrder, SplitType, refineParticlesNbr> basicHierarchy;
    SAMRAI::hier::PatchHierarchy& hierarchy;
    std::shared_ptr<SAMRAI::hier::PatchLevel> level0;
    std::shared_ptr<SAMRAI::hier::PatchLevel> level1;


    int dataID;


    ~aSimpleBasicHierarchyWithTwoLevels()
    {
        auto db = SAMRAI::hier::VariableDatabase::getDatabase();
        db->removeVariable("proton1");
    }



    std::vector<Particle<dimension>> getRefinedL0Particles()
    {
        std::vector<Particle<dimension>> refinedParticles;

        auto split
            = Split<dimension, interpOrder>(Point<int32, dimension>{ratio}, refineParticlesNbr);

        auto geom        = this->hierarchy.getGridGeometry();
        auto domainBoxes = geom->getPhysicalDomain();

        auto domainBox  = domainBoxes.front();
        auto ghostWidth = static_cast<int>(ghostWidthForParticles<interpOrder>());
        auto lowerCell  = domainBox.lower()[0] - ghostWidth;
        auto upperCell  = domainBox.upper()[0] + ghostWidth;


        for (int iCell = lowerCell; iCell <= upperCell; ++iCell)
        {
            auto coarseParticles = loadCell<dimension>(iCell);

            for (auto const& part : coarseParticles)
            {
                auto coarseOnRefinedGrid{part};

                for (auto iDim = 0u; iDim < dimension; ++iDim)
                {
                    auto delta                      = static_cast<int>(part.delta[iDim] * ratio);
                    coarseOnRefinedGrid.iCell[iDim] = part.iCell[iDim] * ratio + delta;
                    coarseOnRefinedGrid.delta[iDim] = part.delta[iDim] * ratio - delta;
                }

                split(coarseOnRefinedGrid, refinedParticles);
            }
        }

        return refinedParticles;
    }
};



template<typename Type>
class levelOneCoarseBoundaries : public aSimpleBasicHierarchyWithTwoLevels<Type>
{
public:
    using BaseType                         = aSimpleBasicHierarchyWithTwoLevels<Type>;
    static constexpr std::size_t dimension = BaseType::dimension;

    std::vector<Particle<dimension>>
    filterLevelGhostParticles(std::vector<Particle<dimension>> const& refinedParticles,
                              std::shared_ptr<SAMRAI::hier::Patch> const& patch)
    {
        std::vector<Particle<dimension>> levelGhosts;

        auto pData    = patch->getPatchData(this->dataID);
        auto partData = std::dynamic_pointer_cast<ParticlesData<dimension>>(pData);

        auto myBox      = partData->getBox();
        auto myGhostBox = partData->getGhostBox();


        // only works with one patch or clearly distinct patches so that all
        // ghost regions are levelGhosts.
        std::copy_if(std::begin(refinedParticles), std::end(refinedParticles),
                     std::back_inserter(levelGhosts), [&myBox, &myGhostBox](auto const& part) {
                         return isInBox(myGhostBox, part) && !isInBox(myBox, part);
                     });


        return levelGhosts;
    }




    auto& getLevelGhosts(std::shared_ptr<SAMRAI::hier::Patch> const& patch)
    {
        auto pDat    = patch->getPatchData(this->dataID);
        auto partDat = std::dynamic_pointer_cast<ParticlesData<dimension>>(pDat);

        return partDat->levelGhostParticles;
    }
};




TYPED_TEST_SUITE_P(levelOneCoarseBoundaries);


template<typename Type>
class levelOneInterior : public aSimpleBasicHierarchyWithTwoLevels<Type>
{
public:
    using BaseType                         = aSimpleBasicHierarchyWithTwoLevels<Type>;
    static constexpr std::size_t dimension = BaseType::dimension;

    std::vector<Particle<dimension>>
    filterInteriorParticles(std::vector<Particle<dimension>> const& refinedParticles,
                            std::shared_ptr<SAMRAI::hier::Patch> const& patch)
    {
        std::vector<Particle<dimension>> interiors;

        auto pData    = patch->getPatchData(this->dataID);
        auto partData = std::dynamic_pointer_cast<ParticlesData<dimension>>(pData);

        auto myBox = partData->getBox();

        std::copy_if(std::begin(refinedParticles), std::end(refinedParticles),
                     std::back_inserter(interiors),
                     [&myBox](auto const& part) { return isInBox(myBox, part); });

        return interiors;
    }




    auto& getInterior(std::shared_ptr<SAMRAI::hier::Patch> const& patch)
    {
        auto pDat    = patch->getPatchData(this->dataID);
        auto partDat = std::dynamic_pointer_cast<ParticlesData<dimension>>(pDat);

        return partDat->domainParticles;
    }
};




TYPED_TEST_SUITE_P(levelOneInterior);



/*
using ::testing::get;
MATCHER_P2(PositionNearEq, epsilon, dimension, "")
{
    bool nearEq = true;

    for (std::size_t iDir = dirX; iDir < dimension; ++iDir)
    {
        nearEq = nearEq && ((get<0>(arg)[iDir] - get<1>(arg)[iDir]) <= epsilon);
    }

    return nearEq;
}
*/




TYPED_TEST_P(levelOneCoarseBoundaries, areCorrectlyFilledByRefinedSchedule)
{
    auto refinedParticles = this->getRefinedL0Particles();

    for (auto const& patch : *this->level1)
    {
        auto expectedParticles = this->filterLevelGhostParticles(refinedParticles, patch);
        auto& actualParticles  = this->getLevelGhosts(patch);

        ASSERT_EQ(expectedParticles.size(), actualParticles.size());
        ASSERT_GT(actualParticles.size(), 0);

        bool allFound = true;
        for (auto const& particle : expectedParticles)
        {
            auto foundActual = std::find_if(
                std::begin(actualParticles), std::end(actualParticles),

                [&particle](auto const& actualPart) //
                {
                    bool sameCell  = true;
                    bool sameDelta = true;

                    for (auto iDim = 0u; iDim < TypeParam::dimension; ++iDim)
                    {
                        sameCell  = sameCell && actualPart.iCell[iDim] == particle.iCell[iDim];
                        sameDelta = sameDelta
                                    && std::fabs(actualPart.delta[iDim] - particle.delta[iDim])
                                           < std::numeric_limits<float>::epsilon();
                    }

                    return sameCell && sameDelta;
                });

            allFound = allFound && (foundActual != std::end(expectedParticles));
        }

        EXPECT_TRUE(allFound);
    }
}


#if 0

TYPED_TEST_P(levelOneInterior, isCorrectlyFilledByRefinedSchedule)
{
    auto refinedParticles = this->getRefinedL0Particles();
    std::cout << this->interpOrder << " " << refinedParticles.size() << "\n";

    for (auto const& patch : *this->level1)
    {
        auto expectedParticles = this->filterInteriorParticles(refinedParticles, patch);
        auto& actualParticles  = this->getInterior(patch);

        ASSERT_EQ(expectedParticles.size(), actualParticles.size());
        ASSERT_GT(actualParticles.size(), 0);

        bool allFound = true;
        for (auto const& particle : expectedParticles)
        {
            auto foundActual = std::find_if(
                std::begin(actualParticles), std::end(actualParticles),

                [&particle](auto const& actualPart) //
                {
                    bool sameCell  = true;
                    bool sameDelta = true;

                    for (auto iDim = 0u; iDim < TypeParam::dimension; ++iDim)
                    {
                        sameCell  = sameCell && actualPart.iCell[iDim] == particle.iCell[iDim];
                        sameDelta = sameDelta
                                    && std::fabs(actualPart.delta[iDim] - particle.delta[iDim])
                                           < std::numeric_limits<float>::epsilon();
                    }

                    return sameCell && sameDelta;
                });

            allFound = allFound && (foundActual != std::end(expectedParticles));
        }

        EXPECT_TRUE(allFound);
    }
}




REGISTER_TYPED_TEST_SUITE_P(levelOneInterior, isCorrectlyFilledByRefinedSchedule);
#endif
REGISTER_TYPED_TEST_SUITE_P(levelOneCoarseBoundaries, areCorrectlyFilledByRefinedSchedule);



// std::size_t constexpr dimension = 1;
// int constexpr ratio               = 2;
// std::size_t constexpr interpOrder = 1;
// int constexpr refinedParticlesNbr = 2;

// dimension , ratio , interpOrder, refinedParticlesNbr
using ParticlesDataSplitDescriptors1Dr2o1ref2C2F
    = ParticlesDataSplitTestDescriptors<1, 2, 1, 2, ParticlesDataSplitType::coarseBoundary>;
using ParticlesDataSplitDescriptors1Dr2o2ref2C2F
    = ParticlesDataSplitTestDescriptors<1, 2, 2, 2, ParticlesDataSplitType::coarseBoundary>;
using ParticlesDataSplitDescriptors1Dr2o3ref2C2F
    = ParticlesDataSplitTestDescriptors<1, 2, 3, 2, ParticlesDataSplitType::coarseBoundary>;

using ParticlesDataSplitDescriptors1Dr2o1ref3C2F
    = ParticlesDataSplitTestDescriptors<1, 2, 1, 3, ParticlesDataSplitType::coarseBoundary>;
using ParticlesDataSplitDescriptors1Dr2o2ref3C2F
    = ParticlesDataSplitTestDescriptors<1, 2, 2, 3, ParticlesDataSplitType::coarseBoundary>;
using ParticlesDataSplitDescriptors1Dr2o3ref3C2F
    = ParticlesDataSplitTestDescriptors<1, 2, 3, 3, ParticlesDataSplitType::coarseBoundary>;


using ParticlesDataSplitDescriptors1Dr2o1ref2Int
    = ParticlesDataSplitTestDescriptors<1, 2, 1, 2, ParticlesDataSplitType::interior>;
using ParticlesDataSplitDescriptors1Dr2o2ref2Int
    = ParticlesDataSplitTestDescriptors<1, 2, 2, 2, ParticlesDataSplitType::interior>;
using ParticlesDataSplitDescriptors1Dr2o3ref2Int
    = ParticlesDataSplitTestDescriptors<1, 2, 3, 2, ParticlesDataSplitType::interior>;
using ParticlesDataSplitDescriptors1Dr2o1ref3Int
    = ParticlesDataSplitTestDescriptors<1, 2, 1, 3, ParticlesDataSplitType::interior>;
using ParticlesDataSplitDescriptors1Dr2o2ref3Int
    = ParticlesDataSplitTestDescriptors<1, 2, 2, 3, ParticlesDataSplitType::interior>;
using ParticlesDataSplitDescriptors1Dr2o3ref3Int
    = ParticlesDataSplitTestDescriptors<1, 2, 3, 3, ParticlesDataSplitType::interior>;




typedef ::testing::Types<
    ParticlesDataSplitDescriptors1Dr2o1ref2C2F, ParticlesDataSplitDescriptors1Dr2o2ref2C2F,
    ParticlesDataSplitDescriptors1Dr2o3ref2C2F, ParticlesDataSplitDescriptors1Dr2o1ref3C2F,
    ParticlesDataSplitDescriptors1Dr2o2ref3C2F, ParticlesDataSplitDescriptors1Dr2o3ref3C2F>
    ParticlesCoarseToFineDataDescriptorsRange;


typedef ::testing::Types<
    ParticlesDataSplitDescriptors1Dr2o1ref2Int, ParticlesDataSplitDescriptors1Dr2o2ref2Int,
    ParticlesDataSplitDescriptors1Dr2o3ref2Int, ParticlesDataSplitDescriptors1Dr2o1ref3Int,
    ParticlesDataSplitDescriptors1Dr2o2ref3Int, ParticlesDataSplitDescriptors1Dr2o3ref3Int>
    ParticlesInteriorDataDescriptorsRange;


// typedef ::testing::Types<ParticlesDataSplitDescriptors1Dr2o1ref2C2F> TestTest;


INSTANTIATE_TYPED_TEST_SUITE_P(TestCoarseToFine, levelOneCoarseBoundaries,
                               ParticlesCoarseToFineDataDescriptorsRange);

#if 0

INSTANTIATE_TYPED_TEST_SUITE_P(TestInterior, levelOneInterior,
                               ParticlesInteriorDataDescriptorsRange);

#endif

// INSTANTIATE_TYPED_TEST_SUITE_P(TestInterior, levelOneCoarseBoundaries, TestTest);
