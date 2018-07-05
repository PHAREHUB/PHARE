#include "data/particles/refine/split.h"
#include "test_basic_hierarchy.h"
#include "test_tag_strategy.h"


#include "gmock/gmock.h"
#include "gtest/gtest.h"

using testing::DoubleNear;
using testing::Eq;

template<std::size_t dimension_, int ratio_, std::size_t interpOrder_, int refineParticlesNbr_>
struct ParticlesDataSplitTestDescriptors
{
    std::size_t constexpr static dimension   = dimension_;
    int constexpr static ratio               = ratio_;
    std::size_t constexpr static interpOrder = interpOrder_;
    int constexpr static refineParticlesNbr  = refineParticlesNbr_;
};

template<typename Type>
class aParticlesDataSplitOperator : public ::testing::Test
{
};

TYPED_TEST_CASE_P(aParticlesDataSplitOperator);




TYPED_TEST_P(aParticlesDataSplitOperator, splitIsExpectedForCoarseBoundary)
{
    //
    std::size_t constexpr dimension   = TypeParam::dimension;
    int constexpr ratio               = TypeParam::ratio;
    std::size_t constexpr interpOrder = TypeParam::interpOrder;
    int constexpr refinedParticlesNbr = TypeParam::refineParticlesNbr;


    BasicHierarchy<dimension, interpOrder, ParticlesDataSplitType::coarseBoundary,
                   refinedParticlesNbr>
        basicHierarchy{ratio};

    Split<dimension, interpOrder> split{{{ratio}}, refinedParticlesNbr};

    auto& hierarchy = basicHierarchy.getHierarchy();


    std::vector<Particle<dimension>> fineParticles;

    std::vector<Point<double, dimension>> fineParticlesPositionFromFill;

    std::vector<Box<double, dimension>> boxes;
    std::vector<Box<double, dimension>> boxesCandidateForSplit;

    auto gridGeom = std::dynamic_pointer_cast<SAMRAI::geom::CartesianGridGeometry const>(
        hierarchy.getGridGeometry());

    auto* coarseDx = gridGeom->getDx();


    // the physical domain start at 0. and finish at 1.0
    // the first refine box is 4 , 15
    // the second refine box is 30 , 50

    // for interpOrder 1 we have one ghost cell to fill
    // for interpOrder 2 and 3 we have 2 ghost cell to fill

    // the border are at : 8 , 32 , 60 , 102

    if constexpr (dimension == 1)
    {
        //

        double maxDistanceX = ((interpOrder) / 2.) * coarseDx[dirX];


        for (int fineIndex : {8, 32, 60, 102})
        {
            Box<double, dimension> box;

            static_assert(interpOrder > 0 && interpOrder < 4,
                          "Error out of range for interpOrder test");

            // ie left border
            if (fineIndex == 8 || fineIndex == 60)
            {
                if constexpr (interpOrder == 1)
                {
                    box.lower[dirX] = (fineIndex - 1) * coarseDx[dirX] / ratio;
                    box.upper[dirX] = fineIndex * coarseDx[dirX] / ratio;
                }
                else if constexpr (interpOrder == 2 || interpOrder == 3)
                {
                    //
                    box.lower[dirX] = (fineIndex - 2) * coarseDx[dirX] / ratio;
                    box.upper[dirX] = fineIndex * coarseDx[dirX] / ratio;
                }
            }
            else // right border
            {
                if constexpr (interpOrder == 1)
                {
                    box.lower[dirX] = fineIndex * coarseDx[dirX] / ratio;
                    box.upper[dirX] = (fineIndex + 1) * coarseDx[dirX] / ratio;
                }
                else if constexpr (interpOrder == 2 || interpOrder == 3)
                {
                    //
                    box.lower[dirX] = fineIndex * coarseDx[dirX] / ratio;
                    box.upper[dirX] = (fineIndex + 2) * coarseDx[dirX] / ratio;
                }
            }

            boxes.push_back(box);

            box.lower[dirX] = box.lower[dirX] - maxDistanceX;
            box.upper[dirX] = box.upper[dirX] + maxDistanceX;

            boxesCandidateForSplit.push_back(box);
        }
    }

    auto level0 = hierarchy.getPatchLevel(0);
    for (auto const& patch : *level0)
    {
        auto specie1It = basicHierarchy.getVariables().find("specie1");
        auto dataId    = specie1It->second;

        auto particlesData
            = std::dynamic_pointer_cast<ParticlesData<dimension>>(patch->getPatchData(dataId));

        auto patchGeom = std::dynamic_pointer_cast<SAMRAI::geom::CartesianPatchGeometry>(
            patch->getPatchGeometry());


        auto* xLower = patchGeom->getXLower();
        auto* dx     = patchGeom->getDx();

        Point<double, dimension> origin;

        for (auto iDir = dirX; iDir < dimension; ++iDir)
        {
            origin[iDir] = xLower[iDir] - particlesData->getGhostCellWidth()[iDir] * dx[iDir];
        }


        auto candidateToSplit = [dx, &origin, &boxesCandidateForSplit](auto const& particle) {
            Point<double, dimension> particlePosition{positionAsPoint(particle, dx, origin)};

            return isIn(particlePosition, boxesCandidateForSplit);
        };

        for (auto const& particle : particlesData->domainParticles)
        {
            //
            if (candidateToSplit(particle))
            {
                Particle<dimension> fineParticle{particle};
                Point<double, dimension> normalizedPosition;

                Point<double, dimension> particlePosition{positionAsPoint(particle, dx, origin)};

                auto ghostCellWidth = particlesData->getGhostCellWidth();

                for (auto iDir = dirX; iDir < dimension; ++iDir)
                {
                    double fineDx = dx[iDir] / ratio;

                    double normalizedPosition = particlePosition[iDir] / fineDx;
                    fineParticle.iCell[iDir]  = static_cast<int>(normalizedPosition);
                    fineParticle.delta[iDir]
                        = static_cast<float>(normalizedPosition - fineParticle.iCell[iDir]);
                }

                std::vector<Particle<dimension>> splittedParticles;


                auto isInSplittedRegion = [&boxes, dx, dimension](auto const& particle) {
                    Point<double, dimension> fineDx;
                    Point<double, dimension> origin;
                    for (auto iDir = dirX; iDir < dimension; ++iDir)
                    {
                        fineDx[iDir] = dx[iDir] / ratio;
                        origin[iDir] = 0.;
                    }

                    Point<double, dimension> particlePosition{
                        positionAsPoint(particle, fineDx, origin)};

                    return isIn(particlePosition, boxes);
                };

                split(fineParticle, splittedParticles);

                std::copy_if(std::begin(splittedParticles), std::end(splittedParticles),
                             std::back_inserter(fineParticles), isInSplittedRegion);
            }
        }
    }

    int const expectedNbrParticles{static_cast<int>(fineParticles.size())};

    int countParticlesInCoarseBoundary{0};

    SAMRAI::hier::IntVector ghostCellWidth{SAMRAI::tbox::Dimension{dimension}};


    auto level1 = hierarchy.getPatchLevel(1);
    for (auto const& patch : *level1)
    {
        auto specie1It = basicHierarchy.getVariables().find("specie1");
        auto dataId    = specie1It->second;

        auto particlesData
            = std::dynamic_pointer_cast<ParticlesData<dimension>>(patch->getPatchData(dataId));

        countParticlesInCoarseBoundary += particlesData->coarseToFineParticles.size();

        auto patchGeom = std::dynamic_pointer_cast<SAMRAI::geom::CartesianPatchGeometry>(
            patch->getPatchGeometry());

        ghostCellWidth = particlesData->getGhostCellWidth();

        auto* xLower = patchGeom->getXLower();
        auto* dx     = patchGeom->getDx();

        auto computePosition = [xLower, dx, &ghostCellWidth, dimension](auto const& particle) {
            Point<double, dimension> origin;

            for (auto iDir = dirX; iDir < dimension; ++iDir)
            {
                origin[iDir] = xLower[iDir] - ghostCellWidth[iDir] * dx[iDir];
            }

            return positionAsPoint(particle, dx, origin);
        };
        std::transform(std::begin(particlesData->coarseToFineParticles),
                       std::end(particlesData->coarseToFineParticles),
                       std::back_inserter(fineParticlesPositionFromFill), computePosition);
    }

    // compute position for the expected particles
    // here the particle is considered on a single grid covering all the domain
    // That come from the fact that we do not use the fine level origin for the expected
    // particles (since we manually split them on the coarse level)
    auto computePosition = [coarseDx, dimension](auto const& particle) {
        Point<double, dimension> fineDx;
        Point<double, dimension> origin;
        for (auto iDir = dirX; iDir < dimension; ++iDir)
        {
            fineDx[iDir] = coarseDx[iDir] / ratio;
            origin[iDir] = 0.;
        }


        return positionAsPoint(particle, fineDx, origin);
    };

    std::vector<Point<double, dimension>> expectedPosition;

    std::transform(std::begin(fineParticles), std::end(fineParticles),
                   std::back_inserter(expectedPosition), computePosition);

    auto comparePosition = [](auto const& lhs, auto const& rhs) { return lhs[dirX] < rhs[dirX]; };

    std::sort(std::begin(expectedPosition), std::end(expectedPosition), comparePosition);
    std::sort(std::begin(fineParticlesPositionFromFill), std::end(fineParticlesPositionFromFill),
              comparePosition);


    EXPECT_EQ(expectedNbrParticles, countParticlesInCoarseBoundary);


    for (std::size_t i = 0u; i < fineParticlesPositionFromFill.size(); ++i)
    {
        for (auto iDir = dirX; iDir < dimension; ++iDir)
        {
            double const epsilon     = 1e-7;
            bool foundAMatch         = false;
            std::size_t indexToMatch = 0;
            for (std::size_t j = 0u; j < expectedPosition.size(); ++j)
            {
                double expectedToFind = fineParticlesPositionFromFill[i][iDir];
                if (expectedToFind >= expectedPosition[j][iDir] - epsilon
                    && expectedToFind <= expectedPosition[j][iDir] + epsilon)
                {
                    foundAMatch  = true;
                    indexToMatch = j;
                    break;
                }
            }
            EXPECT_TRUE(foundAMatch);
            if (foundAMatch)
            {
                EXPECT_NEAR(expectedPosition[indexToMatch][iDir],
                            fineParticlesPositionFromFill[i][iDir], epsilon);
            }
        }
    }

    basicHierarchy.TearDown();
}




TYPED_TEST_P(aParticlesDataSplitOperator, splitIsExpectedForInterior)
{
    //
    std::size_t constexpr dimension   = TypeParam::dimension;
    int constexpr ratio               = TypeParam::ratio;
    std::size_t constexpr interpOrder = TypeParam::interpOrder;
    int constexpr refinedParticlesNbr = TypeParam::refineParticlesNbr;


    BasicHierarchy<dimension, interpOrder, ParticlesDataSplitType::interior, refinedParticlesNbr>
        basicHierarchy{ratio};

    Split<dimension, interpOrder> split{{{ratio}}, refinedParticlesNbr};

    auto& hierarchy = basicHierarchy.getHierarchy();


    std::vector<Particle<dimension>> fineParticles;

    std::vector<Point<double, dimension>> fineParticlesPositionFromFill;

    std::vector<Box<double, dimension>> boxes;
    std::vector<Box<double, dimension>> boxesCandidateForSplit;

    auto gridGeom = std::dynamic_pointer_cast<SAMRAI::geom::CartesianGridGeometry const>(
        hierarchy.getGridGeometry());

    auto* coarseDx = gridGeom->getDx();


    // the physical domain start at 0. and finish at 1.0
    // the first refine box is 4 , 15
    // the second refine box is 30 , 50

    if constexpr (dimension == 1)
    {
        //

        double maxDistanceX = ((interpOrder) / 2.) * coarseDx[dirX];


        for (auto const& fineIndex : std::array<std::array<int, 2>, 2>{{{8, 30}, {60, 100}}})
        {
            Box<double, dimension> box;

            static_assert(interpOrder > 0 && interpOrder < 4,
                          "Error out of range for interpOrder test");


            box.lower[dirX] = fineIndex[0] * coarseDx[dirX] / ratio;
            box.upper[dirX] = (fineIndex[1] + 1) * coarseDx[dirX] / ratio;


            boxes.push_back(box);

            box.lower[dirX] = box.lower[dirX] - maxDistanceX;
            box.upper[dirX] = box.upper[dirX] + maxDistanceX;

            boxesCandidateForSplit.push_back(box);
        }
    }

    auto level0 = hierarchy.getPatchLevel(0);
    for (auto const& patch : *level0)
    {
        auto specie1It = basicHierarchy.getVariables().find("specie1");
        auto dataId    = specie1It->second;

        auto particlesData
            = std::dynamic_pointer_cast<ParticlesData<dimension>>(patch->getPatchData(dataId));

        auto patchGeom = std::dynamic_pointer_cast<SAMRAI::geom::CartesianPatchGeometry>(
            patch->getPatchGeometry());


        auto* xLower = patchGeom->getXLower();
        auto* dx     = patchGeom->getDx();

        Point<double, dimension> origin;

        for (auto iDir = dirX; iDir < dimension; ++iDir)
        {
            origin[iDir] = xLower[iDir] - particlesData->getGhostCellWidth()[iDir] * dx[iDir];
        }


        auto candidateToSplit = [dx, &origin, &boxesCandidateForSplit](auto const& particle) {
            Point<double, dimension> particlePosition{positionAsPoint(particle, dx, origin)};

            return isIn(particlePosition, boxesCandidateForSplit);
        };

        for (auto const& particle : particlesData->domainParticles)
        {
            //
            if (candidateToSplit(particle))
            {
                Particle<dimension> fineParticle{particle};
                Point<double, dimension> normalizedPosition;

                Point<double, dimension> particlePosition{positionAsPoint(particle, dx, origin)};

                auto ghostCellWidth = particlesData->getGhostCellWidth();

                for (auto iDir = dirX; iDir < dimension; ++iDir)
                {
                    double fineDx = dx[iDir] / ratio;

                    double normalizedPosition = particlePosition[iDir] / fineDx;
                    fineParticle.iCell[iDir]  = static_cast<int>(normalizedPosition);
                    fineParticle.delta[iDir]
                        = static_cast<float>(normalizedPosition - fineParticle.iCell[iDir]);
                }

                std::vector<Particle<dimension>> splittedParticles;


                auto isInSplittedRegion = [&boxes, dx, dimension](auto const& particle) {
                    Point<double, dimension> fineDx;
                    Point<double, dimension> origin;
                    for (auto iDir = dirX; iDir < dimension; ++iDir)
                    {
                        fineDx[iDir] = dx[iDir] / ratio;
                        origin[iDir] = 0.;
                    }

                    Point<double, dimension> particlePosition{
                        positionAsPoint(particle, fineDx, origin)};

                    return isIn(particlePosition, boxes);
                };

                split(fineParticle, splittedParticles);

                std::copy_if(std::begin(splittedParticles), std::end(splittedParticles),
                             std::back_inserter(fineParticles), isInSplittedRegion);
            }
        }
    }

    int const expectedNbrParticles{static_cast<int>(fineParticles.size())};

    int countParticlesInInterior{0};

    SAMRAI::hier::IntVector ghostCellWidth{SAMRAI::tbox::Dimension{dimension}};


    auto level1 = hierarchy.getPatchLevel(1);
    for (auto const& patch : *level1)
    {
        auto specie1It = basicHierarchy.getVariables().find("specie1");
        auto dataId    = specie1It->second;

        auto particlesData
            = std::dynamic_pointer_cast<ParticlesData<dimension>>(patch->getPatchData(dataId));

        countParticlesInInterior += particlesData->domainParticles.size();

        auto patchGeom = std::dynamic_pointer_cast<SAMRAI::geom::CartesianPatchGeometry>(
            patch->getPatchGeometry());

        ghostCellWidth = particlesData->getGhostCellWidth();

        auto* xLower = patchGeom->getXLower();
        auto* dx     = patchGeom->getDx();

        auto computePosition = [xLower, dx, &ghostCellWidth, dimension](auto const& particle) {
            Point<double, dimension> origin;

            for (auto iDir = dirX; iDir < dimension; ++iDir)
            {
                origin[iDir] = xLower[iDir] - ghostCellWidth[iDir] * dx[iDir];
            }

            return positionAsPoint(particle, dx, origin);
        };
        std::transform(std::begin(particlesData->coarseToFineParticles),
                       std::end(particlesData->coarseToFineParticles),
                       std::back_inserter(fineParticlesPositionFromFill), computePosition);
    }

    // compute position for the expected particles
    // here the particle is considered on a single grid covering all the domain
    // That come from the fact that we do not use the fine level origin for the expected
    // particles (since we manually split them on the coarse level)
    auto computePosition = [coarseDx, dimension](auto const& particle) {
        Point<double, dimension> fineDx;
        Point<double, dimension> origin;
        for (auto iDir = dirX; iDir < dimension; ++iDir)
        {
            fineDx[iDir] = coarseDx[iDir] / ratio;
            origin[iDir] = 0.;
        }


        return positionAsPoint(particle, fineDx, origin);
    };

    std::vector<Point<double, dimension>> expectedPosition;

    std::transform(std::begin(fineParticles), std::end(fineParticles),
                   std::back_inserter(expectedPosition), computePosition);

    auto comparePosition = [](auto const& lhs, auto const& rhs) { return lhs[dirX] < rhs[dirX]; };

    std::sort(std::begin(expectedPosition), std::end(expectedPosition), comparePosition);
    std::sort(std::begin(fineParticlesPositionFromFill), std::end(fineParticlesPositionFromFill),
              comparePosition);


    EXPECT_EQ(expectedNbrParticles, countParticlesInInterior);


    for (std::size_t i = 0u; i < fineParticlesPositionFromFill.size(); ++i)
    {
        for (auto iDir = dirX; iDir < dimension; ++iDir)
        {
            double const epsilon     = 1e-7;
            bool foundAMatch         = false;
            std::size_t indexToMatch = 0;
            for (std::size_t j = 0u; j < expectedPosition.size(); ++j)
            {
                double expectedToFind = fineParticlesPositionFromFill[i][iDir];
                if (expectedToFind >= expectedPosition[j][iDir] - epsilon
                    && expectedToFind <= expectedPosition[j][iDir] + epsilon)
                {
                    foundAMatch  = true;
                    indexToMatch = j;
                    break;
                }
            }
            EXPECT_TRUE(foundAMatch);
            if (foundAMatch)
            {
                EXPECT_NEAR(expectedPosition[indexToMatch][iDir],
                            fineParticlesPositionFromFill[i][iDir], epsilon);
            }
        }
    }

    basicHierarchy.TearDown();
}


REGISTER_TYPED_TEST_CASE_P(aParticlesDataSplitOperator, splitIsExpectedForCoarseBoundary,
                           splitIsExpectedForInterior);



// std::size_t constexpr dimension = 1;
// int constexpr ratio               = 2;
// std::size_t constexpr interpOrder = 1;
// int constexpr refinedParticlesNbr = 2;

// dimension , ratio , interpOrder, refinedParticlesNbr
using ParticlesDataSplitDescriptors1Dr2o1ref2 = ParticlesDataSplitTestDescriptors<1, 2, 1, 2>;
using ParticlesDataSplitDescriptors1Dr2o2ref2 = ParticlesDataSplitTestDescriptors<1, 2, 2, 2>;
using ParticlesDataSplitDescriptors1Dr2o3ref2 = ParticlesDataSplitTestDescriptors<1, 2, 3, 2>;



typedef ::testing::Types<ParticlesDataSplitDescriptors1Dr2o1ref2,
                         ParticlesDataSplitDescriptors1Dr2o2ref2,
                         ParticlesDataSplitDescriptors1Dr2o3ref2>
    ParticlesDataDescriptorsRange;

INSTANTIATE_TYPED_TEST_CASE_P(TestUsingParticlesPositionThat, aParticlesDataSplitOperator,
                              ParticlesDataDescriptorsRange);
