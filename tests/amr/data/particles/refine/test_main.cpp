#include "data/particles/particles_data.h"
#include "data/particles/particles_variable.h"
#include "data/particles/refine/particles_data_split_on_coarse_boundary.h"
#include "test_basic_hierarchy.h"
#include "tools/amr_utils.h"

#include <SAMRAI/tbox/SAMRAIManager.h>
#include <SAMRAI/tbox/SAMRAI_MPI.h>



#include "gmock/gmock.h"
#include "gtest/gtest.h"


using testing::Eq;

using namespace PHARE;




TEST(AParticlesDataSplitDataOnCoarseBoundaryOperator,
     canBeUsedWithAPatchLevelFillPatternToFillOnlyCoarseBoundary)
{
    //
    std::size_t constexpr dimension{1};
    std::size_t constexpr interpOrder{1};
    int const ratio = 2;
    BasicHierarchy<dimension, interpOrder> basicHierarchy{ratio};


    auto& hierarchy = basicHierarchy.getHierarchy();

    auto level = hierarchy.getPatchLevel(1);

    for (auto& patch : *level)
    {
        for (auto const& speciesId : basicHierarchy.getVariables())
        {
            auto const& dataId = speciesId.second;

            auto particlesData
                = std::dynamic_pointer_cast<ParticlesData<dimension>>(patch->getPatchData(dataId));


            EXPECT_EQ(0, particlesData->domainParticles.size());
            EXPECT_EQ(0, particlesData->ghostParticles.size());

            auto& patchBox = patch->getBox();


            int const lowerFine1 = ratio * 4;
            int const upperFine1 = ratio * 15;

            int const lowerFine2 = ratio * 30;
            int const upperFine2 = ratio * 50;

            int expectedNumberOfParticles{0};

            int const particlePerCell{2};

            auto const& ghostCellWidth = particlesData->getGhostCellWidth();

            auto const& ghostBox = particlesData->getGhostBox();

            constexpr std::size_t nbrBound{2};

            // TODO : hard coded for interpOrder 1

            expectedNumberOfParticles = nbrBound * (particlePerCell / 2 + particlePerCell / 2);


            EXPECT_EQ(expectedNumberOfParticles, particlesData->coarseToFineParticles.size());

            for (auto const& particle : particlesData->coarseToFineParticles)
            {
                //
                EXPECT_GE(particle.iCell[dirX] + ghostBox.lower(dirX),
                          lowerFine1 - ratio * ghostCellWidth(dirX));
                if (patchBox.lower(dirX) >= upperFine1)
                {
                    EXPECT_GE(particle.iCell[dirX] + ghostBox.lower(dirX),
                              lowerFine2 - ratio * ghostCellWidth(dirX));
                    EXPECT_LE(particle.iCell[dirX] + ghostBox.lower(dirX),
                              upperFine2 + ratio * ghostCellWidth(dirX));
                }

                if (patchBox.upper(dirX) <= upperFine1)
                {
                    EXPECT_LE(particle.iCell[dirX] + ghostBox.lower(dirX),
                              upperFine1 + ratio * ghostCellWidth(dirX));
                }
            }
        }
    }

    basicHierarchy.TearDown();
}



TEST(AParticlesDataSplitDataOperator, canBeUsedWithAPatchLevelFillPatternToFillOnlyInterior)
{
    //
    std::size_t constexpr dimension{1};
    std::size_t constexpr interpOrder{1};
    int const ratio = 2;


    bool refineOnlyInterior{true};
    BasicHierarchy<dimension, interpOrder> basicHierarchy{ratio, !refineOnlyInterior};


    auto& hierarchy = basicHierarchy.getHierarchy();

    auto level = hierarchy.getPatchLevel(1);




    for (auto& patch : *level)
    {
        for (auto const& speciesId : basicHierarchy.getVariables())
        {
            auto const& dataId = speciesId.second;

            auto particlesData
                = std::dynamic_pointer_cast<ParticlesData<dimension>>(patch->getPatchData(dataId));


            EXPECT_EQ(0, particlesData->ghostParticles.size());

            auto patchBox = patch->getBox();

            int expectedNumberOfParticles{0};

            int const particlePerCell{2};


            // TODO : hard coded for interpOrder 1

            patchBox.coarsen(SAMRAI::hier::IntVector{SAMRAI::tbox::Dimension{dimension}, ratio});

            expectedNumberOfParticles = particlePerCell * patchBox.numberCells(dirX);


            EXPECT_EQ(expectedNumberOfParticles, particlesData->domainParticles.size());
        }
    }

    basicHierarchy.TearDown();
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);



    SAMRAI::tbox::SAMRAI_MPI::init(&argc, &argv);

    SAMRAI::tbox::SAMRAIManager::initialize();

    SAMRAI::tbox::SAMRAIManager::startup();


    int testResult = RUN_ALL_TESTS();


    // Finalize

    SAMRAI::tbox::SAMRAIManager::shutdown();

    SAMRAI::tbox::SAMRAIManager::finalize();

    SAMRAI::tbox::SAMRAI_MPI::finalize();

    return testResult;
}
