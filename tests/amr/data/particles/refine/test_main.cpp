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



template<std::size_t dimension_, std::size_t interpOrder_, int ratio_>
struct Descriptors
{
    static constexpr std::size_t dimension   = dimension_;
    static constexpr std::size_t interpOrder = interpOrder_;
    static constexpr int ratio               = ratio_;
};

template<typename T>
struct AParticlesDataSplitOnCoarseBoundary : public ::testing::Test
{
    BasicHierarchy<T::dimension, T::interpOrder, ParticlesDataSplitType::coarseBoundary>
        basicHierarchy{T::ratio};
};

TYPED_TEST_CASE_P(AParticlesDataSplitOnCoarseBoundary);

TYPED_TEST_P(AParticlesDataSplitOnCoarseBoundary,
             canBeUsedWithAPatchLevelFillPatternToFillOnlyCoarseBoundary)
{
    auto constexpr dimension   = TypeParam::dimension;
    auto constexpr ratio       = TypeParam::ratio;
    auto constexpr interpOrder = TypeParam::interpOrder;

    auto& hierarchy = this->basicHierarchy.getHierarchy();

    auto level = hierarchy.getPatchLevel(1);

    for (auto& patch : *level)
    {
        for (auto const& speciesId : this->basicHierarchy.getVariables())
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


            auto const& ghostCellWidth = particlesData->getGhostCellWidth();

            auto const& ghostBox = particlesData->getGhostBox();

            constexpr std::size_t nbrBound{2};

            // TODO : hard coded for interpOrder 1

            if constexpr (interpOrder == 1)
            {
                expectedNumberOfParticles = nbrBound * (0 + 1 + 1 + 0);
            }
            else if constexpr (interpOrder == 2)
            {
                expectedNumberOfParticles = nbrBound * (0 + 1 + 1 + 1 + 1);
            }
            else if constexpr (interpOrder == 3)
            {
                expectedNumberOfParticles = nbrBound * (1 + 1 + 1 + 1);
            }

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

    this->basicHierarchy.TearDown();
}

REGISTER_TYPED_TEST_CASE_P(AParticlesDataSplitOnCoarseBoundary,
                           canBeUsedWithAPatchLevelFillPatternToFillOnlyCoarseBoundary);

typedef testing::Types<Descriptors<1, 1, 2>, Descriptors<1, 2, 2>, Descriptors<1, 3, 2>> MyTypes;

INSTANTIATE_TYPED_TEST_CASE_P(TestWithOrderFrom1To3That, AParticlesDataSplitOnCoarseBoundary,
                              MyTypes);



TEST(AParticlesDataSplitDataOperator, canBeUsedWithAPatchLevelFillPatternToFillOnlyInterior)
{
    //
    std::size_t constexpr dimension{1};
    std::size_t constexpr interpOrder{1};
    int const ratio = 2;


    BasicHierarchy<dimension, interpOrder, ParticlesDataSplitType::interior> basicHierarchy{ratio};


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

            // each particles split in two , but on the boundary we have one that left
            // it is compensate by the ghost that split one in
            expectedNumberOfParticles = 2 * particlePerCell * patchBox.numberCells(dirX);


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
