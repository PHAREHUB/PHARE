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



template<std::size_t dimension_, std::size_t interpOrder_, int ratio_,
         std::size_t refinedParticlesNbr_>
struct Descriptors
{
    static constexpr std::size_t dimension           = dimension_;
    static constexpr std::size_t interpOrder         = interpOrder_;
    static constexpr int ratio                       = ratio_;
    static constexpr std::size_t refinedParticlesNbr = refinedParticlesNbr_;
};

template<typename T>
struct AParticlesDataSplitOnCoarseBoundary : public ::testing::Test
{
    BasicHierarchy<T::dimension, T::interpOrder, ParticlesDataSplitType::coarseBoundary,
                   T::refinedParticlesNbr>
        basicHierarchy{T::ratio};
};

TYPED_TEST_CASE_P(AParticlesDataSplitOnCoarseBoundary);

TYPED_TEST_P(AParticlesDataSplitOnCoarseBoundary,
             canBeUsedWithAPatchLevelFillPatternToFillOnlyCoarseBoundary)
{
    auto constexpr dimension           = TypeParam::dimension;
    auto constexpr ratio               = TypeParam::ratio;
    auto constexpr interpOrder         = TypeParam::interpOrder;
    auto constexpr refinedParticlesNbr = TypeParam::refinedParticlesNbr;

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


            if constexpr (refinedParticlesNbr == 2)
            {
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
            }
            else if constexpr (refinedParticlesNbr == 3)
            {
                if constexpr (interpOrder == 1)
                {
                    expectedNumberOfParticles = nbrBound * (0 + 2 + 1 + 0);
                }
                else if constexpr (interpOrder == 2)
                {
                    expectedNumberOfParticles = nbrBound * (0 + 1 + 2 + 2 + 1);
                }
                else if constexpr (interpOrder == 3)
                {
                    expectedNumberOfParticles = nbrBound * (1 + 2 + 2 + 1);
                }
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


// dim , interpOrder, refinementFactor, nbrRefinedParticles
typedef testing::Types<Descriptors<1, 1, 2, 2>, Descriptors<1, 2, 2, 2>, Descriptors<1, 3, 2, 2>,
                       Descriptors<1, 1, 2, 3>, Descriptors<1, 2, 2, 3>, Descriptors<1, 3, 2, 3>>
    MyTypes;

INSTANTIATE_TYPED_TEST_CASE_P(TestWithOrderFrom1To3That, AParticlesDataSplitOnCoarseBoundary,
                              MyTypes);


template<typename T>
struct AParticlesDataSplitOnInterior : public ::testing::Test
{
    BasicHierarchy<T::dimension, T::interpOrder, ParticlesDataSplitType::interior,
                   T::refinedParticlesNbr>
        basicHierarchy{T::ratio};
};

TYPED_TEST_CASE_P(AParticlesDataSplitOnInterior);

TYPED_TEST_P(AParticlesDataSplitOnInterior, canBeUsedWithAPatchLevelFillPatternToFillOnlyInterior)
{
    //
    auto constexpr dimension           = TypeParam::dimension;
    auto constexpr ratio               = TypeParam::ratio;
    auto constexpr refinedParticlesNbr = TypeParam::refinedParticlesNbr;


    auto& hierarchy = this->basicHierarchy.getHierarchy();

    auto level = hierarchy.getPatchLevel(1);




    for (auto& patch : *level)
    {
        for (auto const& speciesId : this->basicHierarchy.getVariables())
        {
            auto const& dataId = speciesId.second;

            auto particlesData
                = std::dynamic_pointer_cast<ParticlesData<dimension>>(patch->getPatchData(dataId));


            EXPECT_EQ(0, particlesData->ghostParticles.size());

            auto patchBox = patch->getBox();

            int expectedNumberOfParticles{0};

            int const particlePerCell{2};



            patchBox.coarsen(SAMRAI::hier::IntVector{SAMRAI::tbox::Dimension{dimension}, ratio});

            if constexpr (refinedParticlesNbr == 2)
            {
                // each particles split in two , but on the boundary we have one that left
                // it is compensate by the ghost that split one in
                expectedNumberOfParticles = 2 * particlePerCell * patchBox.numberCells(dirX);
            }
            else if constexpr (refinedParticlesNbr == 3)
            {
                // each particles split in three , but on the boundary we have one that left
                // it is compensate by the ghost that split one in

                expectedNumberOfParticles = 3 * particlePerCell * patchBox.numberCells(dirX);
            }

            EXPECT_EQ(expectedNumberOfParticles, particlesData->domainParticles.size());
        }
    }

    this->basicHierarchy.TearDown();
}



REGISTER_TYPED_TEST_CASE_P(AParticlesDataSplitOnInterior,
                           canBeUsedWithAPatchLevelFillPatternToFillOnlyInterior);



INSTANTIATE_TYPED_TEST_CASE_P(TestWithOrderFrom1To3That, AParticlesDataSplitOnInterior, MyTypes);

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
