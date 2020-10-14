#include "core/data/particles/particle.h"
#include "amr/data/particles/particles_data.h"
#include "amr/data/particles/refine/split.h"
#include "test_particledata_refine_basic_hierarchy.h"
#include "test_particle_data_refine_tag_strategy.h"


#include "gmock/gmock.h"
#include "gtest/gtest.h"


#include <algorithm>
#include <limits>
#include <vector>


#include <SAMRAI/tbox/SAMRAI_MPI.h>


using testing::DoubleNear;
using testing::Eq;

template<std::size_t dimension_, std::size_t interpOrder_, std::size_t refineParticlesNbr_,
         ParticlesDataSplitType SplitType_>
struct ParticlesDataSplitTestDescriptors
{
    std::size_t constexpr static ratio                = 2;
    std::size_t constexpr static dimension            = dimension_;
    std::size_t constexpr static interpOrder          = interpOrder_;
    std::size_t constexpr static refineParticlesNbr   = refineParticlesNbr_;
    static constexpr ParticlesDataSplitType SplitType = SplitType_;
};




template<typename Type>
class aSimpleBasicHierarchyWithTwoLevels : public ::testing::Test
{
public:
    static auto constexpr dimension                   = Type::dimension;
    static auto constexpr interpOrder                 = Type::interpOrder;
    static auto constexpr refineParticlesNbr          = Type::refineParticlesNbr;
    static auto constexpr ratio                       = Type::ratio;
    static constexpr ParticlesDataSplitType SplitType = Type::SplitType;

    using Splitter = PHARE::amr::Splitter<PHARE::core::DimConst<dimension>,
                                          PHARE::core::InterpConst<interpOrder>,
                                          RefinedParticlesConst<refineParticlesNbr>>;

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
    SAMRAI::tbox::SAMRAI_MPI mpi{MPI_COMM_WORLD};


    ~aSimpleBasicHierarchyWithTwoLevels()
    {
        auto db = SAMRAI::hier::VariableDatabase::getDatabase();
        db->removeVariable("proton1");
    }



    ParticleArray<dimension> getRefinedL0Particles()
    {
        ParticleArray<dimension> refinedParticles;

        Splitter split;

        auto geom        = this->hierarchy.getGridGeometry();
        auto domainBoxes = geom->getPhysicalDomain();

        auto domainBox  = domainBoxes.front();
        auto ghostWidth = static_cast<int>(ghostWidthForParticles<interpOrder>());

        auto lowerXYZ = boxBoundsLower<dimension>(domainBox);
        auto upperXYZ = boxBoundsUpper<dimension>(domainBox);

        for (std::size_t i = 0; i < dimension; i++)
        {
            lowerXYZ[i] -= ghostWidth;
            upperXYZ[i] += ghostWidth;
        }


        for (int iCellX = lowerXYZ[dirX]; iCellX <= upperXYZ[dirX]; ++iCellX)
        {
            for (int iCellY = lowerXYZ[dirY]; iCellY <= upperXYZ[dirY]; ++iCellY)
            {
                for (int iCellZ = lowerXYZ[dirZ]; iCellZ <= upperXYZ[dirZ]; ++iCellZ)
                {
                    auto coarseParticles = loadCell<dimension>(iCellX, iCellY, iCellZ);

                    for (auto const& part : coarseParticles)
                    {
                        auto coarseOnRefinedGrid{part};

                        for (auto iDim = 0u; iDim < dimension; ++iDim)
                        {
                            auto delta = static_cast<int>(part.delta[iDim] * ratio);
                            coarseOnRefinedGrid.iCell[iDim] = part.iCell[iDim] * ratio + delta;
                            coarseOnRefinedGrid.delta[iDim] = part.delta[iDim] * ratio - delta;
                        }

                        refinedParticles.resize(refinedParticles.size() + refineParticlesNbr);
                        split(coarseOnRefinedGrid, refinedParticles,
                              refinedParticles.size() - refineParticlesNbr);
                    }
                }
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
    using ParticlesData_t                  = ParticlesData<ParticleArray<dimension>>;

    ParticleArray<dimension>
    filterLevelGhostParticles(ParticleArray<dimension> const& refinedParticles,
                              std::shared_ptr<SAMRAI::hier::Patch> const& patch)
    {
        ParticleArray<dimension> levelGhosts;

        auto pData    = patch->getPatchData(this->dataID);
        auto partData = std::dynamic_pointer_cast<ParticlesData_t>(pData);

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
        auto partDat = std::dynamic_pointer_cast<ParticlesData_t>(pDat);

        return partDat->levelGhostParticles;
    }
};




TYPED_TEST_SUITE_P(levelOneCoarseBoundaries);




TYPED_TEST_P(levelOneCoarseBoundaries, areCorrectlyFilledByRefinedSchedule)
{
    // this test can only work with one patch on L1
    // since parallel execution will make more
    // disable the test for parallel runs.
    if (this->mpi.getSize() > 1)
    {
        GTEST_SKIP();
    }
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




template<typename Type>
class levelOneInterior : public aSimpleBasicHierarchyWithTwoLevels<Type>
{
public:
    using BaseType                         = aSimpleBasicHierarchyWithTwoLevels<Type>;
    static constexpr std::size_t dimension = BaseType::dimension;
    using ParticlesData_t                  = ParticlesData<ParticleArray<dimension>>;

    ParticleArray<dimension>
    filterInteriorParticles(ParticleArray<dimension> const& refinedParticles,
                            std::shared_ptr<SAMRAI::hier::Patch> const& patch)
    {
        ParticleArray<dimension> interiors;

        auto pData    = patch->getPatchData(this->dataID);
        auto partData = std::dynamic_pointer_cast<ParticlesData_t>(pData);

        auto myBox = partData->getBox();

        std::copy_if(std::begin(refinedParticles), std::end(refinedParticles),
                     std::back_inserter(interiors),
                     [&myBox](auto const& part) { return isInBox(myBox, part); });

        return interiors;
    }




    auto& getInterior(std::shared_ptr<SAMRAI::hier::Patch> const& patch)
    {
        auto pDat    = patch->getPatchData(this->dataID);
        auto partDat = std::dynamic_pointer_cast<ParticlesData_t>(pDat);

        return partDat->domainParticles;
    }
};




TYPED_TEST_SUITE_P(levelOneInterior);




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
REGISTER_TYPED_TEST_SUITE_P(levelOneCoarseBoundaries, areCorrectlyFilledByRefinedSchedule);




// std::size_t constexpr dimension   = 1;
// std::size_t constexpr interpOrder = 1;
// int constexpr refinedParticlesNbr = 2;

// dimension , ratio , interpOrder, refinedParticlesNbr
template<std::size_t dim, std::size_t babies = 2>
using ParticlesDataSplitDescriptorsNDr2o1ref2C2F
    = ParticlesDataSplitTestDescriptors<dim, 1, babies, ParticlesDataSplitType::coarseBoundary>;

template<std::size_t dim, std::size_t babies = 2>
using ParticlesDataSplitDescriptorsNDr2o2ref2C2F
    = ParticlesDataSplitTestDescriptors<dim, 2, babies, ParticlesDataSplitType::coarseBoundary>;

template<std::size_t dim, std::size_t babies = 2>
using ParticlesDataSplitDescriptorsNDr2o3ref2C2F
    = ParticlesDataSplitTestDescriptors<dim, 3, babies, ParticlesDataSplitType::coarseBoundary>;

template<std::size_t dim, std::size_t babies = 3>
using ParticlesDataSplitDescriptorsNDr2o1ref3C2F
    = ParticlesDataSplitTestDescriptors<dim, 1, babies, ParticlesDataSplitType::coarseBoundary>;

template<std::size_t dim, std::size_t babies = 3>
using ParticlesDataSplitDescriptorsNDr2o2ref3C2F
    = ParticlesDataSplitTestDescriptors<dim, 2, babies, ParticlesDataSplitType::coarseBoundary>;

template<std::size_t dim, std::size_t babies = 3>
using ParticlesDataSplitDescriptorsNDr2o3ref3C2F
    = ParticlesDataSplitTestDescriptors<dim, 3, babies, ParticlesDataSplitType::coarseBoundary>;

typedef ::testing::Types<
    ParticlesDataSplitDescriptorsNDr2o1ref2C2F<1>, ParticlesDataSplitDescriptorsNDr2o2ref2C2F<1>,
    ParticlesDataSplitDescriptorsNDr2o3ref2C2F<1>, ParticlesDataSplitDescriptorsNDr2o1ref3C2F<1>,
    ParticlesDataSplitDescriptorsNDr2o2ref3C2F<1>, ParticlesDataSplitDescriptorsNDr2o3ref3C2F<1>,
    ParticlesDataSplitDescriptorsNDr2o1ref2C2F<2, 4>,
    ParticlesDataSplitDescriptorsNDr2o2ref2C2F<2, 4>,
    ParticlesDataSplitDescriptorsNDr2o3ref2C2F<2, 4>,
    ParticlesDataSplitDescriptorsNDr2o1ref3C2F<2, 4>,
    ParticlesDataSplitDescriptorsNDr2o2ref3C2F<2, 4>,
    ParticlesDataSplitDescriptorsNDr2o3ref3C2F<2, 4>>
    ParticlesCoarseToFineDataDescriptorsRange;

template<std::size_t dim, std::size_t babies = 2>
using ParticlesDataSplitDescriptorsNDr2o1ref2Int
    = ParticlesDataSplitTestDescriptors<dim, 1, babies, ParticlesDataSplitType::interior>;

template<std::size_t dim, std::size_t babies = 2>
using ParticlesDataSplitDescriptorsNDr2o2ref2Int
    = ParticlesDataSplitTestDescriptors<dim, 2, babies, ParticlesDataSplitType::interior>;

template<std::size_t dim, std::size_t babies = 2>
using ParticlesDataSplitDescriptorsNDr2o3ref2Int
    = ParticlesDataSplitTestDescriptors<dim, 3, babies, ParticlesDataSplitType::interior>;

template<std::size_t dim, std::size_t babies = 3>
using ParticlesDataSplitDescriptorsNDr2o1ref3Int
    = ParticlesDataSplitTestDescriptors<dim, 1, babies, ParticlesDataSplitType::interior>;

template<std::size_t dim, std::size_t babies = 3>
using ParticlesDataSplitDescriptorsNDr2o2ref3Int
    = ParticlesDataSplitTestDescriptors<dim, 2, babies, ParticlesDataSplitType::interior>;

template<std::size_t dim, std::size_t babies = 3>
using ParticlesDataSplitDescriptorsNDr2o3ref3Int
    = ParticlesDataSplitTestDescriptors<dim, 3, babies, ParticlesDataSplitType::interior>;

typedef ::testing::Types<
    ParticlesDataSplitDescriptorsNDr2o1ref2Int<1>, ParticlesDataSplitDescriptorsNDr2o2ref2Int<1>,
    ParticlesDataSplitDescriptorsNDr2o3ref2Int<1>, ParticlesDataSplitDescriptorsNDr2o1ref3Int<1>,
    ParticlesDataSplitDescriptorsNDr2o2ref3Int<1>, ParticlesDataSplitDescriptorsNDr2o3ref3Int<1>,
    ParticlesDataSplitDescriptorsNDr2o1ref2Int<2, 4>,
    ParticlesDataSplitDescriptorsNDr2o2ref2Int<2, 4>,
    ParticlesDataSplitDescriptorsNDr2o3ref2Int<2, 4>,
    ParticlesDataSplitDescriptorsNDr2o1ref3Int<2, 4>,
    ParticlesDataSplitDescriptorsNDr2o2ref3Int<2, 4>,
    ParticlesDataSplitDescriptorsNDr2o3ref3Int<2, 4>>
    ParticlesInteriorDataDescriptorsRange;



INSTANTIATE_TYPED_TEST_SUITE_P(TestCoarseToFine, levelOneCoarseBoundaries,
                               ParticlesCoarseToFineDataDescriptorsRange);



INSTANTIATE_TYPED_TEST_SUITE_P(TestInterior, levelOneInterior,
                               ParticlesInteriorDataDescriptorsRange);

// INSTANTIATE_TYPED_TEST_SUITE_P(TestInterior, levelOneCoarseBoundaries, TestTest);

namespace
{
template<std::size_t dimension, std::size_t interpOrder, std::size_t refineParticlesNbr>
using Splitter
    = PHARE::amr::Splitter<PHARE::core::DimConst<dimension>, PHARE::core::InterpConst<interpOrder>,
                           RefinedParticlesConst<refineParticlesNbr>>;

template<typename Splitter>
struct SplitterTest : public ::testing::Test
{
    SplitterTest() { Splitter splitter; }
};

using Splitters = testing::Types<Splitter<1, 1, 2>, Splitter<2, 1, 8> /*, Splitter<3, 1, 27>*/>;

TYPED_TEST_SUITE(SplitterTest, Splitters);

TYPED_TEST(SplitterTest, constexpr_init)
{
    constexpr TypeParam param{};
}

} // namespace
