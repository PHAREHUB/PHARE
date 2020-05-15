#include "amr/data/particles/particles_data.h"
#include <SAMRAI/hier/RefineOperator.h>
#include <SAMRAI/tbox/SAMRAIManager.h>
#include <SAMRAI/tbox/SAMRAI_MPI.h>
#include <SAMRAI/xfer/RefineAlgorithm.h>


#include "test_basic_hierarchy.h"


#include "gmock/gmock.h"
#include "gtest/gtest.h"


#include <iostream>

template<size_t dim>
class ALevelWithDomainParticles
{
public:
    ALevelWithDomainParticles()
        : hierarchy{basicHierarchy.getHierarchy()}
        , level{hierarchy.getPatchLevel(0)}
    {
        algo.registerRefine(0, 0, 0, std::shared_ptr<SAMRAI::hier::RefineOperator>{});
        schedule = algo.createSchedule(level);
    }

    BasicHierarchy<dim, 1> basicHierarchy;
    SAMRAI::hier::PatchHierarchy& hierarchy;
    std::shared_ptr<SAMRAI::hier::PatchLevel> level;
    SAMRAI::xfer::RefineAlgorithm algo;
    std::shared_ptr<SAMRAI::xfer::RefineSchedule> schedule;


    ~ALevelWithDomainParticles()
    {
        auto db = SAMRAI::hier::VariableDatabase::getDatabase();
        db->removeVariable("protons");
    }

    void ghostParticleNumberisZeroBeforeScheduleAndCorrectAfterSchedule();
    void hasGhostParticlesInGhostAMRCellsAfterSchedule();
    void hasGhostParticleFilledFromNeighborFirstCell();
};
template<typename Simulator>
struct LevelWithDomainParticlesTest : public ::testing::Test
{
};
using LevelWithDomainParticlesNd
    = testing::Types<ALevelWithDomainParticles<1>, ALevelWithDomainParticles<2>>;
TYPED_TEST_SUITE(LevelWithDomainParticlesTest, LevelWithDomainParticlesNd);

template<size_t dim>
void ALevelWithDomainParticles<
    dim>::ghostParticleNumberisZeroBeforeScheduleAndCorrectAfterSchedule()
{
    ASSERT_TRUE(level->getNumberOfPatches());

    for (auto const& patch : *level)
    {
        auto pdat    = patch->getPatchData(0);
        auto partDat = std::dynamic_pointer_cast<ParticlesData<dim>>(pdat);
        auto& ghosts = partDat->patchGhostParticles;
        EXPECT_EQ(0u, ghosts.size());

        std::cout << "patch has " << ghosts.size() << " ghost particles\n";
    }

    std::cout << "filling ghosts\n";
    schedule->fillData(0.);

    using Cells = size_t;
    using Sides = size_t;

    constexpr size_t PPC                = 3;
    constexpr size_t ghost_particles_1d = PPC * Cells{1} * Sides{2};
    constexpr size_t ghost_particles_2d = PPC * (Cells{6} * Sides{2} + Cells{4} * Sides{2});

    constexpr size_t dim_particles[] = {ghost_particles_1d, ghost_particles_2d};
    constexpr size_t eq              = dim_particles[dim - 1];

    for (auto const& patch : *level)
    {
        auto pdat    = patch->getPatchData(0);
        auto partDat = std::dynamic_pointer_cast<ParticlesData<dim>>(pdat);
        auto& ghosts = partDat->patchGhostParticles;

        // 1d expects 3 particles in each (left/right) ghost cell
        EXPECT_EQ(eq, ghosts.size());
        std::cout << "patch has " << ghosts.size() << " ghost particles\n";
    }
}

TYPED_TEST(LevelWithDomainParticlesTest,
           ghostParticleNumberisZeroBeforeScheduleAndCorrectAfterSchedule)

{
    TypeParam{}.ghostParticleNumberisZeroBeforeScheduleAndCorrectAfterSchedule();
}


template<size_t dim>
void ALevelWithDomainParticles<dim>::hasGhostParticlesInGhostAMRCellsAfterSchedule()
{
    ASSERT_TRUE(level->getNumberOfPatches());

    schedule->fillData(0.);

    for (auto const& patch : *level)
    {
        auto pdat     = patch->getPatchData(0);
        auto partDat  = std::dynamic_pointer_cast<ParticlesData<dim>>(pdat);
        auto ghostBox = partDat->getGhostBox();
        auto lower    = ghostBox.lower();
        auto upper    = ghostBox.upper();

        auto nbrLowerParticles  = 0;
        auto nbrUpperParticles  = 0;
        auto anormalParticleNbr = 0;
        for (auto const& particle : partDat->patchGhostParticles)
        {
            if (particle.iCell[0] == lower[0])
                nbrLowerParticles++;
            else if (particle.iCell[0] == upper[0])
                nbrUpperParticles++;
            else
                anormalParticleNbr++;
        }
        EXPECT_EQ(3, nbrLowerParticles);
        EXPECT_EQ(3, nbrUpperParticles);
        EXPECT_EQ(0, anormalParticleNbr);
    }
}

TYPED_TEST(LevelWithDomainParticlesTest, hasGhostParticlesInGhostAMRCellsAfterSchedule)
{
    TypeParam{}.ghostParticleNumberisZeroBeforeScheduleAndCorrectAfterSchedule();
}


template<size_t dim>
void ALevelWithDomainParticles<dim>::hasGhostParticleFilledFromNeighborFirstCell()
{
    ASSERT_TRUE(level->getNumberOfPatches());

    for (auto const& patch : *level)
    {
        auto pdat       = patch->getPatchData(0);
        auto partDat    = std::dynamic_pointer_cast<ParticlesData<dim>>(pdat);
        auto& particles = partDat->domainParticles;
        auto& ghosts    = partDat->patchGhostParticles;
        std::cout << "patch has " << particles.size() << " domain particles "
                  << "and " << ghosts.size() << " ghost particles\n";
    }

    schedule->fillData(0.);

    for (auto const& patch : *level)
    {
        auto pdat       = patch->getPatchData(0);
        auto partDat    = std::dynamic_pointer_cast<ParticlesData<dim>>(pdat);
        auto& particles = partDat->domainParticles;
        auto& ghosts    = partDat->patchGhostParticles;
        std::cout << "patch has " << particles.size() << " domain particles "
                  << "and " << ghosts.size() << " ghost particles\n";
    }
}

TYPED_TEST(LevelWithDomainParticlesTest, hasGhostParticleFilledFromNeighborFirstCell)
{
    TypeParam{}.hasGhostParticleFilledFromNeighborFirstCell();
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
