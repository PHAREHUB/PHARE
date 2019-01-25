#include "data/particles/particles_data.h"
#include <SAMRAI/hier/RefineOperator.h>
#include <SAMRAI/tbox/SAMRAIManager.h>
#include <SAMRAI/tbox/SAMRAI_MPI.h>
#include <SAMRAI/xfer/RefineAlgorithm.h>


#include "test_basic_hierarchy.h"


#include "gmock/gmock.h"
#include "gtest/gtest.h"



#include <iostream>

class ALevelWithDomainParticles : public ::testing::Test
{
public:
    ALevelWithDomainParticles()
        : hierarchy{basicHierarchy.getHierarchy()}
        , level{hierarchy.getPatchLevel(0)}
    {
        algo.registerRefine(0, 0, 0, std::shared_ptr<SAMRAI::hier::RefineOperator>{});
        schedule = algo.createSchedule(level);
    }

    BasicHierarchy<1, 1> basicHierarchy;
    SAMRAI::hier::PatchHierarchy& hierarchy;
    std::shared_ptr<SAMRAI::hier::PatchLevel> level;
    SAMRAI::xfer::RefineAlgorithm algo;
    std::shared_ptr<SAMRAI::xfer::RefineSchedule> schedule;


    ~ALevelWithDomainParticles()
    {
        auto db = SAMRAI::hier::VariableDatabase::getDatabase();
        db->removeVariable("protons");
    }
};



TEST_F(ALevelWithDomainParticles, ghostParticleNumberisZeroBeforeScheduleAndCorrectAfterSchedule)
{
    for (auto const& patch : *level)
    {
        auto pdat    = patch->getPatchData(0);
        auto partDat = std::dynamic_pointer_cast<ParticlesData<1>>(pdat);
        auto& ghosts = partDat->ghostParticles;
        EXPECT_EQ(0u, ghosts.size());

        std::cout << "patch has " << ghosts.size() << " ghost particles\n";
    }

    std::cout << "filling ghosts\n";
    schedule->fillData(0.);


    for (auto const& patch : *level)
    {
        auto pdat    = patch->getPatchData(0);
        auto partDat = std::dynamic_pointer_cast<ParticlesData<1>>(pdat);
        auto& ghosts = partDat->ghostParticles;
        EXPECT_EQ(6, ghosts.size());
        std::cout << "patch has " << ghosts.size() << " ghost particles\n";
    }
}




TEST_F(ALevelWithDomainParticles, hasGhostParticlesInGhostAMRCellsAfterSchedule)
{
    schedule->fillData(0.);

    for (auto const& patch : *level)
    {
        auto pdat     = patch->getPatchData(0);
        auto partDat  = std::dynamic_pointer_cast<ParticlesData<1>>(pdat);
        auto ghostBox = partDat->getGhostBox();
        auto lower    = ghostBox.lower();
        auto upper    = ghostBox.upper();

        auto nbrLowerParticles  = 0;
        auto nbrUpperParticles  = 0;
        auto anormalParticleNbr = 0;
        for (auto const& particle : partDat->ghostParticles)
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




TEST_F(ALevelWithDomainParticles, hasGhostParticleFilledFromNeighborFirstCell)
{
    for (auto const& patch : *level)
    {
        auto pdat       = patch->getPatchData(0);
        auto partDat    = std::dynamic_pointer_cast<ParticlesData<1>>(pdat);
        auto& particles = partDat->domainParticles;
        auto& ghosts    = partDat->ghostParticles;
        std::cout << "patch has " << particles.size() << " domain particles "
                  << "and " << ghosts.size() << " ghost particles\n";
    }

    schedule->fillData(0.);

    for (auto const& patch : *level)
    {
        auto pdat       = patch->getPatchData(0);
        auto partDat    = std::dynamic_pointer_cast<ParticlesData<1>>(pdat);
        auto& particles = partDat->domainParticles;
        auto& ghosts    = partDat->ghostParticles;
        std::cout << "patch has " << particles.size() << " domain particles "
                  << "and " << ghosts.size() << " ghost particles\n";
    }
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
