
#include "core/def/phare_mpi.hpp"

#include "tests/simulator/per_test.hpp"

using namespace PHARE::core;
using namespace PHARE::amr;
using namespace PHARE::solver;

static constexpr int hybridStartLevel = 0; // levels : hybrid hybrid
static constexpr int maxLevelNbr      = 4;

bool isInHybridRange(int iLevel)
{
    return iLevel >= hybridStartLevel && iLevel < maxLevelNbr;
}
bool isInMHDdRange(int iLevel)
{
    return iLevel >= 0 && iLevel < hybridStartLevel;
}



/** \brief Algorithm purpose is to store an algorithm and the associated schedule for a
 * particular variable
 * TODO : putting refine and coarsen algorithm in the same class is not really good,
 *        given how it is used, never it will contain both a refine and a coarsen schedule.
 *        It just looks like a variant: sometimes it is a refine , and another time it is a coarsen
 *
 */
class Algorithm
{
public:
    explicit Algorithm(SAMRAI::tbox::Dimension const& dimension)
        : coarsen{dimension}
    {
    }

    std::shared_ptr<SAMRAI::hier::RefineOperator> const* refOperator{nullptr};
    std::shared_ptr<SAMRAI::hier::CoarsenOperator> const* coarseOperator{nullptr};

    SAMRAI::xfer::RefineAlgorithm refine;
    SAMRAI::xfer::CoarsenAlgorithm coarsen;

    std::vector<std::shared_ptr<SAMRAI::xfer::RefineSchedule>> refineSchedule;
    std::vector<std::shared_ptr<SAMRAI::xfer::CoarsenSchedule>> coarsenSchedule;
};



// -----------------------------------------------------------------------------
//                          MULTIPHYSICS INTEGRATOR
// -----------------------------------------------------------------------------



TYPED_TEST(SimulatorTest, knowsWhichSolverIsOnAGivenLevel)
{
    TypeParam sim;
    auto& multiphysInteg = *sim.getMultiPhysicsIntegrator();

    for (int iLevel = 0; iLevel < sim.hierarchy->getNumberOfLevels(); ++iLevel)
    {
        if (isInHybridRange(iLevel))
        {
            EXPECT_EQ(std::string{"PPC"}, multiphysInteg.solverName(iLevel));
        }
        else if (isInMHDdRange(iLevel))
        {
            EXPECT_EQ(std::string{"MHDSolver"}, multiphysInteg.solverName(iLevel));
        }
    }
}



TYPED_TEST(SimulatorTest, allocatesModelDataOnAppropriateLevels)
{
    TypeParam sim;
    auto& hierarchy   = *sim.hierarchy;
    auto& hybridModel = *sim.getHybridModel();
    // auto& mhdModel    = *sim.getMHDModel();
    //
    for (int iLevel = 0; iLevel < hierarchy.getNumberOfLevels(); ++iLevel)
    {
        //     if (isInMHDdRange(iLevel))
        //     {
        //         auto Bid = mhdModel.resourcesManager->getIDs(mhdModel.state.B);
        //         auto Vid = mhdModel.resourcesManager->getIDs(mhdModel.state.V);
        //
        //         std::array<std::vector<int> const*, 2> allIDs{{&Bid, &Vid}};
        //
        //         for (auto& idVec : allIDs)
        //         {
        //             for (auto& id : *idVec)
        //             {
        //                 auto level = hierarchy.getPatchLevel(iLevel);
        //                 auto patch = level->begin();
        //                 EXPECT_TRUE(patch->checkAllocated(id));
        //             }
        //         }
        //     }
        /*else*/ if (isInHybridRange(iLevel))
        {
            auto Bid   = hybridModel.resourcesManager->getIDs(hybridModel.state.electromag.B);
            auto Eid   = hybridModel.resourcesManager->getIDs(hybridModel.state.electromag.E);
            auto IonId = hybridModel.resourcesManager->getIDs(hybridModel.state.ions);

            std::array<std::vector<int> const*, 3> allIDs{{&Bid, &Eid, &IonId}};

            for (auto& idVec : allIDs)
            {
                for (auto& id : *idVec)
                {
                    auto level = hierarchy.getPatchLevel(iLevel);
                    auto patch = level->begin();
                    EXPECT_TRUE(patch->checkAllocated(id));
                }
            }
        }
    }
}


TYPED_TEST(SimulatorTest, knowsWhichModelIsSolvedAtAGivenLevel)
{
    TypeParam sim;
    auto& multiphysInteg = *sim.getMultiPhysicsIntegrator();

    auto nbrOfLevels = sim.hierarchy->getNumberOfLevels();
    for (int iLevel = 0; iLevel < nbrOfLevels; ++iLevel)
    {
        if (isInMHDdRange(iLevel))
        {
            EXPECT_EQ(std::string{"MHDModel"}, multiphysInteg.modelName(iLevel));
        }
        else if (isInHybridRange(iLevel))
        {
            EXPECT_EQ(std::string{"HybridModel"}, multiphysInteg.modelName(iLevel));
        }
    }
}




TYPED_TEST(SimulatorTest, returnsCorrectMessengerForEachLevel)
{
    TypeParam sim;
    auto& multiphysInteg = *sim.getMultiPhysicsIntegrator();

    // EXPECT_EQ(std::string{"MHDModel-MHDModel"}, multiphysInteg.messengerName(0));
    // EXPECT_EQ(std::string{"MHDModel-MHDModel"}, multiphysInteg.messengerName(1));
    // EXPECT_EQ(std::string{"MHDModel-HybridModel"}, multiphysInteg.messengerName(2));
    // EXPECT_EQ(std::string{"HybridModel-HybridModel"}, multiphysInteg.messengerName(3));
    for (int i = 0; i < sim.hierarchy->getNumberOfLevels(); i++)
        EXPECT_EQ(std::string{"HybridModel-HybridModel"}, multiphysInteg.messengerName(i));
}




int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    PHARE::SamraiLifeCycle samsam(argc, argv);
    return RUN_ALL_TESTS();
}
