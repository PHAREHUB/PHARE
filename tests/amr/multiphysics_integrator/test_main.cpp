
#include <mpi.h>

#include "tests/simulator/per_test.h"

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


class ParticleSplit : public SAMRAI::hier::RefineOperator
{
public:
    ParticleSplit()
        : SAMRAI::hier::RefineOperator{"ParticlesSplit"}
    {
    }

    virtual ~ParticleSplit() = default;

    int getOperatorPriority() const override { return 0; }

    SAMRAI::hier::IntVector getStencilWidth(SAMRAI::tbox::Dimension const& dim) const override
    {
        return SAMRAI::hier::IntVector::getOne(dim);
    }

    void refine(SAMRAI::hier::Patch& destination, SAMRAI::hier::Patch const& source,
                int const destinationComponent, int const sourceComponent,
                SAMRAI::hier::BoxOverlap const& destinationOverlap,
                SAMRAI::hier::IntVector const& ratio) const override
    {
        // do nothing
    }
};

struct DisabledType
{
};

class Schedule
{
public:
};

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

#if 0
template<typename MHDMessenger, typename MHDToHybrid, typename HybridMessenger,
         typename HybridToFullPIC, typename FullPICMessenger>
class MessengerManager
{
public:
    template<typename Variable>
    void registerVariable(Variable const &variable)
    {
        std::cout << "register variable for future messenger\n";
        // this one will consist at creating the Algorithm for each
        // component.
        // As well as tracking which variable need which algorithm
        // Note that at this moment we do not know the temporal nature
        // of the messenger : so we have two possibilites :
        // create one algorithm for each or adding a parameter to determine
        // if we will need temporal interpolation
    }

private:
    std::vector<MHDMessenger> mhdMessengers_;
    std::vector<HybridToFullPIC> hybridMessengers_;
    std::vector<FullPICMessenger> fullPICMessenger_;

    MHDToHybrid mhdToHybridMessenger_;
    HybridToFullPIC hybridToFullPICMessenger_;
};
#endif


// -----------------------------------------------------------------------------
//                          MULTIPHYSICS INTEGRATOR
// -----------------------------------------------------------------------------


TYPED_TEST(SimulatorTest, knowsWhichSolverisOnAGivenLevel)
{
    TypeParam sim;
    auto& multiphysInteg = *sim.getMultiPhysicsIntegrator();

    for (int iLevel = 0; iLevel < sim.getNumberOfLevels(); ++iLevel)
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
    auto& hierarchy   = *sim.getPrivateHierarchy();
    auto& hybridModel = *sim.getHybridModel();
    auto& mhdModel    = *sim.getMHDModel();

    for (int iLevel = 0; iLevel < hierarchy.getNumberOfLevels(); ++iLevel)
    {
        if (isInMHDdRange(iLevel))
        {
            auto Bid = mhdModel.resourcesManager->getIDs(mhdModel.state.B);
            auto Vid = mhdModel.resourcesManager->getIDs(mhdModel.state.V);

            std::array<std::vector<int> const*, 2> allIDs{{&Bid, &Vid}};

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
        else if (isInHybridRange(iLevel))
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

    auto nbrOfLevels = sim.getNumberOfLevels();
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




TYPED_TEST(SimulatorTest, returnsCorrecMessengerForEachLevel)
{
    TypeParam sim;
    auto& multiphysInteg = *sim.getMultiPhysicsIntegrator();

    // EXPECT_EQ(std::string{"MHDModel-MHDModel"}, multiphysInteg.messengerName(0));
    // EXPECT_EQ(std::string{"MHDModel-MHDModel"}, multiphysInteg.messengerName(1));
    // EXPECT_EQ(std::string{"MHDModel-HybridModel"}, multiphysInteg.messengerName(2));
    // EXPECT_EQ(std::string{"HybridModel-HybridModel"}, multiphysInteg.messengerName(3));
    for (int i = 0; i < sim.getNumberOfLevels(); i++)
        EXPECT_EQ(std::string{"HybridModel-HybridModel"}, multiphysInteg.messengerName(i));
}




/*
#endif

*/

#if 0
TEST(aMessenger, isCreated)
{
    SAMRAI::tbox::Dimension dimension{dim};
    ResourcesManager<GridYee1D> resourcesManager{dimension};


    IonsInit1D ionsInit;
    ionsInit.name   = "Ions";
    ionsInit.masses = {{0.1, 0.3}};

    ionsInit.names.emplace_back("specie1");
    ionsInit.names.emplace_back("specie2");

    ionsInit.nbrPopulations = 2;

    ionsInit.particleInitializers.push_back(std::make_unique<FluidParticleInitializer1D>(
        std::make_unique<ScalarFunction<dim>>(density),
        std::make_unique<VectorFunction<dim>>(bulkVelocity),
        std::make_unique<VectorFunction<dim>>(thermalVelocity), -1., 10));

    ionsInit.particleInitializers.push_back(std::make_unique<FluidParticleInitializer1D>(
        std::make_unique<ScalarFunction<dim>>(density),
        std::make_unique<VectorFunction<dim>>(bulkVelocity),
        std::make_unique<VectorFunction<dim>>(thermalVelocity), -1., 10));


    // HybridState need an ions initializers to create an ions
    HybridState<Electromag1D, Ions1D, IonsInit1D> hybridModel{std::move(ionsInit)};

    SolverPPC<Electromag1D, decltype(hybridModel)> solverppc{hybridModel};

    // we have to register the hybridModel for variable instanciation on  a patch
    resourcesManager.registerResources(hybridModel.electromag.E);
    resourcesManager.registerResources(hybridModel.electromag.B);
    resourcesManager.registerResources(hybridModel.ions);



    // our solver may also have data that need to  be instanciate on a patch
    solverppc.init(resourcesManager);


    auto fieldRefineOp  = std::make_shared<FieldDataLinearRefine<GridImplYee1D, Field1D>>();
    auto fieldCoarsenOp = std::make_shared<FieldDataCoarsen<GridImplYee1D, Field1D>>();

    auto particlesRefineOp = std::make_shared<ParticleSplit>();

    int const hybridStartLevel = 0;

    // then we init the messenger needed
    HybridMessenger hybridMessenger{hybridModel,    resourcesManager,  fieldRefineOp,
                                        fieldCoarsenOp, particlesRefineOp, hybridStartLevel};

    solverppc.initMessenger(hybridMessenger);


    // hierarchy
    auto inputDatabase = SAMRAI::tbox::InputManager::getManager()->parseInputFile(
        inputBase + "input/input_1d_ratio_2.txt");
    auto patchHierarchyDatabase = inputDatabase->getDatabase("PatchHierarchy");

    auto gridGeometry = std::make_shared<SAMRAI::geom::CartesianGridGeometry>(
        dimension, "cartesian", inputDatabase->getDatabase("CartesianGridGeometry"));

    auto hierarchy = std::make_shared<SAMRAI::hier::PatchHierarchy>("PatchHierarchy", gridGeometry,
                                                                    patchHierarchyDatabase);

    auto loadBalancer = std::make_shared<SAMRAI::mesh::ChopAndPackLoadBalancer>(
        dimension, "ChopAndPackLoadBalancer",
        inputDatabase->getDatabase("ChopAndPackLoadBalancer"));


    auto integratorStrategy = std::make_shared<
        MultiPhysicsIntegrator<decltype(hybridModel), decltype(resourcesManager)>>(
        hybridModel, hybridStartLevel, resourcesManager);

    auto standardTag = std::make_shared<SAMRAI::mesh::StandardTagAndInitialize>(
        "StandardTagAndInitialize", integratorStrategy.get(),
        inputDatabase->getDatabase("StandardTagAndInitialize"));

    auto clustering = std::make_shared<SAMRAI::mesh::TileClustering>(
        dimension, inputDatabase->getDatabase("TileClustering"));

    auto gridding = std::make_shared<SAMRAI::mesh::GriddingAlgorithm>(
        hierarchy, "GriddingAlgorithm", inputDatabase->getDatabase("GriddingAlgorithm"),
        standardTag, clustering, loadBalancer);


    auto integrator = std::make_shared<SAMRAI::algs::TimeRefinementIntegrator>(
        "TimeRefinementIntegrator", inputDatabase->getDatabase("TimeRefinementIntegrator"),
        hierarchy, integratorStrategy, gridding);



    // First step is to make the CoarsestLevel, and then to refine
    // as much as needed
    gridding->makeCoarsestLevel(0.0);
    // we refine for tag 0
    // For that we try to refine until we get no more refine boxes or that
    // we have reach the maximum level
    for (int iLevel = 0, hierarchyMaxLevelAllowed = hierarchy.getMaxNumberOfLevels();
         iLevel < hierarchyMaxLevelAllowed - 1; ++iLevel)
    {
        int const cycle   = 0;
        double const time = 0.0;

        SAMRAI::hier::BoxContainer boxes{};
        auto reset = standardTag->getUserSuppliedRefineBoxes(boxes, iLevel, cycle, time);
        NULL_USE(reset);
        if (!boxes.empty())
        {
            gridding->makeFinerLevel(0, true, cycle, time, 0.0);
        }
        else
        {
            break;
        }
    }


    auto level0 = sim.getPatchLevel(0);
    hybridMessenger.setStatus(&hierarchy, &level0, nullptr);

    // then for each level, we give an hybridMessenger that know which level is the current one
    solverppc.advanceLevel(hybridMessenger, hybridMessenger);

    // after advancing each substep of a leaf
    solverppc.syncLevel(hybridMessenger);
}
#endif


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    PHARE_test::SamraiLifeCycle samsam(argc, argv);
    return RUN_ALL_TESTS();
}
