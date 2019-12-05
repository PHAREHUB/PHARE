

#include "tests/simulator/per_test.h"

#include "test_basichierarchy.h"
#include "test_integrator_strat.h"
#include "test_tag_strategy.h"


using namespace PHARE::core;
using namespace PHARE::amr;
using namespace PHARE::solver;

static constexpr std::size_t dim         = 1;
static constexpr std::size_t interpOrder = 1;

using Simulator         = PHARE::Simulator<dim, interpOrder>;
using HybridModelT      = Simulator::HybridModel;
using MHDModelT         = Simulator::MHDModel;
using ResourcesManagerT = typename HybridModelT::resources_manager_type;


double density(double x)
{
    return /*x * +*/ 2.;
}

double vx(double /*x*/)
{
    return 1.;
}


double vy(double /*x*/)
{
    return 1.;
}


double vz(double /*x*/)
{
    return 1.;
}


double vthx(double /*x*/)
{
    return 1.;
}


double vthy(double /*x*/)
{
    return 1.;
}


double vthz(double /*x*/)
{
    return 1.;
}


double bx(double x)
{
    return x /* + 1.*/;
}

double by(double x)
{
    return x /* + 2.*/;
}

double bz(double x)
{
    return x /*+ 3.*/;
}

double ex(double x)
{
    return x /* + 4.*/;
}

double ey(double x)
{
    return x /* + 5.*/;
}

double ez(double x)
{
    return x /* + 6.*/;
}



using ScalarFunctionT = PHARE::initializer::ScalarFunction<1>;

PHARE::initializer::PHAREDict createIonsDict()
{
    PHARE::initializer::PHAREDict dict;
    dict["ions"]["name"]                                    = std::string{"ions"};
    dict["ions"]["nbrPopulations"]                          = int{2};
    dict["ions"]["pop0"]["name"]                            = std::string{"protons"};
    dict["ions"]["pop0"]["mass"]                            = 1.;
    dict["ions"]["pop0"]["particle_initializer"]["name"]    = std::string{"maxwellian"};
    dict["ions"]["pop0"]["particle_initializer"]["density"] = static_cast<ScalarFunctionT>(density);

    dict["ions"]["pop0"]["particle_initializer"]["bulk_velocity_x"]
        = static_cast<ScalarFunctionT>(vx);

    dict["ions"]["pop0"]["particle_initializer"]["bulk_velocity_y"]
        = static_cast<ScalarFunctionT>(vy);

    dict["ions"]["pop0"]["particle_initializer"]["bulk_velocity_z"]
        = static_cast<ScalarFunctionT>(vz);


    dict["ions"]["pop0"]["particle_initializer"]["thermal_velocity_x"]
        = static_cast<ScalarFunctionT>(vthx);

    dict["ions"]["pop0"]["particle_initializer"]["thermal_velocity_y"]
        = static_cast<ScalarFunctionT>(vthy);

    dict["ions"]["pop0"]["particle_initializer"]["thermal_velocity_z"]
        = static_cast<ScalarFunctionT>(vthz);


    dict["ions"]["pop0"]["particle_initializer"]["nbr_part_per_cell"] = int{100};
    dict["ions"]["pop0"]["particle_initializer"]["charge"]            = -1.;
    dict["ions"]["pop0"]["particle_initializer"]["basis"]             = std::string{"cartesian"};

    dict["ions"]["pop1"]["name"]                            = std::string{"alpha"};
    dict["ions"]["pop1"]["mass"]                            = 1.;
    dict["ions"]["pop1"]["particle_initializer"]["name"]    = std::string{"maxwellian"};
    dict["ions"]["pop1"]["particle_initializer"]["density"] = static_cast<ScalarFunctionT>(density);

    dict["ions"]["pop1"]["particle_initializer"]["bulk_velocity_x"]
        = static_cast<ScalarFunctionT>(vx);

    dict["ions"]["pop1"]["particle_initializer"]["bulk_velocity_y"]
        = static_cast<ScalarFunctionT>(vy);

    dict["ions"]["pop1"]["particle_initializer"]["bulk_velocity_z"]
        = static_cast<ScalarFunctionT>(vz);


    dict["ions"]["pop1"]["particle_initializer"]["thermal_velocity_x"]
        = static_cast<ScalarFunctionT>(vthx);

    dict["ions"]["pop1"]["particle_initializer"]["thermal_velocity_y"]
        = static_cast<ScalarFunctionT>(vthy);

    dict["ions"]["pop1"]["particle_initializer"]["thermal_velocity_z"]
        = static_cast<ScalarFunctionT>(vthz);


    dict["ions"]["pop1"]["particle_initializer"]["nbr_part_per_cell"] = int{100};
    dict["ions"]["pop1"]["particle_initializer"]["charge"]            = -1.;
    dict["ions"]["pop1"]["particle_initializer"]["basis"]             = std::string{"cartesian"};

    dict["electromag"]["name"]             = std::string{"EM"};
    dict["electromag"]["electric"]["name"] = std::string{"E"};
    dict["electromag"]["magnetic"]["name"] = std::string{"B"};

    dict["electromag"]["electric"]["initializer"]["x_component"] = static_cast<ScalarFunctionT>(ex);
    dict["electromag"]["electric"]["initializer"]["y_component"] = static_cast<ScalarFunctionT>(ey);
    dict["electromag"]["electric"]["initializer"]["z_component"] = static_cast<ScalarFunctionT>(ez);

    dict["electromag"]["magnetic"]["initializer"]["x_component"] = static_cast<ScalarFunctionT>(bx);
    dict["electromag"]["magnetic"]["initializer"]["y_component"] = static_cast<ScalarFunctionT>(by);
    dict["electromag"]["magnetic"]["initializer"]["z_component"] = static_cast<ScalarFunctionT>(bz);

    return dict;
}


TEST(MessengerDescriptors, areObtainedFromAModelList)
{
    auto modelList   = std::vector<std::string>{"MHDModel", "HybridModel"};
    auto descriptors = makeDescriptors(modelList);

    EXPECT_EQ(3, descriptors.size());
    EXPECT_EQ("MHDModel", descriptors[0].coarseModel);
    EXPECT_EQ("MHDModel", descriptors[0].fineModel);

    EXPECT_EQ("MHDModel", descriptors[1].coarseModel);
    EXPECT_EQ("HybridModel", descriptors[1].fineModel);

    EXPECT_EQ("HybridModel", descriptors[2].coarseModel);
    EXPECT_EQ("HybridModel", descriptors[2].fineModel);

    modelList   = std::vector<std::string>{"HybridModel"};
    descriptors = makeDescriptors(modelList);

    EXPECT_EQ(1, descriptors.size());
    EXPECT_EQ("HybridModel", descriptors[0].coarseModel);
    EXPECT_EQ("HybridModel", descriptors[0].fineModel);
}



// ----------------------------------------------------------------------------
// The tests below test that hybrid messengers (with either MHDHybrid or HybridHybrid
// strategies) can take quantities to communicate from models and solvers
// ----------------------------------------------------------------------------


class HybridMessengers : public ::testing::Test
{
    std::vector<MessengerDescriptor> descriptors{
        {"MHDModel", "MHDModel"}, {"MHDModel", "HybridModel"}, {"HybridModel", "HybridModel"}};
    MessengerFactory<MHDModelT, HybridModelT, IPhysicalModel<SAMRAI_Types>> messengerFactory{
        descriptors};


public:
    std::vector<std::unique_ptr<IMessenger<IPhysicalModel<SAMRAI_Types>>>> messengers;
    std::vector<std::unique_ptr<IPhysicalModel<SAMRAI_Types>>> models;

    HybridMessengers()
    {
        auto resourcesManagerHybrid = std::make_shared<ResourcesManagerT>();
        auto resourcesManagerMHD    = std::make_shared<ResourcesManagerT>();

        auto hybridModel = std::make_unique<HybridModelT>(createIonsDict(), resourcesManagerHybrid);
        auto mhdModel    = std::make_unique<MHDModelT>(resourcesManagerMHD);

        hybridModel->resourcesManager->registerResources(hybridModel->state.electromag);
        hybridModel->resourcesManager->registerResources(hybridModel->state.ions);

        mhdModel->resourcesManager->registerResources(mhdModel->state.B);
        mhdModel->resourcesManager->registerResources(mhdModel->state.V);

        models.push_back(std::move(mhdModel));
        models.push_back(std::move(hybridModel));


        auto mhdmhdMessenger{
            messengerFactory.create("MHDModel-MHDModel", *models[0], *models[0], 0)};
        auto mhdHybridMessenger{
            messengerFactory.create("MHDModel-HybridModel", *models[0], *models[1], 2)};
        auto hybridHybridMessenger{
            messengerFactory.create("HybridModel-HybridModel", *models[1], *models[1], 3)};

        messengers.push_back(std::move(mhdmhdMessenger));
        messengers.push_back(std::move(mhdHybridMessenger));
        messengers.push_back(std::move(hybridHybridMessenger));
    }
};



TEST_F(HybridMessengers, receiveQuantitiesFromMHDHybridModelsAndHybridSolver)
{
    auto hybridSolver = std::make_unique<SolverPPC<HybridModelT, SAMRAI_Types>>();

    MessengerRegistration::registerQuantities(*messengers[1], *models[0], *models[1],
                                              *hybridSolver);
}



TEST_F(HybridMessengers, receiveQuantitiesFromMHDHybridModelsAndMHDSolver)
{
    auto mhdSolver = std::make_unique<SolverMHD<MHDModelT, SAMRAI_Types>>();
    MessengerRegistration::registerQuantities(*messengers[0], *models[0], *models[0], *mhdSolver);
}



TEST_F(HybridMessengers, receiveQuantitiesFromHybridModelsOnlyAndHybridSolver)
{
    auto hybridSolver = std::make_unique<SolverPPC<HybridModelT, SAMRAI_Types>>();
    MessengerRegistration::registerQuantities(*messengers[2], *models[1], *models[1],
                                              *hybridSolver);
}



TEST_F(HybridMessengers, throwsIfGivenAnIncompatibleFineModel)
{
    auto hybridSolver = std::make_unique<SolverPPC<HybridModelT, SAMRAI_Types>>();

    auto& hybridhybridMessenger = *messengers[2];
    auto& mhdModel              = *models[0];
    auto& hybridModel           = *models[1];
    EXPECT_ANY_THROW(MessengerRegistration::registerQuantities(hybridhybridMessenger, mhdModel,
                                                               hybridModel, *hybridSolver));
}


TEST_F(HybridMessengers, throwsIfGivenAnIncompatibleCoarseModel)
{
    auto hybridSolver = std::make_unique<SolverPPC<HybridModelT, SAMRAI_Types>>();

    auto& hybridhybridMessenger = *messengers[2];
    auto& mhdModel              = *models[0];
    auto& hybridModel           = *models[1];
    EXPECT_ANY_THROW(MessengerRegistration::registerQuantities(hybridhybridMessenger, hybridModel,
                                                               mhdModel, *hybridSolver));
}




TEST_F(HybridMessengers, areNamedByTheirStrategyName)
{
    EXPECT_EQ(std::string{"MHDModel-MHDModel"}, messengers[0]->name());
    EXPECT_EQ(std::string{"MHDModel-HybridModel"}, messengers[1]->name());
    EXPECT_EQ(std::string{"HybridModel-HybridModel"}, messengers[2]->name());
}



// ----------------------------------------------------------------------------
//
// ----------------------------------------------------------------------------

// level 0 doesn't match due to periodicity
TYPED_TEST(SimulatorTest, initializesFieldsOnRefinedLevels)
{
    TypeParam sim;
    auto& hybridModel = *sim.getHybridModel();
    using GridLayout  = typename TypeParam::PHARETypes::GridLayout_t;

    auto visit = [&](GridLayout& layout, std::string patchID, size_t iLevel) {
        auto& Ex = hybridModel.state.electromag.E.getComponent(Component::X);
        auto& Ey = hybridModel.state.electromag.E.getComponent(Component::Y);
        auto& Ez = hybridModel.state.electromag.E.getComponent(Component::Z);
        auto& Bx = hybridModel.state.electromag.B.getComponent(Component::X);
        auto& By = hybridModel.state.electromag.B.getComponent(Component::Y);
        auto& Bz = hybridModel.state.electromag.B.getComponent(Component::Z);

        auto checkMyField = [&](auto const& field, auto const& func) {
            auto iStart = layout.physicalStartIndex(field, Direction::X);
            auto iEnd   = layout.physicalEndIndex(field, Direction::X);

            for (auto ix = iStart; ix <= iEnd; ++ix)
            {
                auto origin   = layout.origin();
                auto x        = layout.fieldNodeCoordinates(field, origin, ix);
                auto expected = func(x[0]);

                EXPECT_DOUBLE_EQ(expected, field(ix));
            }
        };

        checkMyField(Ex, ex);
        checkMyField(Ey, ey);
        checkMyField(Ez, ez);
        checkMyField(Bx, bx);
        checkMyField(By, by);
        checkMyField(Bz, bz);
    };

    sim.visitHierarchy(visit, 1, sim.getNumberOfLevels(), hybridModel);
}


// This test needs explicit access to the Hierarchy
//  TODO if/when boundary information is accessible without it, refactor
TYPED_TEST(SimulatorTest, initializesParticlesOnRefinedLevels)
{
    TypeParam sim;
    auto& hierarchy   = *sim.getPrivateHierarchy();
    auto& hybridModel = *sim.getHybridModel();
    auto& rm          = hybridModel.resourcesManager;

    for (auto iLevel = 0; iLevel < hierarchy.getNumberOfLevels(); ++iLevel)
    {
        auto const& level = hierarchy.getPatchLevel(iLevel);

        for (auto& patch : *level)
        {
            auto onPatch = rm->setOnPatch(*patch, hybridModel.state.ions);
            auto& ions   = hybridModel.state.ions;

            for (auto& pop : ions)
            {
                // domain particles
                EXPECT_GT(pop.nbrParticles(), 0);
                EXPECT_GT(pop.patchGhostParticles().size(), 0);

                // here we expect to have level border particles only for
                // refined levels and for a patch that has boundaries//
                // note that here "boundaries" refers to both physical boundaries and
                // coarse to fine boundaries. So technically a patch could be on a refined
                // level and have 'boundaries' without having any coarse-to-fine ones but only
                // physical and the test would fail. However here since we are periodic, there are
                // no physical boundaries wo we're ok.
                auto const& boundaries = patch->getPatchGeometry()->getPatchBoundaries();
                if (iLevel > 0 && boundaries[0].size() > 0)
                {
                    EXPECT_GT(pop.levelGhostParticlesOld().size(), 0);
                }
            }
        }
    }
}


#if 0
TEST_F(HybridHybridMessenger, initializesNewLevelDuringRegrid)
{
    auto tagStrat   = std::make_shared<TagStrategy<HybridModelT>>(hybridModel, solver, messenger);
    int const ratio = 2;
    short unsigned const dimension = 1;

    auto integratorStrat = std::make_shared<TestIntegratorStrat>();


    BasicHierarchy basicHierarchy{ratio, dimension, tagStrat.get(), integratorStrat};

    auto& integrator = basicHierarchy.integrator;

    // regrid all > 0

    double rootDt = 0.1;
    integrator->advanceHierarchy(rootDt);

    auto& hierarchy = basicHierarchy.getHierarchy();

    for (auto iLevel = 0; iLevel < hierarchy.getNumberOfLevels(); ++iLevel)
    {
        auto const& level = hierarchy.getPatchLevel(iLevel);

        std::cout << "iLevel = " << iLevel << "\n";

        for (auto& patch : *level)
        {
            auto _
                = hybridModel->resourcesManager->setOnPatch(*patch, hybridModel->state.electromag);

            auto layout = PHARE::layoutFromPatch<typename HybridModelT::gridLayout_type>(*patch);

            auto& Ex = hybridModel->state.electromag.E.getComponent(PHARE::Component::X);
            auto& Ey = hybridModel->state.electromag.E.getComponent(PHARE::Component::Y);
            auto& Ez = hybridModel->state.electromag.E.getComponent(PHARE::Component::Z);

            auto& Bx = hybridModel->state.electromag.B.getComponent(PHARE::Component::X);
            auto& By = hybridModel->state.electromag.B.getComponent(PHARE::Component::Y);
            auto& Bz = hybridModel->state.electromag.B.getComponent(PHARE::Component::Z);


            auto checkMyField = [&layout](auto const& field, auto const& func) //
            {
                auto iStart = layout.physicalStartIndex(field, PHARE::Direction::X);
                auto iEnd   = layout.physicalEndIndex(field, PHARE::Direction::X);

                for (auto ix = iStart; ix <= iEnd; ++ix)
                {
                    auto origin   = layout.origin();
                    auto x        = layout.fieldNodeCoordinates(field, origin, ix);
                    auto expected = func(x[0]);
                    EXPECT_DOUBLE_EQ(expected, field(ix));
                }
            };


            checkMyField(Ex, TagStrategy<HybridModelT>::fillEx);
            checkMyField(Ey, TagStrategy<HybridModelT>::fillEy);
            checkMyField(Ez, TagStrategy<HybridModelT>::fillEz);

            checkMyField(Bx, TagStrategy<HybridModelT>::fillBx);
            checkMyField(By, TagStrategy<HybridModelT>::fillBy);
            checkMyField(Bz, TagStrategy<HybridModelT>::fillBz);
        }
    }
}



TEST_F(HybridHybridMessenger, initializesNewFinestLevelAfterRegrid)
{
    auto tagStrat   = std::make_shared<TagStrategy<HybridModelT>>(hybridModel, solver, messenger);
    int const ratio = 2;
    short unsigned const dimension = 1;

    auto integratorStrat = std::make_shared<TestIntegratorStrat>();


    //   BasicHierarchy hierarchy{ratio, dimension, tagStrat.get(),integratorStrat};
}
#endif


struct AfullHybridBasicHierarchy : public ::testing::Test
{
    int const firstHybLevel{0};
    int const ratio{2};
    short unsigned const dimension = 1;

    using HybridHybridT = HybridHybridMessengerStrategy<HybridModelT, IPhysicalModel<SAMRAI_Types>>;

    SAMRAI::tbox::SAMRAI_MPI mpi{MPI_COMM_WORLD};

    std::shared_ptr<ResourcesManagerT> resourcesManagerHybrid{
        std::make_shared<ResourcesManagerT>()};

    std::shared_ptr<HybridModelT> hybridModel{
        std::make_shared<HybridModelT>(createIonsDict(), resourcesManagerHybrid)};


    std::unique_ptr<HybridMessengerStrategy<HybridModelT, IPhysicalModel<SAMRAI_Types>>>
        hybhybStrat{std::make_unique<HybridHybridT>(resourcesManagerHybrid, firstHybLevel)};

    std::shared_ptr<HybridMessenger<HybridModelT, IPhysicalModel<SAMRAI_Types>>> messenger{
        std::make_shared<HybridMessenger<HybridModelT, IPhysicalModel<SAMRAI_Types>>>(
            std::move(hybhybStrat))};

    std::shared_ptr<SolverPPC<HybridModelT, SAMRAI_Types>> solver{
        std::make_shared<SolverPPC<HybridModelT, SAMRAI_Types>>()};

    std::shared_ptr<TagStrategy<HybridModelT>> tagStrat;


    std::shared_ptr<TestIntegratorStrat> integrator;

    std::shared_ptr<BasicHierarchy> basicHierarchy;

    AfullHybridBasicHierarchy()
    {
        hybridModel->resourcesManager->registerResources(hybridModel->state.electromag);
        hybridModel->resourcesManager->registerResources(hybridModel->state.ions);
        solver->registerResources(*hybridModel);

        tagStrat   = std::make_shared<TagStrategy<HybridModelT>>(hybridModel, solver, messenger);
        integrator = std::make_shared<TestIntegratorStrat>();
        basicHierarchy
            = std::make_shared<BasicHierarchy>(ratio, dimension, tagStrat.get(), integrator);
    }
};

TEST_F(AfullHybridBasicHierarchy, fillsRefinedLevelFieldGhosts)
{
    if (mpi.getSize() > 1)
    {
        GTEST_SKIP() << "Test Broken for // execution, SHOULD BE FIXED";
    }
    auto newTime       = 1.;
    auto& hierarchy    = basicHierarchy->getHierarchy();
    auto const& level0 = hierarchy.getPatchLevel(0);
    auto const& level1 = hierarchy.getPatchLevel(1);
    auto& rm           = hybridModel->resourcesManager;


    // this prepareStep copies the current model EM into messenger EM
    messenger->prepareStep(*hybridModel, *level0);


    // here we set the level 0 at t=1, this simulates the advanceLevel
    for (auto& patch : *level0)
    {
        auto dataOnPatch = rm->setOnPatch(*patch, hybridModel->state.electromag);
        rm->setTime(hybridModel->state.electromag, *patch, newTime);
    }


    // this simulates a substep of level 1 to an intermediate time t=0.5
    for (auto& patch : *level1)
    {
        rm->setTime(hybridModel->state.electromag, *patch, 0.5);
    }


    // now we want to fill ghosts on level 1
    // this will need the space/time interpolation of level0 EM fields between
    // t=0 and t=1. The Model on level0 is at t=1 (above set time) and thanks
    // to the call to prepareStep() the messenger holds the copy of level 0 EM fields
    // at t=0. So at this point the ghosts should be filled OK at t=0.5.
    messenger->fillMagneticGhosts(hybridModel->state.electromag.B, 1, 0.5);
    messenger->fillElectricGhosts(hybridModel->state.electromag.E, 1, 0.5);



    auto iPatch = 0;
    for (auto patch : *level1)
    {
        auto exOldId = hybridModel->resourcesManager->getID("HybridModel-HybridModel_EM_old_E_x");
        auto exId    = hybridModel->resourcesManager->getID("EM_E_x");

        ASSERT_TRUE(exOldId);
        ASSERT_TRUE(exId);

        EXPECT_TRUE(patch->checkAllocated(*exOldId));
        EXPECT_TRUE(patch->checkAllocated(*exId));

        auto exOldData = patch->getPatchData(*exOldId);
        auto exData    = patch->getPatchData(*exId);

        auto dataOnPatch = rm->setOnPatch(*patch, hybridModel->state.electromag);


        EXPECT_DOUBLE_EQ(0., patch->getPatchData(*exOldId)->getTime());
        EXPECT_DOUBLE_EQ(0.5, patch->getPatchData(*exId)->getTime());

        auto layout = layoutFromPatch<typename HybridModelT::gridLayout_type>(*patch);

        auto& Ex = hybridModel->state.electromag.E.getComponent(Component::X);
        auto& Ey = hybridModel->state.electromag.E.getComponent(Component::Y);
        auto& Ez = hybridModel->state.electromag.E.getComponent(Component::Z);

        auto& Bx = hybridModel->state.electromag.B.getComponent(Component::X);
        auto& By = hybridModel->state.electromag.B.getComponent(Component::Y);
        auto& Bz = hybridModel->state.electromag.B.getComponent(Component::Z);



        // since we have not changed the fields on level0 between time t=0 and t=1
        // but just changed the time, the time interpolation at t=0.5 on level 1
        // should be 0.5*(FieldAtT0 + FieldAtT1) = 0.5*(2*FieldAtT0) = FieldAtT0
        // moreoever, since the level0 fields are linear function of space
        // the spatial interpolation on level 1 should be equal to the result of the function
        // that defined the field on level0.
        // As a consequence, if the space/time interpolation worked the field on level1
        // should be equal to the outcome of the function used on level0
        auto checkMyField = [&layout, &iPatch](auto const& field, auto const& func) //
        {
            auto iGhostStart = layout.ghostStartIndex(field, Direction::X);
            auto iStart      = layout.physicalStartIndex(field, Direction::X);
            auto iEnd        = layout.physicalEndIndex(field, Direction::X);
            auto iGhostEnd   = layout.ghostEndIndex(field, Direction::X);

            for (auto ix = iGhostStart; ix < iStart; ++ix)
            {
                auto origin   = layout.origin();
                auto x        = layout.fieldNodeCoordinates(field, origin, ix);
                auto expected = func(x[0]);
                std::cout << iPatch << " " << ix << " " << expected << " " << field(ix)
                          << expected - field(ix) << "\n";
                EXPECT_DOUBLE_EQ(expected, field(ix));
            }


            for (auto ix = iEnd; ix < iGhostEnd; ++ix)
            {
                auto origin   = layout.origin();
                auto x        = layout.fieldNodeCoordinates(field, origin, ix);
                auto expected = func(x[0]);
                std::cout << iPatch << " " << ix << " " << expected << " " << field(ix) << " "
                          << expected - field(ix) << "\n";
                EXPECT_DOUBLE_EQ(expected, field(ix));
            }
        };

        checkMyField(Bx, bx);
        checkMyField(By, by);
        checkMyField(Bz, bz);

        checkMyField(Ex, ex);
        checkMyField(Ey, ey);
        checkMyField(Ez, ez);

        iPatch++;
    }
}




#if 0
TEST_F(AfullHybridBasicHierarchy, fillsRefinedLevelGhostsAfterRegrid)
{
    auto tagStrat = std::make_shared<TagStrategy<HybridModelT>>(hybridModel, solver, messenger);

    int const ratio                = 2;
    short unsigned const dimension = 1;

    auto integratorStrat = std::make_shared<TestIntegratorStrat>();


    //   BasicHierarchy hierarchy{ratio, dimension, tagStrat.get(), integratorStrat};
}
#endif


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    PHARE_test::SamraiLifeCycle samsam(argc, argv);
    return RUN_ALL_TESTS();
}
