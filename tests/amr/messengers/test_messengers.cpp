#include "amr/types/amr_types.h"
#include "data/electromag/electromag.h"
#include "data/grid/gridlayout.h"
#include "data/grid/gridlayout_impl.h"
#include "data/grid/gridlayoutimplyee.h"
#include "data/ions/ion_population/ion_population.h"
#include "data/ions/ions.h"
#include "data/ions/particle_initializers/maxwellian_particle_initializer.h"
#include "data/ndarray/ndarray_vector.h"
#include "data/particles/particle_array.h"
#include "data/vecfield/vecfield.h"
#include "data_provider.h"
#include "hybrid/hybrid_quantities.h"
#include "messenger_registration.h"
#include "messengers/hybrid_messenger.h"
#include "messengers/messenger_factory.h"
#include "physical_models/hybrid_model.h"
#include "physical_models/mhd_model.h"
#include "resources_manager/resources_manager.h"
#include "solvers/solver_mhd.h"
#include "solvers/solver_ppc.h"
#include "test_basichierarchy.h"
#include "test_integrator_strat.h"
#include "test_tag_strategy.h"



#include "amr/types/amr_types.h"
#include <SAMRAI/algs/TimeRefinementIntegrator.h>
#include <SAMRAI/algs/TimeRefinementLevelStrategy.h>
#include <SAMRAI/geom/CartesianGridGeometry.h>
#include <SAMRAI/hier/CoarsenOperator.h>
#include <SAMRAI/hier/RefineOperator.h>
#include <SAMRAI/mesh/ChopAndPackLoadBalancer.h>
#include <SAMRAI/mesh/GriddingAlgorithm.h>
#include <SAMRAI/mesh/StandardTagAndInitStrategy.h>
#include <SAMRAI/mesh/StandardTagAndInitialize.h>
#include <SAMRAI/mesh/TileClustering.h>
#include <SAMRAI/tbox/InputManager.h>
#include <SAMRAI/tbox/Logger.h>
#include <SAMRAI/tbox/SAMRAIManager.h>
#include <SAMRAI/tbox/SAMRAI_MPI.h>
#include <SAMRAI/xfer/CoarsenAlgorithm.h>
#include <SAMRAI/xfer/RefineAlgorithm.h>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

using namespace PHARE::core;
using namespace PHARE::amr;
using namespace PHARE::amr;

static constexpr std::size_t dim         = 1;
static constexpr std::size_t interpOrder = 1;
using GridImplYee1D                      = GridLayoutImplYee<dim, interpOrder>;
using GridYee1D                          = GridLayout<GridImplYee1D>;
using VecField1D                         = VecField<NdArrayVector1D<>, HybridQuantity>;
using MaxwellianParticleInitializer1D
    = MaxwellianParticleInitializer<ParticleArray<dim>, GridYee1D>;
using IonsPop1D         = IonPopulation<ParticleArray<dim>, VecField1D, GridYee1D>;
using Ions1D            = Ions<IonsPop1D, GridYee1D>;
using Electromag1D      = Electromag<VecField1D>;
using HybridModelT      = HybridModel<GridYee1D, Electromag1D, Ions1D, SAMRAI_Types>;
using MHDModelT         = MHDModel<GridYee1D, VecField1D, SAMRAI_Types>;
using ResourcesManagerT = ResourcesManager<GridYee1D>;
using hybhybStratT      = HybridHybridMessengerStrategy<HybridModelT, IPhysicalModel<SAMRAI_Types>>;
using mhdhybStratT
    = MHDHybridMessengerStrategy<MHDModelT, HybridModelT, IPhysicalModel<SAMRAI_Types>>;
using HybridMessengerT = HybridMessenger<HybridModelT, IPhysicalModel<SAMRAI_Types>>;



double density(double x)
{
    return x * +2.;
}

double vx(double x)
{
    (void)x;
    return 1.;
}


double vy(double x)
{
    (void)x;
    return 1.;
}


double vz(double x)
{
    (void)x;
    return 1.;
}


double vthx(double x)
{
    (void)x;
    return 1.;
}


double vthy(double x)
{
    (void)x;
    return 1.;
}



double vthz(double x)
{
    (void)x;
    return 1.;
}



double bx(double x)
{
    return x + 2.;
}

double by(double x)
{
    return 2; // x + 2.;
}

double bz(double x)
{
    return x + 4.;
}

double ex(double x)
{
    return x + 5.;
}

double ey(double x)
{
    return x + 6.;
}

double ez(double x)
{
    return x + 7.;
}



using ScalarFunctionT = PHARE::initializer::ScalarFunction<1>;

PHARE::initializer::PHAREDict<1> createIonsDict()
{
    PHARE::initializer::PHAREDict<1> dict;
    dict["ions"]["name"]           = std::string{"ions"};
    dict["ions"]["nbrPopulations"] = std::size_t{2};
    dict["ions"]["pop0"]["name"]   = std::string{"protons"};
    dict["ions"]["pop0"]["mass"]   = 1.;
    dict["ions"]["pop0"]["ParticleInitializer"]["name"]
        = std::string{"MaxwellianParticleInitializer"};
    dict["ions"]["pop0"]["ParticleInitializer"]["density"] = static_cast<ScalarFunctionT>(density);

    dict["ions"]["pop0"]["ParticleInitializer"]["bulk_velocity_x"]
        = static_cast<ScalarFunctionT>(vx);

    dict["ions"]["pop0"]["ParticleInitializer"]["bulk_velocity_y"]
        = static_cast<ScalarFunctionT>(vy);

    dict["ions"]["pop0"]["ParticleInitializer"]["bulk_velocity_z"]
        = static_cast<ScalarFunctionT>(vz);


    dict["ions"]["pop0"]["ParticleInitializer"]["thermal_velocity_x"]
        = static_cast<ScalarFunctionT>(vthx);

    dict["ions"]["pop0"]["ParticleInitializer"]["thermal_velocity_y"]
        = static_cast<ScalarFunctionT>(vthy);

    dict["ions"]["pop0"]["ParticleInitializer"]["thermal_velocity_z"]
        = static_cast<ScalarFunctionT>(vthz);


    dict["ions"]["pop0"]["ParticleInitializer"]["nbrPartPerCell"] = std::size_t{100};
    dict["ions"]["pop0"]["ParticleInitializer"]["charge"]         = -1.;
    dict["ions"]["pop0"]["ParticleInitializer"]["basis"]          = std::string{"Cartesian"};

    dict["ions"]["pop1"]["name"] = std::string{"alpha"};
    dict["ions"]["pop1"]["mass"] = 1.;
    dict["ions"]["pop1"]["ParticleInitializer"]["name"]
        = std::string{"MaxwellianParticleInitializer"};
    dict["ions"]["pop1"]["ParticleInitializer"]["density"] = static_cast<ScalarFunctionT>(density);

    dict["ions"]["pop1"]["ParticleInitializer"]["bulk_velocity_x"]
        = static_cast<ScalarFunctionT>(vx);

    dict["ions"]["pop1"]["ParticleInitializer"]["bulk_velocity_y"]
        = static_cast<ScalarFunctionT>(vy);

    dict["ions"]["pop1"]["ParticleInitializer"]["bulk_velocity_z"]
        = static_cast<ScalarFunctionT>(vz);


    dict["ions"]["pop1"]["ParticleInitializer"]["thermal_velocity_x"]
        = static_cast<ScalarFunctionT>(vthx);

    dict["ions"]["pop1"]["ParticleInitializer"]["thermal_velocity_y"]
        = static_cast<ScalarFunctionT>(vthy);

    dict["ions"]["pop1"]["ParticleInitializer"]["thermal_velocity_z"]
        = static_cast<ScalarFunctionT>(vthz);


    dict["ions"]["pop1"]["ParticleInitializer"]["nbrPartPerCell"] = std::size_t{100};
    dict["ions"]["pop1"]["ParticleInitializer"]["charge"]         = -1.;
    dict["ions"]["pop1"]["ParticleInitializer"]["basis"]          = std::string{"Cartesian"};

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



struct AfullHybridBasicHierarchy : public ::testing::Test
{
    int const firstHybLevel{0};
    int const ratio{2};
    short unsigned const dimension = 1;

    using HybridHybridT = HybridHybridMessengerStrategy<HybridModelT, IPhysicalModel<SAMRAI_Types>>;



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


TEST_F(AfullHybridBasicHierarchy, initializesFieldsOnRefinedLevels)
{
    auto& hierarchy = basicHierarchy->getHierarchy();
    auto& rm        = hybridModel->resourcesManager;

    for (auto iLevel = 0; iLevel < hierarchy.getNumberOfLevels(); ++iLevel)
    {
        auto const& level = hierarchy.getPatchLevel(iLevel);

        for (auto& patch : *level)
        {
            auto onPatch
                = rm->setOnPatch(*patch, hybridModel->state.electromag, hybridModel->state.ions);

            auto layout = layoutFromPatch<typename HybridModelT::gridLayout_type>(*patch);

            auto& Ex = hybridModel->state.electromag.E.getComponent(Component::X);
            auto& Ey = hybridModel->state.electromag.E.getComponent(Component::Y);
            auto& Ez = hybridModel->state.electromag.E.getComponent(Component::Z);
            auto& Bx = hybridModel->state.electromag.B.getComponent(Component::X);
            auto& By = hybridModel->state.electromag.B.getComponent(Component::Y);
            auto& Bz = hybridModel->state.electromag.B.getComponent(Component::Z);

            auto checkMyField = [&layout](auto const& field, auto const& func) //
            {
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
        }
    }
}


TEST_F(AfullHybridBasicHierarchy, initializesParticlesOnRefinedLevels)
{
    auto& hierarchy = basicHierarchy->getHierarchy();
    auto& rm        = hybridModel->resourcesManager;

    for (auto iLevel = 0; iLevel < hierarchy.getNumberOfLevels(); ++iLevel)
    {
        auto const& level = hierarchy.getPatchLevel(iLevel);

        std::cout << "iLevel = " << iLevel << "\n";

        for (auto& patch : *level)
        {
            auto onPatch = rm->setOnPatch(*patch, hybridModel->state.ions);
            auto& ions   = hybridModel->state.ions;

            for (auto& pop : ions)
            {
                auto& patchGhosts          = pop.patchGhostParticles();
                auto& oldLevelBorderGhosts = pop.levelGhostParticlesOld();

                auto geom              = patch->getPatchGeometry();
                auto const& boundaries = geom->getPatchBoundaries();
                EXPECT_GT(patchGhosts.size(), 0);

                // domain particles
                EXPECT_GT(pop.nbrParticles(), 0);

                // here we expect to have level border particles only for
                // refined levels and for a patch that has boundaries//
                // note that here "boundaries" refers to both physical boundaries and
                // coarse to fine boundaries. So technically a patch could be on a refined
                // level and have 'boundaries' without having any coarse-to-fine ones but only
                // physical and the test would fail. However here since we are periodic, there are
                // no physical boundaries wo we're ok.
                if (iLevel > 0 && boundaries[0].size() > 0)
                {
                    EXPECT_GT(oldLevelBorderGhosts.size(), 0);
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




TEST_F(AfullHybridBasicHierarchy, fillsRefinedLevelFieldGhosts)
{
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

class StreamAppender : public SAMRAI::tbox::Logger::Appender
{
public:
    StreamAppender(std::ostream* stream) { d_stream = stream; }
    void logMessage(const std::string& message, const std::string& filename, const int line)
    {
        (*d_stream) << "At :" << filename << " line :" << line << " message: " << message
                    << std::endl;
    }

private:
    std::ostream* d_stream;
};



int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    SAMRAI::tbox::SAMRAI_MPI::init(&argc, &argv);
    SAMRAI::tbox::SAMRAIManager::initialize();
    SAMRAI::tbox::SAMRAIManager::startup();

    int world_size, world_rank, name_len;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Get_processor_name(processor_name, &name_len);

    SAMRAI::tbox::SAMRAI_MPI::setCallAbortInParallelInsteadOfMPIAbort(true);

    std::shared_ptr<SAMRAI::tbox::Logger::Appender> appender
        = std::make_shared<StreamAppender>(StreamAppender{&std::cout});
    SAMRAI::tbox::Logger::getInstance()->setWarningAppender(appender);

    int testResult = RUN_ALL_TESTS();

    // Finalize
    SAMRAI::tbox::SAMRAIManager::shutdown();
    SAMRAI::tbox::SAMRAIManager::finalize();
    SAMRAI::tbox::SAMRAI_MPI::finalize();

    return testResult;
}
