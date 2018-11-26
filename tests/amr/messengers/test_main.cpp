#include "data/electromag/electromag.h"
#include "data/grid/gridlayout.h"
#include "data/grid/gridlayout_impl.h"
#include "data/grid/gridlayoutimplyee.h"
#include "data/ions/ion_population/ion_population.h"
#include "data/ions/ions.h"
#include "data/ions/particle_initializers/fluid_particle_initializer.h"
#include "data/ndarray/ndarray_vector.h"
#include "data/particles/particle_array.h"
#include "data/vecfield/vecfield.h"
#include "evolution/messengers/hybrid_messenger.h"
#include "evolution/messengers/messenger_factory.h"
#include "evolution/messengers/messenger_initializer.h"
#include "evolution/solvers/solver_mhd.h"
#include "evolution/solvers/solver_ppc.h"
#include "hybrid/hybrid_quantities.h"
#include "physical_models/hybrid_model.h"
#include "test_basichierarchy.h"
#include "test_integrator_strat.h"
#include "test_tag_strategy.h"
#include "tools/resources_manager.h"



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
#include <SAMRAI/tbox/SAMRAIManager.h>
#include <SAMRAI/tbox/SAMRAI_MPI.h>
#include <SAMRAI/xfer/CoarsenAlgorithm.h>
#include <SAMRAI/xfer/RefineAlgorithm.h>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <iostream>
#include <map>
#include <memory>
#include <vector>



using namespace PHARE;


static constexpr std::size_t dim         = 1;
static constexpr std::size_t interpOrder = 1;
using VecField1D                         = VecField<NdArrayVector1D<>, HybridQuantity>;
using IonsPop1D                          = IonPopulation<ParticleArray<dim>, VecField1D>;
using GridImplYee1D                      = GridLayoutImplYee<dim, interpOrder>;
using GridYee1D                          = GridLayout<GridImplYee1D>;
using Ions1D                             = Ions<IonsPop1D, GridYee1D>;
using Electromag1D                       = Electromag<VecField1D>;
using IonsInit1D                         = IonsInitializer<ParticleArray<dim>, GridYee1D>;
using FluidParticleInitializer1D         = FluidParticleInitializer<ParticleArray<dim>, GridYee1D>;
using HybridModelT                       = HybridModel<GridYee1D, Electromag1D, Ions1D, IonsInit1D>;
using MHDModelT                          = MHDModel<GridYee1D, VecField1D>;
using ResourcesManagerT                  = ResourcesManager<GridYee1D>;
using hybhybStratT                       = HybridHybridMessengerStrategy<HybridModelT>;
using mhdhybStratT                       = MHDHybridMessengerStrategy<MHDModelT, HybridModelT>;
using HybridMessengerT                   = HybridMessenger<HybridModelT>;



double density(double x)
{
    return x * x + 2.;
}

std::array<double, 3> bulkVelocity(double x)
{
    return std::array<double, 3>{{1.0, 0.0, 0.0}};
}


std::array<double, 3> thermalVelocity(double x)
{
    return std::array<double, 3>{{0.5, 0.0, 0.0}};
}

#if 1

TEST(anHybridMessenger, hasACorrectName)
{
    auto getIonsInit_ = []() {
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

        return ionsInit;
    };

    auto hyb_rm = std::make_shared<ResourcesManagerT>();
    auto mhd_rm = std::make_shared<ResourcesManagerT>();

    HybridModelT hybridModel{getIonsInit_(), hyb_rm};
    MHDModelT mhdModel{hyb_rm};

    auto hybhybstrat = std::make_unique<hybhybStratT>(hyb_rm, 0);

    std::unique_ptr<IMessenger> ht = std::make_unique<HybridMessengerT>(std::move(hybhybstrat));
    EXPECT_EQ(std::string{"HybridModel-HybridModel"}, ht->name());


    auto mhdstrat = std::make_unique<mhdhybStratT>(mhd_rm, hyb_rm, 0);
    ht            = std::make_unique<HybridMessenger<decltype(hybridModel)>>(std::move(mhdstrat));
    EXPECT_EQ(std::string{"MHDModel-HybridModel"}, ht->name());

    auto db = SAMRAI::hier::VariableDatabase::getDatabase();


    db->removeVariable("HybridModel-HybridModel_EM_old_E_x");
    db->removeVariable("HybridModel-HybridModel_EM_old_E_y");
    db->removeVariable("HybridModel-HybridModel_EM_old_E_z");

    db->removeVariable("MHDModel-HybridModel_EM_old_E_x");
    db->removeVariable("MHDModel-HybridModel_EM_old_E_y");
    db->removeVariable("MHDModel-HybridModel_EM_old_E_z");

    db->removeVariable("HybridModel-HybridModel_EM_old_B_x");
    db->removeVariable("HybridModel-HybridModel_EM_old_B_y");
    db->removeVariable("HybridModel-HybridModel_EM_old_B_z");

    db->removeVariable("MHDModel-HybridModel_EM_old_B_x");
    db->removeVariable("MHDModel-HybridModel_EM_old_B_y");
    db->removeVariable("MHDModel-HybridModel_EM_old_B_z");
}

#endif



auto getIonsInit() // TODO refactor this getIonInit used in several tests
{
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

    return ionsInit;
}



class AMessengerInitializerWithHybridAndMHDModels : public ::testing::Test
{
    std::vector<MessengerDescriptor> descriptors{
        {"MHDModel", "MHDModel"}, {"MHDModel", "HybridModel"}, {"HybridModel", "HybridModel"}};
    MessengerFactory<MHDModelT, HybridModelT> messengerFactory{descriptors};


public:
    std::vector<std::unique_ptr<IMessenger>> messengers;
    std::vector<std::unique_ptr<IPhysicalModel>> models;

    AMessengerInitializerWithHybridAndMHDModels()
    {
        auto resourcesManagerHybrid = std::make_shared<ResourcesManagerT>();
        auto resourcesManagerMHD    = std::make_shared<ResourcesManagerT>();

        auto hybridModel = std::make_unique<HybridModelT>(getIonsInit(), resourcesManagerHybrid);
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

    MessengerInitializer initializer;

    virtual ~AMessengerInitializerWithHybridAndMHDModels()
    {
        auto db = SAMRAI::hier::VariableDatabase::getDatabase();

        auto& mhdModel    = dynamic_cast<MHDModelT&>(*models[0]);
        auto& hybridModel = dynamic_cast<HybridModelT&>(*models[1]);

        for (auto vecFieldProperty : mhdModel.state.B.getFieldNamesAndQuantities())
        {
            db->removeVariable(vecFieldProperty.name);
        }
        for (auto vecFieldProperty : mhdModel.state.V.getFieldNamesAndQuantities())
        {
            db->removeVariable(vecFieldProperty.name);
        }

        for (auto vecFieldProperty : hybridModel.state.electromag.B.getFieldNamesAndQuantities())
        {
            db->removeVariable(vecFieldProperty.name);
        }
        for (auto vecFieldProperty : hybridModel.state.electromag.E.getFieldNamesAndQuantities())
        {
            db->removeVariable(vecFieldProperty.name);
        }

        auto rhoName = extractNames(hybridModel.state.ions);
        for (auto const& name : rhoName)
            db->removeVariable(name);

        auto bulktuple = hybridModel.state.ions.getCompileTimeResourcesUserList();
        auto& bulk     = std::get<0>(bulktuple);

        for (auto vecFieldProperty : bulk.getFieldNamesAndQuantities())
        {
            db->removeVariable(vecFieldProperty.name);
        }


        for (auto& pop : hybridModel.state.ions)
        {
            auto namesQties = pop.getFieldNamesAndQuantities();
            for (auto nameQty : namesQties)
            {
                auto name = nameQty.name;
                db->removeVariable(name);
            }

            auto bulktuple = pop.getCompileTimeResourcesUserList();
            auto& bulk     = std::get<0>(bulktuple);

            for (auto vecFieldProperty : bulk.getFieldNamesAndQuantities())
            {
                db->removeVariable(vecFieldProperty.name);
            }


            auto namesParticleArrays = pop.getParticleArrayNames();
            for (auto partProp : namesParticleArrays)
            {
                db->removeVariable(partProp.name);
            }

            db->removeVariable("HybridModel-HybridModel_EM_old_E_x");
            db->removeVariable("HybridModel-HybridModel_EM_old_E_y");
            db->removeVariable("HybridModel-HybridModel_EM_old_E_z");

            db->removeVariable("MHDModel-HybridModel_EM_old_E_x");
            db->removeVariable("MHDModel-HybridModel_EM_old_E_y");
            db->removeVariable("MHDModel-HybridModel_EM_old_E_z");

            db->removeVariable("HybridModel-HybridModel_EM_old_B_x");
            db->removeVariable("HybridModel-HybridModel_EM_old_B_y");
            db->removeVariable("HybridModel-HybridModel_EM_old_B_z");

            db->removeVariable("MHDModel-HybridModel_EM_old_B_x");
            db->removeVariable("MHDModel-HybridModel_EM_old_B_y");
            db->removeVariable("MHDModel-HybridModel_EM_old_B_z");
        }
    }
};



TEST_F(AMessengerInitializerWithHybridAndMHDModels, canInitializeMHDHybridMessengers)
{
    auto hybridSolver = std::make_unique<SolverPPC<HybridModelT>>();
    auto mhdSolver    = std::make_unique<SolverMHD<MHDModelT>>();

    MessengerInitializer::setup(*messengers[1], *models[0], *models[1], *hybridSolver);
}



TEST_F(AMessengerInitializerWithHybridAndMHDModels, canInitializeMHDMessengers)
{
    auto mhdSolver = std::make_unique<SolverMHD<MHDModelT>>();
    MessengerInitializer::setup(*messengers[0], *models[0], *models[0], *mhdSolver);
}



TEST_F(AMessengerInitializerWithHybridAndMHDModels, canInitializeHybridHybridMessengers)
{
    auto hybridSolver = std::make_unique<SolverPPC<HybridModelT>>();
    MessengerInitializer::setup(*messengers[2], *models[1], *models[1], *hybridSolver);
}




struct ABasicHierarchyWithHybridMessenger : public ::testing::Test
{
    int const firstHybLevel{0};
    using HybridHybridT = HybridHybridMessengerStrategy<HybridModelT>;



    std::shared_ptr<ResourcesManagerT> resourcesManagerHybrid{
        std::make_shared<ResourcesManagerT>()};

    std::shared_ptr<HybridModelT> hybridModel{
        std::make_shared<HybridModelT>(getIonsInit(), resourcesManagerHybrid)};


    std::unique_ptr<HybridMessengerStrategy<HybridModelT>> hybhybStrat{
        std::make_unique<HybridHybridT>(resourcesManagerHybrid, firstHybLevel)};

    std::shared_ptr<HybridMessenger<HybridModelT>> messenger{
        std::make_shared<HybridMessenger<HybridModelT>>(std::move(hybhybStrat))};

    std::shared_ptr<SolverPPC<HybridModelT>> solver{std::make_shared<SolverPPC<HybridModelT>>()};




    ABasicHierarchyWithHybridMessenger()
    {
        hybridModel->resourcesManager->registerResources(hybridModel->state.electromag);
        hybridModel->resourcesManager->registerResources(hybridModel->state.ions);
        solver->registerResources(*hybridModel);

        auto fromCoarserInfo = messenger->emptyInfoFromCoarser();
        auto fromFinerInfo   = messenger->emptyInfoFromFiner();

        hybridModel->fillMessengerInfo(fromCoarserInfo);
        hybridModel->fillMessengerInfo(fromFinerInfo);

        messenger->registerQuantities(std::move(fromCoarserInfo), std::move(fromFinerInfo));
    }


    ~ABasicHierarchyWithHybridMessenger()
    {
        auto db = SAMRAI::hier::VariableDatabase::getDatabase();

        for (auto vecFieldProperty : hybridModel->state.electromag.B.getFieldNamesAndQuantities())
        {
            db->removeVariable(vecFieldProperty.name);
        }
        for (auto vecFieldProperty : hybridModel->state.electromag.E.getFieldNamesAndQuantities())
        {
            db->removeVariable(vecFieldProperty.name);
        }

        auto rhoName = extractNames(hybridModel->state.ions);
        for (auto const& name : rhoName)
            db->removeVariable(name);

        auto bulktuple = hybridModel->state.ions.getCompileTimeResourcesUserList();
        auto& bulk     = std::get<0>(bulktuple);

        for (auto vecFieldProperty : bulk.getFieldNamesAndQuantities())
        {
            db->removeVariable(vecFieldProperty.name);
        }


        for (auto& pop : hybridModel->state.ions)
        {
            auto namesQties = pop.getFieldNamesAndQuantities();
            for (auto nameQty : namesQties)
            {
                auto name = nameQty.name;
                db->removeVariable(name);
            }

            auto bulktuple = pop.getCompileTimeResourcesUserList();
            auto& bulk     = std::get<0>(bulktuple);

            for (auto vecFieldProperty : bulk.getFieldNamesAndQuantities())
            {
                db->removeVariable(vecFieldProperty.name);
            }


            auto namesParticleArrays = pop.getParticleArrayNames();
            for (auto partProp : namesParticleArrays)
            {
                db->removeVariable(partProp.name);
            }

            db->removeVariable("HybridModel-HybridModel_EM_old_E_x");
            db->removeVariable("HybridModel-HybridModel_EM_old_E_y");
            db->removeVariable("HybridModel-HybridModel_EM_old_E_z");

            db->removeVariable("HybridModel-HybridModel_EM_old_B_x");
            db->removeVariable("HybridModel-HybridModel_EM_old_B_y");
            db->removeVariable("HybridModel-HybridModel_EM_old_B_z");

            db->removeVariable("EMPred_E_x");
            db->removeVariable("EMPred_E_y");
            db->removeVariable("EMPred_E_z");
            db->removeVariable("EMPred_B_x");
            db->removeVariable("EMPred_B_y");
            db->removeVariable("EMPred_B_z");

            db->removeVariable("EMAvg_E_x");
            db->removeVariable("EMAvg_E_y");
            db->removeVariable("EMAvg_E_z");
            db->removeVariable("EMAvg_B_x");
            db->removeVariable("EMAvg_B_y");
            db->removeVariable("EMAvg_B_z");
        }
    }
};




TEST_F(ABasicHierarchyWithHybridMessenger, initializesRefinedLevels)
{
    auto tagStrat   = std::make_shared<TagStrategy<HybridModelT>>(hybridModel, solver, messenger);
    int const ratio = 2;
    short unsigned const dimension = 1;

    auto integratorStrat = std::make_shared<TestIntegratorStrat>();

    BasicHierarchy basicHierarchy{ratio, dimension, tagStrat.get(), integratorStrat};

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




TEST_F(ABasicHierarchyWithHybridMessenger, initializesNewLevelDuringRegrid)
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




TEST_F(ABasicHierarchyWithHybridMessenger, initializesNewFinestLevelAfterRegrid)
{
    auto tagStrat   = std::make_shared<TagStrategy<HybridModelT>>(hybridModel, solver, messenger);
    int const ratio = 2;
    short unsigned const dimension = 1;

    auto integratorStrat = std::make_shared<TestIntegratorStrat>();


    //   BasicHierarchy hierarchy{ratio, dimension, tagStrat.get(),integratorStrat};
}




TEST_F(ABasicHierarchyWithHybridMessenger, fillsRefinedLevelGhosts)
{
    auto tagStrat   = std::make_shared<TagStrategy<HybridModelT>>(hybridModel, solver, messenger);
    int const ratio = 2;
    short unsigned const dimension = 1;
    auto newTime                   = 1.;
    auto integratorStrat           = std::make_shared<TestIntegratorStrat>();


    BasicHierarchy basicHierarchy{ratio, dimension, tagStrat.get(), integratorStrat};



#if 0
    auto& hierarchy = basicHierarchy.getHierarchy();

    auto const& level = hierarchy.getPatchLevel(2);

    for (auto& patch : *level)
    {
        auto _ = hybridModel->resourcesManager->makeResourcesGuard(*patch,
                                                                   hybridModel->state.electromag);

        auto layout = PHARE::layoutFromPatch<typename HybridModelT::gridLayout_type>(*patch);

        auto& Ex = hybridModel->state.electromag.E.getComponent(PHARE::Component::X);
        auto& Ey = hybridModel->state.electromag.E.getComponent(PHARE::Component::Y);
        auto& Ez = hybridModel->state.electromag.E.getComponent(PHARE::Component::Z);

        auto& Bx = hybridModel->state.electromag.B.getComponent(PHARE::Component::X);
        auto& By = hybridModel->state.electromag.B.getComponent(PHARE::Component::Y);
        auto& Bz = hybridModel->state.electromag.B.getComponent(PHARE::Component::Z);


        auto pseudoAdvanceMyField = [&layout](auto& field, auto const& func) //
        {
            auto iStart = layout.physicalStartIndex(field, PHARE::Direction::X);
            auto iEnd   = layout.physicalEndIndex(field, PHARE::Direction::X);

            for (auto ix = iStart; ix <= iEnd; ++ix)
            {
                auto origin = layout.origin();
                auto x      = layout.fieldNodeCoordinates(field, origin, ix);
                field(ix)   = func(x[0]) + 2;
                // auto expected = func(x[0]);
                // EXPECT_DOUBLE_EQ(expected, field(ix));
            }
        };

        pseudoAdvanceMyField(Ex, TagStrategy<HybridModelT>::fillInt);
        pseudoAdvanceMyField(Ey, TagStrategy<HybridModelT>::fillInt);
        pseudoAdvanceMyField(Ez, TagStrategy<HybridModelT>::fillInt);
        pseudoAdvanceMyField(Bx, TagStrategy<HybridModelT>::fillInt);
        pseudoAdvanceMyField(By, TagStrategy<HybridModelT>::fillInt);
        pseudoAdvanceMyField(Bz, TagStrategy<HybridModelT>::fillInt);

        hybridModel->resourcesManager->setTime(hybridModel->state.electromag, *patch, newTime);
    }


    auto const& level3 = hierarchy.getPatchLevel(3);
    auto const& level2 = hierarchy.getPatchLevel(2);

    for (auto patch : *level2)
    {
        auto exOldId = hybridModel->resourcesManager->getID("EM_old_E_x");
        auto exId    = hybridModel->resourcesManager->getID("EM_E_x");


        ASSERT_TRUE(exOldId);
        ASSERT_TRUE(exId);

        EXPECT_TRUE(patch->checkAllocated(*exOldId));
        EXPECT_TRUE(patch->checkAllocated(*exId));

        ASSERT_TRUE(patch->checkAllocated(*exOldId));
        ASSERT_TRUE(patch->checkAllocated(*exId));


        auto exOldData = patch->getPatchData(*exOldId);
        auto exData    = patch->getPatchData(*exId);



        EXPECT_DOUBLE_EQ(0., patch->getPatchData(*exOldId)->getTime());
        EXPECT_DOUBLE_EQ(1., patch->getPatchData(*exId)->getTime());
    }


    /*
        for (auto& patch : *level3)
        {
            hybridModel->resourcesManager->setTime(hybridModel->state.electromag, *patch,
                                                   newTime * 0.5);
        }
    */

    messenger->fillElectricGhosts(hybridModel->state.electromag.E, 3, newTime * 0.5);
    messenger->fillMagneticGhosts(hybridModel->state.electromag.B, 3, newTime * 0.5);


    auto iLevel = 3;



    std::cout << "iLevel = " << iLevel << "\n";




    for (auto& patch : *level3)
    {
        auto _ = hybridModel->resourcesManager->makeResourcesGuard(*patch,
                                                                   hybridModel->state.electromag);

        auto layout = PHARE::layoutFromPatch<typename HybridModelT::gridLayout_type>(*patch);

        auto& Ex = hybridModel->state.electromag.E.getComponent(PHARE::Component::X);
        auto& Ey = hybridModel->state.electromag.E.getComponent(PHARE::Component::Y);
        auto& Ez = hybridModel->state.electromag.E.getComponent(PHARE::Component::Z);

        auto& Bx = hybridModel->state.electromag.B.getComponent(PHARE::Component::X);
        auto& By = hybridModel->state.electromag.B.getComponent(PHARE::Component::Y);
        auto& Bz = hybridModel->state.electromag.B.getComponent(PHARE::Component::Z);



        auto patchBox = patch->getBox();


        auto checkMyFieldGhosts = [&layout](auto const& field, auto const& func) //
        {



            auto iStart = layout.ghostStartIndex(
                field, PHARE::Direction::X); // physicalStartIndex(field, PHARE::Direction::X);
            auto iEnd = layout.physicalStartIndex(field, PHARE::Direction::X);

            for (auto ix = iStart; ix < iEnd; ++ix)
            {
                auto origin   = layout.origin();
                auto x        = layout.fieldNodeCoordinates(field, origin, ix);
                auto expected = 0.5 * (1 * func(x[0]) + 2.);
                std::cout << "original  = " << func(x[0]) << " actual = " << field(ix)
                          << " expected = " << expected << "\n";
                EXPECT_DOUBLE_EQ(expected, field(ix));
            }


            iStart = layout.physicalEndIndex(
                field, PHARE::Direction::X); // physicalStartIndex(field, PHARE::Direction::X);
            iEnd = layout.ghostEndIndex(field, PHARE::Direction::X);

            for (auto ix = iStart + 1; ix <= iEnd; ++ix)
            {
                auto origin   = layout.origin();
                auto x        = layout.fieldNodeCoordinates(field, origin, ix);
                auto expected = 0.5 * (1 * func(x[0]) + 2.);
                EXPECT_DOUBLE_EQ(expected, field(ix));
            }
        };


        checkMyFieldGhosts(Ex, TagStrategy<HybridModelT>::fillInt);
        checkMyFieldGhosts(Ey, TagStrategy<HybridModelT>::fillInt);
        checkMyFieldGhosts(Ez, TagStrategy<HybridModelT>::fillInt);

        checkMyFieldGhosts(Bx, TagStrategy<HybridModelT>::fillInt);
        checkMyFieldGhosts(By, TagStrategy<HybridModelT>::fillInt);
        checkMyFieldGhosts(Bz, TagStrategy<HybridModelT>::fillInt);
    }

#endif
}




TEST_F(ABasicHierarchyWithHybridMessenger, fillsRefinedLevelGhostsAfterRegrid)
{
    auto tagStrat = std::make_shared<TagStrategy<HybridModelT>>(hybridModel, solver, messenger);

    int const ratio                = 2;
    short unsigned const dimension = 1;

    auto integratorStrat = std::make_shared<TestIntegratorStrat>();


    //   BasicHierarchy hierarchy{ratio, dimension, tagStrat.get(), integratorStrat};
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
