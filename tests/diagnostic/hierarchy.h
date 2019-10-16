
#ifndef PHARE_TEST_DIAGNOSTIC_HIERARCHY
#define PHARE_TEST_DIAGNOSTIC_HIERARCHY

class BasicHierarchy
{
public:
    BasicHierarchy(
        short unsigned const dimension, SAMRAI::mesh::StandardTagAndInitStrategy* tagStrat,
        std::shared_ptr<SAMRAI::algs::TimeRefinementLevelStrategy> const& integratorStrat,
        std::string input)
        : inputDatabase_{SAMRAI::tbox::InputManager::getManager()->parseInputFile(input)}
        , patchHierarchyDatabase_{inputDatabase_->getDatabase("PatchHierarchy")}
        , dimension_{dimension}

        , gridGeometry_{std::make_shared<SAMRAI::geom::CartesianGridGeometry>(
              dimension_, "cartesian", inputDatabase_->getDatabase("CartesianGridGeometry"))}
        , hierarchy_{std::make_shared<SAMRAI::hier::PatchHierarchy>("PatchHierarchy", gridGeometry_,
                                                                    patchHierarchyDatabase_)}
        , loadBalancer_{std::make_shared<SAMRAI::mesh::ChopAndPackLoadBalancer>(
              dimension_, "ChopAndPackLoadBalancer",
              inputDatabase_->getDatabase("ChopAndPackLoadBalancer"))}
        , standardTag_{std::make_shared<SAMRAI::mesh::StandardTagAndInitialize>(
              "StandardTagAndInitialize", tagStrat,
              inputDatabase_->getDatabase("StandardTagAndInitialize"))}

        , clustering_{std::make_shared<SAMRAI::mesh::TileClustering>(
              dimension_, inputDatabase_->getDatabase("TileClustering"))}
        , gridding{std::make_shared<SAMRAI::mesh::GriddingAlgorithm>(
              hierarchy_, "GriddingAlgorithm", inputDatabase_->getDatabase("GriddingAlgorithm"),
              standardTag_, clustering_, loadBalancer_)}
        , integrator{std::make_shared<SAMRAI::algs::TimeRefinementIntegrator>(
              "TimeRefinementIntegrator", inputDatabase_->getDatabase("TimeRefinementIntegrator"),
              hierarchy_, integratorStrat, gridding)}
    {
        integrator->initializeHierarchy();
    }

    SAMRAI::hier::PatchHierarchy& getHierarchy() { return *hierarchy_; }

    ~BasicHierarchy() {}

private:
    std::shared_ptr<SAMRAI::tbox::Database> inputDatabase_;
    std::shared_ptr<SAMRAI::tbox::Database> patchHierarchyDatabase_;

    SAMRAI::tbox::Dimension dimension_;

    std::shared_ptr<SAMRAI::geom::CartesianGridGeometry> gridGeometry_;

    std::shared_ptr<SAMRAI::hier::PatchHierarchy> hierarchy_;
    std::shared_ptr<SAMRAI::mesh::ChopAndPackLoadBalancer> loadBalancer_;

    std::shared_ptr<SAMRAI::mesh::StandardTagAndInitialize> standardTag_;
    std::shared_ptr<SAMRAI::mesh::TileClustering> clustering_;

public:
    std::shared_ptr<SAMRAI::mesh::GriddingAlgorithm> gridding;
    std::shared_ptr<SAMRAI::algs::TimeRefinementIntegrator> integrator;
};

struct AfullHybridBasicHierarchy : public ::testing::Test
{
    uint8_t const dimension = 1;
    int const firstHybLevel{0};
    int const ratio{2};
    std::string input
        = "tests/amr/messengers/input/input_1d_ratio_" + std::to_string(ratio) + ".txt";

    std::shared_ptr<ResourcesManagerT> resourcesManagerHybrid{
        std::make_shared<ResourcesManagerT>()};

    std::shared_ptr<HybridModelT> hybridModel{
        std::make_shared<HybridModelT>(CREATE_IONS_DICT(), resourcesManagerHybrid)};


    std::unique_ptr<HybridMessengerStrategy<HybridModelT, IPhysicalModel<SAMRAI_Types>>>
        hybhybStrat{std::make_unique<HybridHybridT>(resourcesManagerHybrid, firstHybLevel)};

    std::shared_ptr<HybridMessenger<HybridModelT, IPhysicalModel<SAMRAI_Types>>> messenger{
        std::make_shared<HybridMessenger<HybridModelT, IPhysicalModel<SAMRAI_Types>>>(
            std::move(hybhybStrat))};

    std::shared_ptr<SolverPPC<HybridModelT, SAMRAI_Types>> solver{
        std::make_shared<SolverPPC<HybridModelT, SAMRAI_Types>>()};

    std::shared_ptr<TagStrategy<HybridModelT>> tagStrat;

    std::shared_ptr<IntegratorStrat> integrator;

    std::shared_ptr<BasicHierarchy> basicHierarchy;

    AfullHybridBasicHierarchy()
    {
        hybridModel->resourcesManager->registerResources(hybridModel->state.electromag);
        hybridModel->resourcesManager->registerResources(hybridModel->state.ions);
        solver->registerResources(*hybridModel);

        tagStrat   = std::make_shared<TagStrategy<HybridModelT>>(hybridModel, solver, messenger);
        integrator = std::make_shared<IntegratorStrat>();
        basicHierarchy
            = std::make_shared<BasicHierarchy>(dimension, tagStrat.get(), integrator, input);
    }
    ~AfullHybridBasicHierarchy() {}
};

#endif /*PHARE_TEST_DIAGNOSTIC_HIERARCHY*/
