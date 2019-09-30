
#ifndef PHARE_TEST_DIAGNOSTIC_HIERARCHY
#define PHARE_TEST_DIAGNOSTIC_HIERARCHY

class BasicHierarchy
{
public:
    /**
     * @brief Construct the hierarchy consisting of two levels
     *
     * Only the ratio argument is required, we will use one input file per ratio.
     * Here we do multiple things, first we parse the inputFile to create the inputDatabase
     * that will be used by others SAMRAI's objects (GridGeometry, PatchHierarchy, LoadBalancer,
     * BoxGenerator, GriddingAlgorithm, StandardTagAndInitialize), they contain information for
     * each of them.
     *
     * We then store a view to the VariableDatabase singleton , that we will use to register all
     * our Variables. For each hybridQuantities we register a variable.
     *
     * We instantiate each objects needed by the GriddingAlgorithm, and then instantiate
     * a GriddingAlgorithm with the previous objects
     *
     */
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

    std::unique_ptr<HybridMessengerStrategy<HybridModelT, IPhysicalModel>> hybhybStrat{
        std::make_unique<HybridHybridT>(resourcesManagerHybrid, firstHybLevel)};

    std::shared_ptr<HybridMessengerT> messenger{
        std::make_shared<HybridMessengerT>(std::move(hybhybStrat))};

    std::shared_ptr<SolverPPC<HybridModelT>> solver{std::make_shared<SolverPPC<HybridModelT>>()};

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
};

#endif /*PHARE_TEST_DIAGNOSTIC_HIERARCHY*/
