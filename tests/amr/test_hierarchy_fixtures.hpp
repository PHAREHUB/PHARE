#ifndef PHARE_TEST_AMR_HIERARCHY_FIXTURES_HPP
#define PHARE_TEST_AMR_HIERARCHY_FIXTURES_HPP

#include "phare_solver.hpp"
#include "amr/types/amr_types.hpp"

#include "tests/initializer/init_functions.hpp"
#include "tests/amr/messengers/test_integrator_strat.hpp"
#include "tests/amr/messengers/test_messenger_tag_strategy.hpp"


template<uint8_t dim>
using InitFunctionT = PHARE::initializer::InitFunction<dim>;


template<uint8_t dimension>
struct DimDict
{
};

template<>
struct DimDict<1>
{
    static constexpr uint8_t dim = 1;
    static void set(PHARE::initializer::PHAREDict& dict)
    {
        using namespace PHARE::initializer::test_fn::func_1d; // density/etc are here

        dict["ions"]["pop0"]["particle_initializer"]["density"]
            = static_cast<InitFunctionT<dim>>(density);

        dict["ions"]["pop0"]["particle_initializer"]["bulk_velocity_x"]
            = static_cast<InitFunctionT<dim>>(vx);

        dict["ions"]["pop0"]["particle_initializer"]["bulk_velocity_y"]
            = static_cast<InitFunctionT<dim>>(vy);

        dict["ions"]["pop0"]["particle_initializer"]["bulk_velocity_z"]
            = static_cast<InitFunctionT<dim>>(vz);


        dict["ions"]["pop0"]["particle_initializer"]["thermal_velocity_x"]
            = static_cast<InitFunctionT<dim>>(vthx);

        dict["ions"]["pop0"]["particle_initializer"]["thermal_velocity_y"]
            = static_cast<InitFunctionT<dim>>(vthy);

        dict["ions"]["pop0"]["particle_initializer"]["thermal_velocity_z"]
            = static_cast<InitFunctionT<dim>>(vthz);

        dict["ions"]["pop1"]["particle_initializer"]["density"]
            = static_cast<InitFunctionT<dim>>(density);

        dict["ions"]["pop1"]["particle_initializer"]["bulk_velocity_x"]
            = static_cast<InitFunctionT<dim>>(vx);

        dict["ions"]["pop1"]["particle_initializer"]["bulk_velocity_y"]
            = static_cast<InitFunctionT<dim>>(vy);

        dict["ions"]["pop1"]["particle_initializer"]["bulk_velocity_z"]
            = static_cast<InitFunctionT<dim>>(vz);


        dict["ions"]["pop1"]["particle_initializer"]["thermal_velocity_x"]
            = static_cast<InitFunctionT<dim>>(vthx);

        dict["ions"]["pop1"]["particle_initializer"]["thermal_velocity_y"]
            = static_cast<InitFunctionT<dim>>(vthy);

        dict["ions"]["pop1"]["particle_initializer"]["thermal_velocity_z"]
            = static_cast<InitFunctionT<dim>>(vthz);

        dict["electromag"]["magnetic"]["initializer"]["x_component"]
            = static_cast<InitFunctionT<dim>>(bx);
        dict["electromag"]["magnetic"]["initializer"]["y_component"]
            = static_cast<InitFunctionT<dim>>(by);
        dict["electromag"]["magnetic"]["initializer"]["z_component"]
            = static_cast<InitFunctionT<dim>>(bz);

        dict["simulation"]["algo"]["ion_updater"]["pusher"]["name"] = std::string{"modified_boris"};
    }
};




template<>
struct DimDict<2>
{
    static constexpr uint8_t dim = 2;
    static void set(PHARE::initializer::PHAREDict& dict)
    {
        using namespace PHARE::initializer::test_fn::func_2d; // density/etc are here
        dict["simulation"]["algo"]["pusher"]["name"] = std::string{"modified_boris"};

        dict["ions"]["pop0"]["particle_initializer"]["density"]
            = static_cast<InitFunctionT<dim>>(density);

        dict["ions"]["pop0"]["particle_initializer"]["bulk_velocity_x"]
            = static_cast<InitFunctionT<dim>>(vx);

        dict["ions"]["pop0"]["particle_initializer"]["bulk_velocity_y"]
            = static_cast<InitFunctionT<dim>>(vy);

        dict["ions"]["pop0"]["particle_initializer"]["bulk_velocity_z"]
            = static_cast<InitFunctionT<dim>>(vz);

        dict["ions"]["pop0"]["particle_initializer"]["thermal_velocity_x"]
            = static_cast<InitFunctionT<dim>>(vthx);

        dict["ions"]["pop0"]["particle_initializer"]["thermal_velocity_y"]
            = static_cast<InitFunctionT<dim>>(vthy);

        dict["ions"]["pop0"]["particle_initializer"]["thermal_velocity_z"]
            = static_cast<InitFunctionT<dim>>(vthz);

        dict["ions"]["pop1"]["particle_initializer"]["density"]
            = static_cast<InitFunctionT<dim>>(density);

        dict["ions"]["pop1"]["particle_initializer"]["bulk_velocity_x"]
            = static_cast<InitFunctionT<dim>>(vx);

        dict["ions"]["pop1"]["particle_initializer"]["bulk_velocity_y"]
            = static_cast<InitFunctionT<dim>>(vy);

        dict["ions"]["pop1"]["particle_initializer"]["bulk_velocity_z"]
            = static_cast<InitFunctionT<dim>>(vz);

        dict["ions"]["pop1"]["particle_initializer"]["thermal_velocity_x"]
            = static_cast<InitFunctionT<dim>>(vthx);

        dict["ions"]["pop1"]["particle_initializer"]["thermal_velocity_y"]
            = static_cast<InitFunctionT<dim>>(vthy);

        dict["ions"]["pop1"]["particle_initializer"]["thermal_velocity_z"]
            = static_cast<InitFunctionT<dim>>(vthz);

        dict["electromag"]["magnetic"]["initializer"]["x_component"]
            = static_cast<InitFunctionT<dim>>(bx);
        dict["electromag"]["magnetic"]["initializer"]["y_component"]
            = static_cast<InitFunctionT<dim>>(by);
        dict["electromag"]["magnetic"]["initializer"]["z_component"]
            = static_cast<InitFunctionT<dim>>(bz);

        dict["simulation"]["algo"]["ion_updater"]["pusher"]["name"] = std::string{"modified_boris"};
    }
};

template<uint8_t dimension>
PHARE::initializer::PHAREDict createDict()
{
    PHARE::initializer::PHAREDict dict;

    dict["simulation"]["algo"]["pusher"]["name"] = std::string{"modified_boris"};

    dict["simulation"]["algo"]["ohm"]["resistivity"]       = 0.0;
    dict["simulation"]["algo"]["ohm"]["hyper_resistivity"] = 0.0001;

    dict["ions"]["nbrPopulations"]                       = std::size_t{2};
    dict["ions"]["pop0"]["name"]                         = std::string{"protons"};
    dict["ions"]["pop0"]["mass"]                         = 1.;
    dict["ions"]["pop0"]["particle_initializer"]["name"] = std::string{"maxwellian"};

    dict["ions"]["pop0"]["particle_initializer"]["nbr_part_per_cell"] = int{100};
    dict["ions"]["pop0"]["particle_initializer"]["charge"]            = -1.;
    dict["ions"]["pop0"]["particle_initializer"]["basis"]             = std::string{"cartesian"};

    dict["ions"]["pop1"]["name"]                         = std::string{"alpha"};
    dict["ions"]["pop1"]["mass"]                         = 1.;
    dict["ions"]["pop1"]["particle_initializer"]["name"] = std::string{"maxwellian"};

    dict["ions"]["pop1"]["particle_initializer"]["nbr_part_per_cell"] = int{100};
    dict["ions"]["pop1"]["particle_initializer"]["charge"]            = -1.;
    dict["ions"]["pop1"]["particle_initializer"]["basis"]             = std::string{"cartesian"};

    dict["electromag"]["name"]             = std::string{"EM"};
    dict["electromag"]["electric"]["name"] = std::string{"E"};
    dict["electromag"]["magnetic"]["name"] = std::string{"B"};

    dict["electrons"]["pressure_closure"]["name"] = std::string{"isothermal"};
    dict["electrons"]["pressure_closure"]["Te"]   = 0.12;


    DimDict<dimension>::set(dict);

    return dict;
}


class BasicHierarchy
{
public:
    BasicHierarchy(
        int const ratio, short unsigned const dimension,
        SAMRAI::mesh::StandardTagAndInitStrategy* tagStrat,
        std::shared_ptr<SAMRAI::algs::TimeRefinementLevelStrategy> const& integratorStrat,
        std::string const& inputConfigFile)
        : inputDatabase_{SAMRAI::tbox::InputManager::getManager()->parseInputFile(inputConfigFile)}
        , patchHierarchyDatabase_{inputDatabase_->getDatabase("PatchHierarchy")}
        , dimension_{dimension}

        , gridGeometry_{std::make_shared<SAMRAI::geom::CartesianGridGeometry>(
              dimension_, "cartesian", inputDatabase_->getDatabase("CartesianGridGeometry"))}
        , hierarchy_{std::make_shared<SAMRAI::hier::PatchHierarchy>("PatchHierarchy", gridGeometry_,
                                                                    patchHierarchyDatabase_)}
        , loadBalancer_{std::make_shared<SAMRAI::mesh::TreeLoadBalancer>(
              dimension_, "LoadBalancer", inputDatabase_->getDatabase("LoadBalancer"))}
        , standardTag_{std::make_shared<SAMRAI::mesh::StandardTagAndInitialize>(
              "StandardTagAndInitialize", tagStrat,
              inputDatabase_->getDatabase("StandardTagAndInitialize"))}
        , clustering_{std::make_shared<SAMRAI::mesh::BergerRigoutsos>(
              dimension_, inputDatabase_->getDatabaseWithDefault(
                              "BergerRigoutsos", std::shared_ptr<SAMRAI::tbox::Database>()))}
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
    auto& hierarchy() { return hierarchy_; }

private:
    std::shared_ptr<SAMRAI::tbox::Database> inputDatabase_;
    std::shared_ptr<SAMRAI::tbox::Database> patchHierarchyDatabase_;

    SAMRAI::tbox::Dimension dimension_;

    std::shared_ptr<SAMRAI::geom::CartesianGridGeometry> gridGeometry_;

    std::shared_ptr<SAMRAI::hier::PatchHierarchy> hierarchy_;
    std::shared_ptr<SAMRAI::mesh::TreeLoadBalancer> loadBalancer_;

    std::shared_ptr<SAMRAI::mesh::StandardTagAndInitialize> standardTag_;
    std::shared_ptr<SAMRAI::mesh::BergerRigoutsos> clustering_;

public:
    std::shared_ptr<SAMRAI::mesh::GriddingAlgorithm> gridding;
    std::shared_ptr<SAMRAI::algs::TimeRefinementIntegrator> integrator;
};


template<std::uint8_t dim, std::size_t interpOrder, std::size_t nbRefinePart>
struct AfullHybridBasicHierarchy
{
    using Phare_core_Types   = PHARE::core::PHARE_Types<dim, interpOrder>;
    using Phare_amr_Types    = PHARE::amr::PHARE_Types<dim, interpOrder, nbRefinePart>;
    using Phare_solver_Types = PHARE::solver::PHARE_Types<dim, interpOrder, nbRefinePart>;
    using HybridModelT       = Phare_solver_Types::HybridModel_t;
    using MHDModelT          = Phare_solver_Types::MHDModel_t;
    using ResourcesManagerT  = typename HybridModelT::resources_manager_type;

    int const firstHybLevel{0};
    int const ratio{2};

    using HybridHybridT
        = HybridHybridMessengerStrategy<HybridModelT,
                                        typename Phare_solver_Types::RefinementParams>;

    AfullHybridBasicHierarchy(std::string const& inputConfigFile)
    {
        hybridModel->resourcesManager->registerResources(hybridModel->state);

        solver->registerResources(*hybridModel);

        tagStrat   = std::make_shared<TagStrategy<HybridModelT>>(hybridModel, solver, messenger);
        integrator = std::make_shared<TestIntegratorStrat>();
        basicHierarchy = std::make_shared<BasicHierarchy>(ratio, dim, tagStrat.get(), integrator,
                                                          inputConfigFile);
    }

    PHARE::initializer::PHAREDict dict{createDict<dim>()};
    SAMRAI::tbox::SAMRAI_MPI mpi{MPI_COMM_WORLD};

    std::shared_ptr<ResourcesManagerT> resourcesManagerHybrid{
        std::make_shared<ResourcesManagerT>()};

    std::shared_ptr<HybridModelT> hybridModel{
        std::make_shared<HybridModelT>(dict, resourcesManagerHybrid)};

    std::shared_ptr<HybridMessenger<HybridModelT>> messenger{
        std::make_shared<HybridMessenger<HybridModelT>>(
            std::make_unique<HybridHybridT>(resourcesManagerHybrid, firstHybLevel))};

    std::shared_ptr<SolverPPC<HybridModelT, SAMRAI_Types>> solver{
        std::make_shared<SolverPPC<HybridModelT, SAMRAI_Types>>(dict["simulation"]["algo"])};

    std::shared_ptr<TagStrategy<HybridModelT>> tagStrat;
    std::shared_ptr<TestIntegratorStrat> integrator;
    std::shared_ptr<BasicHierarchy> basicHierarchy;
};

#endif /*PHARE_TEST_AMR_HIERARCHY_FIXTURES_HPP*/
