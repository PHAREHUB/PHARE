#ifndef PHARE_SIMULATOR_H
#define PHARE_SIMULATOR_H

#include "python_data_provider.h"
// intended blank HAVE_SYS_TIMES_His defined by samrai
#include "amr/types/amr_types.h"
#include "data/electromag/electromag.h"
#include "data/grid/gridlayout.h"
#include "data/grid/gridlayoutimplyee.h"
#include "data/ions/ion_population/ion_population.h"
#include "data/ions/ions.h"
#include "data/ions/particle_initializers/maxwellian_particle_initializer.h"
#include "data/ndarray/ndarray_vector.h"
#include "data/particles/particle_array.h"
#include "data/vecfield/vecfield.h"
#include "data_provider.h"
#include "messengers/messenger_factory.h"
#include "models/physical_state.h"
#include "multiphysics_integrator.h"
#include "physical_models/hybrid_model.h"
#include "physical_models/mhd_model.h"
#include "physical_models/physical_model.h"
#include "solvers/solver.h"
#include "solvers/solver_mhd.h"
#include "solvers/solver_ppc.h"
#include "types/amr_types.h"
#include "utilities/algorithm.h"
#include "utilities/meta/meta_utilities.h"

#include <SAMRAI/algs/TimeRefinementIntegrator.h>
#include <SAMRAI/geom/CartesianGridGeometry.h>
#include <SAMRAI/hier/Box.h>
#include <SAMRAI/hier/BoxContainer.h>
#include <SAMRAI/hier/IntVector.h>
#include <SAMRAI/hier/PatchHierarchy.h>
#include <SAMRAI/mesh/BergerRigoutsos.h>
#include <SAMRAI/mesh/GriddingAlgorithm.h>
#include <SAMRAI/mesh/StandardTagAndInitialize.h>
#include <SAMRAI/mesh/TreeLoadBalancer.h>
#include <SAMRAI/tbox/Database.h>
#include <SAMRAI/tbox/DatabaseBox.h>
#include <SAMRAI/tbox/InputManager.h>
#include <SAMRAI/tbox/MemoryDatabase.h>


#include <algorithm>
#include <iterator>
#include <memory>
#include <sstream>




template<std::size_t dimension>
void getDomainCoords(PHARE::initializer::PHAREDict& grid, float lower[dimension],
                     float upper[dimension])
{
    static_assert(dimension > 0 and dimension <= 3, "invalid dimension should be >0 and <=3");

    if constexpr (dimension >= 1)
    {
        auto nx  = grid["nbr_cells"]["x"].template to<int>();
        auto dx  = grid["meshsize"]["x"].template to<double>();
        lower[0] = static_cast<float>(grid["origin"]["x"].template to<double>());
        upper[0] = static_cast<float>(lower[0] + nx * dx);

        if constexpr (dimension >= 2)
        {
            auto ny  = grid["nbr_cells"]["y"].template to<int>();
            auto dy  = grid["meshsize"]["y"].template to<double>();
            lower[1] = static_cast<float>(grid["origin"]["y"].template to<double>());
            upper[1] = static_cast<float>(lower[1] + ny * dy);

            if constexpr (dimension == 3)
            {
                auto nz  = grid["nbr_cells"]["z"].template to<int>();
                auto dz  = grid["meshsize"]["z"].template to<double>();
                lower[2] = static_cast<float>(grid["origin"]["z"].template to<double>());
                upper[2] = static_cast<float>(lower[2] + nz * dz);
            }
        }
    }
}



template<std::size_t dimension>
auto griddingAlgorithmDatabase(PHARE::initializer::PHAREDict& grid)
{
    static_assert(dimension > 0 and dimension <= 3, "invalid dimension should be >0 and <=3");
    auto samraiDim = SAMRAI::tbox::Dimension{dimension};
    auto db        = std::make_shared<SAMRAI::tbox::MemoryDatabase>("griddingAlgoDB");

    int lowerCell[dimension];
    int upperCell[dimension];
    float lowerCoord[dimension];
    float upperCoord[dimension];
    int periodicity[dimension];
    std::vector<SAMRAI::tbox::DatabaseBox> dbBoxes;

    if constexpr (dimension == 1)
    {
        lowerCell[0]   = 0;
        upperCell[0]   = grid["nbr_cells"]["x"].template to<int>() - 1;
        periodicity[0] = 1; // 1==periodic, hardedcoded for all dims for now.
    }


    if constexpr (dimension == 2)
    {
        lowerCell[0]   = 0;
        upperCell[0]   = grid["nbr_cells"]["x"].template to<int>();
        lowerCell[1]   = 0;
        upperCell[1]   = grid["nbr_cells"]["y"].template to<int>();
        periodicity[0] = 1;
        periodicity[1] = 1;
    }
    if constexpr (dimension == 3)
    {
        lowerCell[0]   = 0;
        upperCell[0]   = grid["nbr_cells"]["x"].template to<int>();
        lowerCell[1]   = 0;
        upperCell[1]   = grid["nbr_cells"]["y"].template to<int>();
        lowerCell[2]   = 0;
        upperCell[2]   = grid["nbr_cells"]["z"].template to<int>();
        periodicity[0] = 1;
        periodicity[1] = 1;
        periodicity[2] = 1;
    }


    dbBoxes.push_back(SAMRAI::tbox::DatabaseBox(samraiDim, lowerCell, upperCell));
    db->putDatabaseBoxVector("domain_boxes", dbBoxes);
    getDomainCoords<dimension>(grid, lowerCoord, upperCoord);
    db->putFloatArray("x_lo", lowerCoord, dimension);
    db->putFloatArray("x_up", upperCoord, dimension);
    db->putIntegerArray("periodic_dimension", periodicity, dimension);
    return db;
}



/*
// Required input: maximum number of levels in patch hierarchy
max_levels = 4
// Required input: vector ratio between each finer level and next coarser
ratio_to_coarser {
   level_1 = 2, 2, 2
   level_2 = 2, 2, 2
   level_3 = 4, 4, 4
}
// Optional input: int vector for largest patch size on each level.
largest_patch_size {
   level_0 = 40, 40, 40
   level_1 = 30, 30, 30
   // all finer levels will use same values as level_1...
}
// Optional input: int vector for smallest patch size on each level.
smallest_patch_size {
   level_0 = 16, 16, 16
   // all finer levels will use same values as level_0...
}
// Optional input:  buffer of one cell used on each level
proper_nesting_buffer = 1
 */
template<std::size_t dimension>
auto patchHierarchyDatabase(PHARE::initializer::PHAREDict& amr)
{
    auto maxLevelNumber = amr["max_nbr_levels"].template to<int>();

    auto hierDB = std::make_shared<SAMRAI::tbox::MemoryDatabase>("HierarchyDB");

    hierDB->putInteger("max_levels", maxLevelNumber);

    auto smallestPatchSizeDB = hierDB->putDatabase("smallest_patch_size");
    auto ratioToCoarserDB    = hierDB->putDatabase("ratio_to_coarser");

    constexpr int smallestPatchSizeValue = 10;
    constexpr int ratio                  = 2;

    for (auto iLevel = 0; iLevel < maxLevelNumber; ++iLevel)
    {
        if constexpr (dimension == 1)
        {
            if (iLevel > 0)
            {
                int ratioToCoarser[dimension] = {ratio};
                ratioToCoarserDB->putIntegerArray("level_" + std::to_string(iLevel), ratioToCoarser,
                                                  dimension);
            }
            int smallestPatchSize[dimension] = {smallestPatchSizeValue};
            smallestPatchSizeDB->putIntegerArray("level_" + std::to_string(iLevel),
                                                 smallestPatchSize, dimension);
        }
        else if constexpr (dimension == 2)
        {
            if (iLevel > 0)
            {
                int ratioToCoarser[dimension] = {ratio, ratio};
                ratioToCoarserDB->putIntegerArray("level_" + std::to_string(iLevel), ratioToCoarser,
                                                  dimension);
            }
            int smallestPatchSize[dimension] = {smallestPatchSizeValue, smallestPatchSizeValue};
            smallestPatchSizeDB->putIntegerArray("smallest_patch_size", smallestPatchSize,
                                                 dimension);
        }
        if constexpr (dimension == 3)
        {
            if (iLevel > 0)
            {
                int ratioToCoarser[dimension] = {2, 2, 2};
                ratioToCoarserDB->putIntegerArray("level_" + std::to_string(iLevel), ratioToCoarser,
                                                  dimension);
            }
            int smallestPatchSize[dimension] = {10, 10, 10};
            smallestPatchSizeDB->putIntegerArray("smallest_patch_size", smallestPatchSize,
                                                 dimension);
        }
    }
    return hierDB;
}



template<std::size_t dimension>
std::shared_ptr<SAMRAI::tbox::MemoryDatabase>
getUserRefinementBoxesDatabase(PHARE::initializer::PHAREDict& amr)
{
    try
    {
        [[maybe_unused]] auto& refDict = amr["refinement_boxes"];
    }
    catch (...)
    {
        return nullptr;
    }
    auto maxLevelNumber = amr["max_nbr_levels"].template to<int>();
    auto& refDict       = amr["refinement_boxes"];


    std::shared_ptr<SAMRAI::tbox::MemoryDatabase> refinementBoxesDatabase
        = std::make_shared<SAMRAI::tbox::MemoryDatabase>("StandardTagAndInitialize");

    // user refinement boxes are always defined at t=0
    auto at0db = refinementBoxesDatabase->putDatabase("at_0");
    at0db->putInteger("cycle", 0);
    auto tag0db = at0db->putDatabase("tag_0");
    tag0db->putString("tagging_method", "REFINE_BOXES");


    for (int levelNumber = 0; levelNumber < maxLevelNumber; ++levelNumber)
    {
        // not all levels are necessarily specified for refinement
        // cppdict will throw when trying to access key L{i} with i = levelNumber
        try
        {
            auto& levelDict = refDict["L" + std::to_string(levelNumber)];
            auto samraiDim  = SAMRAI::tbox::Dimension{dimension};
            auto nbrBoxes   = levelDict["nbr_boxes"].template to<int>();
            auto levelDB    = tag0db->putDatabase("level_" + std::to_string(levelNumber));

            std::vector<SAMRAI::tbox::DatabaseBox> dbBoxes;
            for (int iBox = 0; iBox < nbrBoxes; ++iBox)
            {
                int lower[dimension];
                int upper[dimension];
                auto& boxDict = levelDict["B" + std::to_string(iBox)];

                if constexpr (dimension == 1)
                {
                    lower[0] = boxDict["lower"]["x"].template to<int>();
                    upper[0] = boxDict["upper"]["x"].template to<int>();
                }

                if constexpr (dimension == 2)
                {
                    lower[0] = boxDict["lower"]["x"].template to<int>();
                    lower[1] = boxDict["lower"]["y"].template to<int>();

                    upper[0] = boxDict["upper"]["x"].template to<int>();
                    upper[1] = boxDict["upper"]["y"].template to<int>();
                }

                if constexpr (dimension == 3)
                {
                    lower[0] = boxDict["lower"]["x"].template to<int>();
                    lower[1] = boxDict["lower"]["y"].template to<int>();
                    lower[2] = boxDict["lower"]["z"].template to<int>();

                    upper[0] = boxDict["upper"]["x"].template to<int>();
                    upper[1] = boxDict["upper"]["y"].template to<int>();
                    upper[2] = boxDict["upper"]["z"].template to<int>();
                }


                dbBoxes.push_back(SAMRAI::tbox::DatabaseBox(samraiDim, lower, upper));
            }
            levelDB->putDatabaseBoxVector("boxes", dbBoxes);
        }
        catch (...)
        {
        }
    } // end loop on levels
    return refinementBoxesDatabase;
}




namespace PHARE
{
template<std::size_t dimension, std::size_t interp_order>
struct PHARE_Types
{
    using Array_t
        = decltype(PHARE::core::makeNdArray<dimension>(std::array<std::uint32_t, dimension>{}));
    using VecField_t      = PHARE::core::VecField<Array_t, PHARE::core::HybridQuantity>;
    using Field_t         = PHARE::core::Field<Array_t, PHARE::core::HybridQuantity::Scalar>;
    using Electromag_t    = PHARE::core::Electromag<VecField_t>;
    using YeeLayout_t     = PHARE::core::GridLayoutImplYee<dimension, interp_order>;
    using GridLayout_t    = PHARE::core::GridLayout<YeeLayout_t>;
    using ParticleArray_t = PHARE::core::ParticleArray<dimension>;
    using MaxwellianParticleInitializer_t
        = PHARE::core::MaxwellianParticleInitializer<ParticleArray_t, GridLayout_t>;
    using IonPopulation_t = PHARE::core::IonPopulation<ParticleArray_t, VecField_t, GridLayout_t>;
    using Ions_t          = PHARE::core::Ions<IonPopulation_t, GridLayout_t>;
    using HybridModel_t
        = PHARE::solver::HybridModel<GridLayout_t, Electromag_t, Ions_t, PHARE::amr::SAMRAI_Types>;
    using MHDModel_t  = PHARE::solver::MHDModel<GridLayout_t, VecField_t, PHARE::amr::SAMRAI_Types>;
    using SolverPPC_t = PHARE::solver::SolverPPC<HybridModel_t, PHARE::amr::SAMRAI_Types>;
    using SolverMHD_t = PHARE::solver::SolverMHD<MHDModel_t, PHARE::amr::SAMRAI_Types>;
};



class ISimulator
{
public:
    virtual void initialize()    = 0;
    virtual double startTime()   = 0;
    virtual double endTime()     = 0;
    virtual double timeStep()    = 0;
    virtual void advance()       = 0;
    virtual std::string to_str() = 0;
    virtual ~ISimulator() {}
};




template<std::size_t dimension, std::size_t interp_order>
class Simulator : public ISimulator
{
private:
    using SAMRAITypes = PHARE::amr::SAMRAI_Types;
    using PHARETypes  = PHARE_Types<dimension, interp_order>;

    using IPhysicalModel = PHARE::solver::IPhysicalModel<SAMRAITypes>;
    using HybridModel    = typename PHARETypes::HybridModel_t;
    using MHDModel       = typename PHARETypes::MHDModel_t;

    using SolverMHD = typename PHARETypes::SolverMHD_t;
    using SolverPPC = typename PHARETypes::SolverPPC_t;

    using MessengerFactory = PHARE::amr::MessengerFactory<MHDModel, HybridModel, IPhysicalModel>;

    using MultiPhysicsIntegrator
        = PHARE::solver::MultiPhysicsIntegrator<MessengerFactory, SAMRAITypes>;

public:
    Simulator(PHARE::initializer::PHAREDict dict)
        : modelNames_{"HybridModel"}
        , descriptors_{PHARE::amr::makeDescriptors(modelNames_)}
        , messengerFactory_{descriptors_}
        , maxLevelNumber_{dict["simulation"]["AMR"]["max_nbr_levels"].template to<int>()}
        , dt_{dict["simulation"]["time_step"].template to<double>()}
        , timeStepNbr_{dict["simulation"]["time_step_nbr"].template to<int>()}
        , finalTime_{dt_ * timeStepNbr_}
        , multiphysInteg_{std::make_shared<MultiPhysicsIntegrator>(
              dict["simulation"]["AMR"]["max_nbr_levels"].template to<int>())}
    {
        if (find_model("HybridModel"))
        {
            hybridModel_ = std::make_shared<HybridModel>(
                dict["simulation"],
                std::make_shared<typename HybridModel::resources_manager_type>());


            hybridModel_->resourcesManager->registerResources(hybridModel_->state);

            // we register the hybrid model for all possible levels in the hierarchy
            // since for now it is the only model available
            // same for the solver
            multiphysInteg_->registerModel(0, maxLevelNumber_ - 1, hybridModel_);
            multiphysInteg_->registerAndInitSolver(0, maxLevelNumber_ - 1,
                                                   std::make_unique<SolverPPC>());
            multiphysInteg_->registerAndSetupMessengers(messengerFactory_);

            auto samraiDim    = SAMRAI::tbox::Dimension{dimension};
            auto gridGeometry = std::make_shared<SAMRAI::geom::CartesianGridGeometry>(
                samraiDim, "CartesianGridGeom",
                griddingAlgorithmDatabase<dimension>(dict["simulation"]["grid"]));


            auto hierDB = patchHierarchyDatabase<dimension>(dict["simulation"]["AMR"]);
            hierarchy_  = std::make_shared<SAMRAI::hier::PatchHierarchy>("PHARE_hierarchy",
                                                                        gridGeometry, hierDB);


            auto loadBalancer = std::make_shared<SAMRAI::mesh::TreeLoadBalancer>(
                SAMRAI::tbox::Dimension{dimension}, "LoadBalancer");

            auto refineDB    = getUserRefinementBoxesDatabase<dimension>(dict["simulation"]["AMR"]);
            auto standardTag = std::make_shared<SAMRAI::mesh::StandardTagAndInitialize>(
                "StandardTagAndInitialize", multiphysInteg_.get(), refineDB);


            auto clustering = std::make_shared<SAMRAI::mesh::BergerRigoutsos>(
                SAMRAI::tbox::Dimension{dimension});

            auto gridding = std::make_shared<SAMRAI::mesh::GriddingAlgorithm>(
                hierarchy_, "GriddingAlgorithm", std::shared_ptr<SAMRAI::tbox::Database>{},
                standardTag, clustering, loadBalancer);

            std::shared_ptr<SAMRAI::tbox::Database> db
                = std::make_shared<SAMRAI::tbox::MemoryDatabase>("TRIdb");

            db->putDouble("start_time", 0.);
            db->putDouble("end_time", dt_ * timeStepNbr_);
            db->putInteger("max_integrator_steps", 1000000);
            db->putInteger("regrid_interval", 1);


            timeRefIntegrator_ = std::make_shared<SAMRAI::algs::TimeRefinementIntegrator>(
                "TimeRefinementIntegrator", db, hierarchy_, multiphysInteg_, gridding);
        }
        else
            throw std::runtime_error("unsupported model");
    }

    virtual std::string to_str() override
    {
        std::stringstream ss;
        ss << "PHARE SIMULATOR\n";
        ss << "------------------------------------\n";
        ss << "interpolation order  : " << interp_order << "\n";
        ss << "dimension            : " << dimension << "\n";
        ss << "time step            : " << dt_ << "\n";
        ss << "number of time steps : " << timeStepNbr_ << "\n";
        ss << core::to_str(hybridModel_->state);
        return ss.str();
    }


    virtual void initialize() override
    {
        std::cerr << "initializing hierarchy...";
        timeRefIntegrator_->initializeHierarchy();
        std::cerr << "done !\n";
    }
    virtual double startTime() override { return 0.; }
    virtual double endTime() override { return 100.; }
    virtual double timeStep() override { return 0.01; }
    virtual void advance() override {}


private:
    auto find_model(std::string name)
    {
        return std::find(std::begin(modelNames_), std::end(modelNames_), name)
               != std::end(modelNames_);
    }




    std::vector<std::string> modelNames_;
    std::vector<PHARE::amr::MessengerDescriptor> descriptors_;
    MessengerFactory messengerFactory_;

    float x_lo_[dimension];
    float x_up_[dimension];

    int maxLevelNumber_;
    double dt_;
    int timeStepNbr_;
    double finalTime_;


    // physical models that can be used
    std::shared_ptr<HybridModel> hybridModel_;
    std::shared_ptr<MHDModel> mhdModel_;


    std::shared_ptr<MultiPhysicsIntegrator> multiphysInteg_;
    std::shared_ptr<SAMRAI::algs::TimeRefinementIntegrator> timeRefIntegrator_;
    std::shared_ptr<SAMRAI::hier::PatchHierarchy> hierarchy_;
};


struct makeSimulator
{
    template<typename Dimension, typename InterpOrder>
    std::unique_ptr<ISimulator> operator()(std::size_t userDim, std::size_t userInterpOrder,
                                           Dimension dimension, InterpOrder interp_order)
    {
        if (userDim == dimension() and userInterpOrder == interp_order())
        {
            std::size_t constexpr d  = dimension();
            std::size_t constexpr io = interp_order();

            PHARE::initializer::PHAREDict& theDict
                = PHARE::initializer::PHAREDictHandler::INSTANCE().dict();
            return std::make_unique<Simulator<d, io>>(theDict);
        }
        else
        {
            return nullptr;
        }
    }
};


std::unique_ptr<PHARE::ISimulator> getSimulator()
{
    PHARE::initializer::PHAREDict& theDict
        = PHARE::initializer::PHAREDictHandler::INSTANCE().dict();
    auto dim         = theDict["simulation"]["dimension"].template to<int>();
    auto interpOrder = theDict["simulation"]["interp_order"].template to<int>();
    return core::makeAtRuntime<makeSimulator>(dim, interpOrder, makeSimulator{});
}




} // namespace PHARE
#endif
