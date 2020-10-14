#ifndef PHARE_TEST_BASIC_HIERARCHY_H
#define PHARE_TEST_BASIC_HIERARCHY_H


#include "amr/data/field/coarsening/field_coarsen_operator.h"
#include "amr/data/field/field_variable.h"
#include "test_tag_strategy.h"


#include <SAMRAI/geom/CartesianGridGeometry.h>
#include <SAMRAI/hier/PatchHierarchy.h>
#include <SAMRAI/hier/VariableContext.h>
#include <SAMRAI/mesh/ChopAndPackLoadBalancer.h>
#include <SAMRAI/mesh/GriddingAlgorithm.h>
#include <SAMRAI/mesh/StandardTagAndInitialize.h>
#include <SAMRAI/mesh/TileClustering.h>
#include <SAMRAI/tbox/InputManager.h>
#include <SAMRAI/xfer/CoarsenAlgorithm.h>
#include <SAMRAI/xfer/CoarsenSchedule.h>

#include <memory>
#include <string>


#include "input_config.h"

using namespace PHARE::core;
using namespace PHARE::amr;

/**
 * @brief  In this class we will create a hierarchy of two levels
 * the refined level will be determined by the refined box defined
 * in the input files. The creation of the patch is managed directly
 * by SAMRAI ( thus the use of a GriddingAlgorithm ). To be able to
 * use a GriddingAlgorithm, we have to implement a TagStrategy.
 * The one implemented does only one thing : allocate data for each hybridQuantities
 * on each patch of a given level
 */
template<typename GridLayoutT, typename FieldT,
         typename PhysicalQuantity = decltype(std::declval<FieldT>().physicalQuantity())>
class BasicHierarchy
{
public:
    static constexpr std::size_t dimension    = GridLayoutT::dimension;
    static constexpr std::size_t interp_order = GridLayoutT::interp_order;

    /**
     * @brief Construct the hierarchy consisting of two levels
     *
     * Only the ratio argument is required, we will use one input file per ratio.
     * Here we do multiple things, first we parse the inputFile to create the inputDatabase
     * that will be used by others SAMRAI's objects (GridGeometry, PatchHierarchy, LoadBalancer,
     * BoxGenerator, GriddingAlgorithm, StandardTagAndInitialize), they contain information for each
     * of them.
     *
     * We then store a view to the VariableDatabase singleton , that we will use to register all our
     * Variables. For each hybridQuantities we register a variable.
     *
     * We instantiate each objects needed by the GriddingAlgorithm, and then instantiate
     * a GriddingAlgorithm with the previous objects
     *
     */
    explicit BasicHierarchy(int ratio)
        : inputDatabase_{SAMRAI::tbox::InputManager::getManager()->parseInputFile(
            inputBase + "input/input_" + std::to_string(dimension) + "d_ratio_"
            + std::to_string(ratio) + ".txt")}
        , patchHierarchyDatabase_{inputDatabase_->getDatabase("PatchHierarchy")}
        , variableDatabase_{SAMRAI::hier::VariableDatabase::getDatabase()} //
        /* Ex Ey Ez */
        , ex_{std::make_shared<FieldVariable<GridLayoutT, FieldT>>("Ex", PhysicalQuantity::Ex)}
        , ey_{std::make_shared<FieldVariable<GridLayoutT, FieldT>>("Ey", PhysicalQuantity::Ey)}
        , ez_{std::make_shared<FieldVariable<GridLayoutT, FieldT>>("Ez", PhysicalQuantity::Ez)}
        // Bx By Bz
        , bx_{std::make_shared<FieldVariable<GridLayoutT, FieldT>>("Bx", PhysicalQuantity::Bx)}
        , by_{std::make_shared<FieldVariable<GridLayoutT, FieldT>>("By", PhysicalQuantity::By)}
        , bz_{std::make_shared<FieldVariable<GridLayoutT, FieldT>>("Bz", PhysicalQuantity::Bz)}
        // Jx Jy Jz
        , jx_{std::make_shared<FieldVariable<GridLayoutT, FieldT>>("Jx", PhysicalQuantity::Jx)}
        , jy_{std::make_shared<FieldVariable<GridLayoutT, FieldT>>("Jy", PhysicalQuantity::Jy)}
        , jz_{std::make_shared<FieldVariable<GridLayoutT, FieldT>>("Jz", PhysicalQuantity::Jz)}
        // Vx Vy Vz
        , vx_{std::make_shared<FieldVariable<GridLayoutT, FieldT>>("Vx", PhysicalQuantity::Vx)}
        , vy_{std::make_shared<FieldVariable<GridLayoutT, FieldT>>("Vy", PhysicalQuantity::Vy)}
        , vz_{std::make_shared<FieldVariable<GridLayoutT, FieldT>>("Vz", PhysicalQuantity::Vz)} //
                                                                                                //
        /* rho */
        , rho_{std::make_shared<FieldVariable<GridLayoutT, FieldT>>("Rho", PhysicalQuantity::rho)}

        /* P */
        , p_{std::make_shared<FieldVariable<GridLayoutT, FieldT>>("P", PhysicalQuantity::P)}

        , context_{variableDatabase_->getContext("context")}
        , variablesIds_{getVariablesIds_()}

        , gridGeometry_{std::make_shared<SAMRAI::geom::CartesianGridGeometry>(
              dimension_, "cartesian", inputDatabase_->getDatabase("CartesianGridGeometry"))}
        , hierarchy_{std::make_shared<SAMRAI::hier::PatchHierarchy>("PatchHierarchy", gridGeometry_,
                                                                    patchHierarchyDatabase_)}
        , loadBalancer_{std::make_shared<SAMRAI::mesh::ChopAndPackLoadBalancer>(
              dimension_, "ChopAndPackLoadBalancer",
              inputDatabase_->getDatabase("ChopAndPackLoadBalancer"))}
        , tagStrategy_{std::make_shared<TagStrategy>(variablesIds_)}
        , standardTag_{std::make_shared<SAMRAI::mesh::StandardTagAndInitialize>(
              "StandardTagAndInitialize", tagStrategy_.get(),
              inputDatabase_->getDatabase("StandardTagAndInitialize"))}

        , clustering_{std::make_shared<SAMRAI::mesh::TileClustering>(
              dimension_, inputDatabase_->getDatabase("TileClustering"))}
        , gridding_{std::make_shared<SAMRAI::mesh::GriddingAlgorithm>(
              hierarchy_, "GriddingAlgorithm", inputDatabase_->getDatabase("GriddingAlgorithm"),
              standardTag_, clustering_, loadBalancer_)}
    {
        // First step is to make the CoarsestLevel, and then to refine
        // as much as needed
        gridding_->makeCoarsestLevel(0.0);
        // we refine for tag 0
        // For that we try to refine until we get no more refine boxes or that
        // we have reach the maximum level
        for (int iLevel = 0, hierarchyMaxLevelAllowed = hierarchy_->getMaxNumberOfLevels();
             iLevel < hierarchyMaxLevelAllowed - 1; ++iLevel)
        {
            int const cycle   = 0;
            double const time = 0.0;

            SAMRAI::hier::BoxContainer boxes{};
            auto reset = standardTag_->getUserSuppliedRefineBoxes(boxes, iLevel, cycle, time);
            NULL_USE(reset);
            if (!boxes.empty())
            {
                gridding_->makeFinerLevel(0, true, cycle, time, 0.0);
            }
            else
            {
                break;
            }
        }

        // We have our hierarchy setup, now is time to register the coarsenOperator
        // that we will use

        auto fieldVariableTypeName = typeid(FieldVariable<GridLayoutT, FieldT>).name();

        auto coarsenOperator = std::make_shared<FieldCoarsenOperator<GridLayoutT, FieldT>>();

        gridGeometry_->addCoarsenOperator(fieldVariableTypeName, coarsenOperator);

        // register variable for coarse algorithm
        for (auto const& nameToIds : variablesIds_)
        {
            coarseAlg_.registerCoarsen(nameToIds.second, nameToIds.second, coarsenOperator);
        }

        // create schedule
        coarseSchedule_
            = coarseAlg_.createSchedule(hierarchy_->getPatchLevel(0), hierarchy_->getPatchLevel(1));
    }



    SAMRAI::hier::PatchHierarchy& getHierarchy() { return *hierarchy_; }


    std::map<std::string, int> const& getVariables() const { return variablesIds_; }

    void coarsify() { coarseSchedule_->coarsenData(); }

    void TearDown()
    {
        for (auto const& namesToId : variablesIds_)
        {
            variableDatabase_->removeVariable(namesToId.first);
        }
    }

private:
    std::map<std::string, int> getVariablesIds_()
    {
        std::map<std::string, int> variablesIds;

        SAMRAI::hier::IntVector ghostWidth{dimension_, 5};

        variablesIds.try_emplace(
            "Ex", variableDatabase_->registerVariableAndContext(ex_, context_, ghostWidth));
        variablesIds.try_emplace(
            "Ey", variableDatabase_->registerVariableAndContext(ey_, context_, ghostWidth));
        variablesIds.try_emplace(
            "Ez", variableDatabase_->registerVariableAndContext(ez_, context_, ghostWidth));
        variablesIds.try_emplace(
            "Bx", variableDatabase_->registerVariableAndContext(bx_, context_, ghostWidth));
        variablesIds.try_emplace(
            "By", variableDatabase_->registerVariableAndContext(by_, context_, ghostWidth));
        variablesIds.try_emplace(
            "Bz", variableDatabase_->registerVariableAndContext(bz_, context_, ghostWidth));
        variablesIds.try_emplace(
            "Jx", variableDatabase_->registerVariableAndContext(jx_, context_, ghostWidth));
        variablesIds.try_emplace(
            "Jy", variableDatabase_->registerVariableAndContext(jy_, context_, ghostWidth));
        variablesIds.try_emplace(
            "Jz", variableDatabase_->registerVariableAndContext(jz_, context_, ghostWidth));
        variablesIds.try_emplace(
            "Rho", variableDatabase_->registerVariableAndContext(rho_, context_, ghostWidth));
        variablesIds.try_emplace(
            "P", variableDatabase_->registerVariableAndContext(p_, context_, ghostWidth));

        return variablesIds;
    }

    std::shared_ptr<SAMRAI::tbox::Database> inputDatabase_;
    std::shared_ptr<SAMRAI::tbox::Database> patchHierarchyDatabase_;

    SAMRAI::tbox::Dimension dimension_{dimension};

    SAMRAI::hier::VariableDatabase* variableDatabase_;

    // Each FieldVariable for each hybridQuantities
    std::shared_ptr<FieldVariable<GridLayoutT, FieldT>> ex_;
    std::shared_ptr<FieldVariable<GridLayoutT, FieldT>> ey_;
    std::shared_ptr<FieldVariable<GridLayoutT, FieldT>> ez_;


    std::shared_ptr<FieldVariable<GridLayoutT, FieldT>> bx_;
    std::shared_ptr<FieldVariable<GridLayoutT, FieldT>> by_;
    std::shared_ptr<FieldVariable<GridLayoutT, FieldT>> bz_;


    std::shared_ptr<FieldVariable<GridLayoutT, FieldT>> jx_;
    std::shared_ptr<FieldVariable<GridLayoutT, FieldT>> jy_;
    std::shared_ptr<FieldVariable<GridLayoutT, FieldT>> jz_;


    std::shared_ptr<FieldVariable<GridLayoutT, FieldT>> vx_;
    std::shared_ptr<FieldVariable<GridLayoutT, FieldT>> vy_;
    std::shared_ptr<FieldVariable<GridLayoutT, FieldT>> vz_;

    std::shared_ptr<FieldVariable<GridLayoutT, FieldT>> rho_;
    std::shared_ptr<FieldVariable<GridLayoutT, FieldT>> p_;


    std::shared_ptr<SAMRAI::hier::VariableContext> context_;

    std::map<std::string, int> variablesIds_;

    std::shared_ptr<SAMRAI::geom::CartesianGridGeometry> gridGeometry_;

    std::shared_ptr<SAMRAI::hier::PatchHierarchy> hierarchy_;
    std::shared_ptr<SAMRAI::mesh::ChopAndPackLoadBalancer> loadBalancer_;
    std::shared_ptr<TagStrategy> tagStrategy_;
    std::shared_ptr<SAMRAI::mesh::StandardTagAndInitialize> standardTag_;
    std::shared_ptr<SAMRAI::mesh::TileClustering> clustering_;

    std::shared_ptr<SAMRAI::mesh::GriddingAlgorithm> gridding_;

    SAMRAI::xfer::CoarsenAlgorithm coarseAlg_{dimension_};

    std::shared_ptr<SAMRAI::xfer::CoarsenSchedule> coarseSchedule_;
};

#endif
