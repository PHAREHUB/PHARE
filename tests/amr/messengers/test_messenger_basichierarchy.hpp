#ifndef PHARE_TEST_BASIC_HIERARCHY_HPP
#define PHARE_TEST_BASIC_HIERARCHY_HPP




#include <SAMRAI/algs/TimeRefinementIntegrator.h>
#include <SAMRAI/algs/TimeRefinementLevelStrategy.h>
#include <SAMRAI/geom/CartesianGridGeometry.h>
#include <SAMRAI/hier/PatchHierarchy.h>
#include <SAMRAI/hier/VariableContext.h>
#include <SAMRAI/mesh/ChopAndPackLoadBalancer.h>
#include <SAMRAI/mesh/GriddingAlgorithm.h>
#include <SAMRAI/mesh/StandardTagAndInitStrategy.h>
#include <SAMRAI/mesh/StandardTagAndInitialize.h>
#include <SAMRAI/mesh/TileClustering.h>

#include <SAMRAI/mesh/BergerRigoutsos.h>
#include <SAMRAI/mesh/TreeLoadBalancer.h>
#include <SAMRAI/tbox/InputManager.h>


#include <memory>
#include <string>


#include "input_config.h"


/**
 * @brief  In this class we will create a hierarchy of two levels
 * the refined level will be determined by the refined box defined
 * in the input files. The creation of the patch is managed directly
 * by SAMRAI ( thus the use of a GriddingAlgorithm ). To be able to
 * use a GriddingAlgorithm, we have to implement a TagStrategy.
 * The one implemented does only one thing : allocate data for each hybridQuantities
 * on each patch of a given level
 */

class BasicHierarchy
{
public:
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
    BasicHierarchy(
        int const ratio, short unsigned const dimension,
        SAMRAI::mesh::StandardTagAndInitStrategy* tagStrat,
        std::shared_ptr<SAMRAI::algs::TimeRefinementLevelStrategy> const& integratorStrat)
        : inputDatabase_{SAMRAI::tbox::InputManager::getManager()->parseInputFile(
              inputBase + "input/input_" + std::to_string(dimension) + "d_ratio_"
              + std::to_string(ratio) + ".txt")}
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

#endif
