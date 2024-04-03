#ifndef INTEGRATOR_HPP
#define INTEGRATOR_HPP

#include "core/def/phare_mpi.hpp"

#include <SAMRAI/algs/TimeRefinementIntegrator.h>
#include <SAMRAI/geom/CartesianGridGeometry.h>
#include <SAMRAI/hier/Box.h>
#include <SAMRAI/hier/BoxContainer.h>
#include <SAMRAI/hier/IntVector.h>
#include <SAMRAI/hier/PatchHierarchy.h>
#include <SAMRAI/mesh/TileClustering.h>
#include <SAMRAI/mesh/BergerRigoutsos.h>
#include <SAMRAI/mesh/GriddingAlgorithm.h>
#include <SAMRAI/mesh/StandardTagAndInitialize.h>
#include <SAMRAI/mesh/TreeLoadBalancer.h>
#include <SAMRAI/tbox/Database.h>
#include <SAMRAI/tbox/DatabaseBox.h>
#include <SAMRAI/tbox/InputManager.h>
#include <SAMRAI/tbox/MemoryDatabase.h>



#include "initializer/data_provider.hpp"


namespace PHARE::amr
{
template<std::size_t _dimension>
class Integrator
{
public:
    static constexpr std::size_t dimension = _dimension;

    double advance(double dt) { return timeRefIntegrator_->advanceHierarchy(dt); }

    void initialize() { timeRefIntegrator_->initializeHierarchy(); }


    Integrator(PHARE::initializer::PHAREDict const& dict,
               std::shared_ptr<PHARE::amr::Hierarchy> hierarchy,
               std::shared_ptr<SAMRAI::algs::TimeRefinementLevelStrategy> timeRefLevelStrategy,
               std::shared_ptr<SAMRAI::mesh::StandardTagAndInitStrategy> tagAndInitStrategy,
               double startTime, double endTime);

private:
    std::shared_ptr<SAMRAI::algs::TimeRefinementIntegrator> timeRefIntegrator_;
};




//-----------------------------------------------------------------------------
//                           Definitions
//-----------------------------------------------------------------------------


template<std::size_t dimension>
std::shared_ptr<SAMRAI::tbox::MemoryDatabase>
getUserRefinementBoxesDatabase(PHARE::initializer::PHAREDict const& amr, bool isFromRestart);




template<std::size_t _dimension>
Integrator<_dimension>::Integrator(
    PHARE::initializer::PHAREDict const& dict, std::shared_ptr<PHARE::amr::Hierarchy> hierarchy,
    std::shared_ptr<SAMRAI::algs::TimeRefinementLevelStrategy> timeRefLevelStrategy,
    std::shared_ptr<SAMRAI::mesh::StandardTagAndInitStrategy> tagAndInitStrategy, double startTime,
    double endTime)
{
    auto loadBalancer = std::make_shared<SAMRAI::mesh::TreeLoadBalancer>(
        SAMRAI::tbox::Dimension{dimension}, "LoadBalancer");

    auto refineDB    = getUserRefinementBoxesDatabase<dimension>(dict["simulation"]["AMR"],
                                                              hierarchy->isFromRestart());
    auto standardTag = std::make_shared<SAMRAI::mesh::StandardTagAndInitialize>(
        "StandardTagAndInitialize", tagAndInitStrategy.get(), refineDB);


    auto clustering = [&]() -> std::shared_ptr<SAMRAI::mesh::BoxGeneratorStrategy> {
        if (!dict["simulation"]["AMR"].contains("clustering"))
            throw std::runtime_error(std::string{"clustering type not specificed"});

        auto clustering_type = dict["simulation"]["AMR"]["clustering"].template to<std::string>();

        if (clustering_type == "berger")
        {
            std::shared_ptr<SAMRAI::tbox::Database> bergerDB
                = std::make_shared<SAMRAI::tbox::MemoryDatabase>("Bergerdb");
            bergerDB->putIntegerVector("max_box_size", std::vector<int>(dimension, 10));
            return std::make_shared<SAMRAI::mesh::BergerRigoutsos>(
                SAMRAI::tbox::Dimension{dimension}, bergerDB);
        }

        if (clustering_type == "tile")
            return std::make_shared<SAMRAI::mesh::TileClustering>(
                SAMRAI::tbox::Dimension{dimension});

        throw std::runtime_error(std::string{"Unknown clustering type "} + clustering_type);
    }();

    auto gridding = std::make_shared<SAMRAI::mesh::GriddingAlgorithm>(
        hierarchy, "GriddingAlgorithm", std::shared_ptr<SAMRAI::tbox::Database>{}, standardTag,
        clustering, loadBalancer);

    std::shared_ptr<SAMRAI::tbox::Database> db
        = std::make_shared<SAMRAI::tbox::MemoryDatabase>("TRIdb");


    db->putDouble("start_time", startTime);
    db->putDouble("end_time", endTime);
    db->putInteger("max_integrator_steps", 1000000);
    db->putIntegerVector(
        "tag_buffer", std::vector<int>(hierarchy->getMaxNumberOfLevels(),
                                       dict["simulation"]["AMR"]["tag_buffer"].template to<int>()));


    timeRefIntegrator_ = std::make_shared<SAMRAI::algs::TimeRefinementIntegrator>(
        "TimeRefinementIntegrator", db, hierarchy, timeRefLevelStrategy, gridding);
}




template<std::size_t dimension>
std::shared_ptr<SAMRAI::tbox::MemoryDatabase>
getUserRefinementBoxesDatabase(PHARE::initializer::PHAREDict const& amr, bool isFromRestart)
{
    auto const& refinement = amr["refinement"];
    auto maxLevelNumber    = amr["max_nbr_levels"].template to<int>();

    auto standardTagInitDB
        = std::make_shared<SAMRAI::tbox::MemoryDatabase>("StandardTagAndInitialize");

    if (refinement.contains("boxes"))
    {
        auto& refDict = refinement["boxes"];

        // user refinement boxes are always defined at t=0
        // auto at0db = refinementBoxesDatabase->putDatabase("at_0");
        // at0db->putInteger("cycle", 0);
        // auto tag0db = at0db->putDatabase("tag_0");
        std::cout << "tagging method is set to REFINE_BOXES\n";

        if (isFromRestart)
            standardTagInitDB->putString("tagging_method", "GRADIENT_DETECTOR");
        else // restarts are special
            standardTagInitDB->putString("tagging_method", "REFINE_BOXES");

        for (int levelNumber = 0; levelNumber < maxLevelNumber; ++levelNumber)
        {
            // not all levels are necessarily specified for refinement
            // cppdict will throw when trying to access key L{i} with i = levelNumber
            std::string levelString{"L" + std::to_string(levelNumber)};
            if (refDict.contains(levelString))
            {
                auto& levelDict = refDict[levelString];
                auto samraiDim  = SAMRAI::tbox::Dimension{dimension};
                auto nbrBoxes   = levelDict["nbr_boxes"].template to<int>();
                auto levelDB
                    = standardTagInitDB->putDatabase("level_" + std::to_string(levelNumber));

                std::vector<SAMRAI::tbox::DatabaseBox> dbBoxes;
                for (int iBox = 0; iBox < nbrBoxes; ++iBox)
                {
                    int lower[dimension];
                    int upper[dimension];
                    auto& boxDict = levelDict["B" + std::to_string(iBox)];

                    lower[0] = boxDict["lower"]["x"].template to<int>();
                    upper[0] = boxDict["upper"]["x"].template to<int>();

                    if constexpr (dimension >= 2)
                    {
                        lower[1] = boxDict["lower"]["y"].template to<int>();
                        upper[1] = boxDict["upper"]["y"].template to<int>();
                    }

                    if constexpr (dimension == 3)
                    {
                        lower[2] = boxDict["lower"]["z"].template to<int>();
                        upper[2] = boxDict["upper"]["z"].template to<int>();
                    }

                    dbBoxes.push_back(SAMRAI::tbox::DatabaseBox(samraiDim, lower, upper));
                }
                levelDB->putDatabaseBoxVector("boxes", dbBoxes);
            }
        } // end loop on levels
        return standardTagInitDB;
    }
    else if (refinement.contains("tagging"))
    {
        standardTagInitDB->putString("tagging_method", "GRADIENT_DETECTOR");
        return standardTagInitDB;
    }
    return nullptr;
}

} // namespace PHARE::amr


#endif // INTEGRATOR_HPP
