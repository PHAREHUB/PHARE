#ifndef PHARE_TESTS_AMR_TOOLS_BUFFER_BASIC_HIERARCHY_HPP
#define PHARE_TESTS_AMR_TOOLS_BUFFER_BASIC_HIERARCHY_HPP

#include "core/def/phare_mpi.hpp"


#include <SAMRAI/geom/CartesianGridGeometry.h>
#include <SAMRAI/hier/PatchHierarchy.h>
#include <SAMRAI/pdat/CellVariable.h>
#include <SAMRAI/tbox/InputManager.h>

#include <memory>

class BasicHierarchy
{
public:
    explicit BasicHierarchy(std::string const& inputFile);

    std::shared_ptr<SAMRAI::tbox::Database> inputDatabase;

    SAMRAI::tbox::Dimension dimension;

    std::shared_ptr<SAMRAI::tbox::Database> patchHierarchyDatabase;
    std::shared_ptr<SAMRAI::geom::CartesianGridGeometry> gridGeometry;
    std::shared_ptr<SAMRAI::hier::PatchHierarchy> hierarchy;

    void init();
};


BasicHierarchy::BasicHierarchy(std::string const& inputFile)
    : inputDatabase{SAMRAI::tbox::InputManager::getManager()->parseInputFile(inputFile)}
    , dimension{static_cast<unsigned short int>(
          inputDatabase->getDatabase("Main")->getInteger("dim"))}
    , patchHierarchyDatabase{inputDatabase->getDatabase("PatchHierarchy")}
    , gridGeometry{std::make_shared<SAMRAI::geom::CartesianGridGeometry>(
          dimension, "cartesian", inputDatabase->getDatabase("CartesianGridGeometry"))}
    , hierarchy{std::make_shared<SAMRAI::hier::PatchHierarchy>("PatchHierarchy", gridGeometry,
                                                               patchHierarchyDatabase)}

{
}




void BasicHierarchy::init()
{
    auto dim
        = static_cast<unsigned short int>(inputDatabase->getDatabase("Main")->getInteger("dim"));
    std::vector<int> lowerIndex, upperIndex;
    for (std::size_t i = 0; i < dim; i++)
    {
        lowerIndex.emplace_back(0);
        upperIndex.emplace_back(64);
    }

    SAMRAI::hier::Box boxLevel0{dimension};

    int const ownerRank  = SAMRAI::tbox::SAMRAI_MPI::getSAMRAIWorld().getRank();
    static int counterId = 0;

    auto idLevel0 = SAMRAI::hier::GlobalId{SAMRAI::hier::LocalId{counterId}, ownerRank};
    ++counterId;

    boxLevel0.setBlockId(SAMRAI::hier::BlockId{0});
    boxLevel0.setId(SAMRAI::hier::BoxId{idLevel0});

    boxLevel0.setLower(SAMRAI::hier::Index(lowerIndex));
    boxLevel0.setUpper(SAMRAI::hier::Index(upperIndex));

    SAMRAI::hier::BoxContainer level0Container;

    level0Container.push_back(boxLevel0);

    SAMRAI::hier::BoxLevel level0{level0Container, SAMRAI::hier::IntVector::getOne(dimension),
                                  gridGeometry};
    hierarchy->makeNewPatchLevel(0, level0);
}


#endif
