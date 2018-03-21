#include "basic_hierarchy.h"

BasicHierarchy::BasicHierarchy(std::string const &inputFile)
    : inputDatabase{SAMRAI::tbox::InputManager::getManager()->parseInputFile(inputFile)}
    , patchHierarchyDatabase{inputDatabase->getDatabase("PatchHierarchy")}
    , dimension{static_cast<unsigned short int>(
          inputDatabase->getDatabase("Main")->getInteger("dim"))}
    , gridGeometry{std::make_shared<SAMRAI::geom::CartesianGridGeometry>(
          dimension, "cartesian", inputDatabase->getDatabase("CartesianGridGeometry"))}
    , hierarchy{std::make_shared<SAMRAI::hier::PatchHierarchy>("PatchHierarchy", gridGeometry,
                                                               patchHierarchyDatabase)}

{
}




void BasicHierarchy::init()
{
    SAMRAI::hier::Box boxLevel0{dimension};

    int const ownerRank  = SAMRAI::tbox::SAMRAI_MPI::getSAMRAIWorld().getRank();
    static int counterId = 0;

    auto idLevel0 = SAMRAI::hier::GlobalId{SAMRAI::hier::LocalId{counterId}, ownerRank};
    ++counterId;

    boxLevel0.setBlockId(SAMRAI::hier::BlockId{0});
    boxLevel0.setId(SAMRAI::hier::BoxId{idLevel0});

    boxLevel0.setLower(SAMRAI::hier::Index(0, 0, 0));
    boxLevel0.setUpper(SAMRAI::hier::Index(64, 64, 64));

    SAMRAI::hier::BoxContainer level0Container;

    level0Container.push_back(boxLevel0);

    SAMRAI::hier::BoxLevel level0{level0Container, SAMRAI::hier::IntVector::getOne(dimension),
                                  gridGeometry};
    hierarchy->makeNewPatchLevel(0, level0);
}
