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

#endif
