#ifndef PHARE_AMR_INCLUDE_HPP
#define PHARE_AMR_INCLUDE_HPP


#include "phare_mpi.hpp" // IWYU pragma: keep

#include "phare_core.hpp"
#include "amr/wrappers/hierarchy.hpp"
#include "amr/data/particles/refine/splitter.hpp"
#include "amr/data/particles/refine/particles_data_split.hpp"


#include <SAMRAI/hier/Box.h>
#include <SAMRAI/tbox/Database.h>
#include <SAMRAI/hier/IntVector.h>
#include <SAMRAI/tbox/DatabaseBox.h>
#include <SAMRAI/hier/BoxContainer.h>
#include <SAMRAI/tbox/InputManager.h>
#include <SAMRAI/tbox/MemoryDatabase.h>
#include <SAMRAI/hier/PatchHierarchy.h>
#include <SAMRAI/mesh/BergerRigoutsos.h>
#include <SAMRAI/mesh/GriddingAlgorithm.h>
#include <SAMRAI/mesh/TreeLoadBalancer.h>
#include <SAMRAI/geom/CartesianGridGeometry.h>
#include <SAMRAI/mesh/StandardTagAndInitialize.h>
#include <SAMRAI/algs/TimeRefinementIntegrator.h>




namespace PHARE::amr
{
template<SimOpts opts>
struct PHARE_Types
{
    auto static constexpr dimension = opts.dimension;

    using hierarchy_t = PHARE::amr::Hierarchy;
};

} // namespace PHARE::amr


#endif // PHARE_AMR_INCLUDE_HPP
