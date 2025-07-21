#ifndef PHARE_AMR_INCLUDE_HPP
#define PHARE_AMR_INCLUDE_HPP


#include "core/def/phare_mpi.hpp" // IWYU pragma: keep

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
    auto static constexpr dimension     = opts.dimension;
    auto static constexpr interp_order  = opts.interp_order;
    auto static constexpr nbRefinedPart = opts.nbRefinedPart;

    using core_types = PHARE::core::PHARE_Types<opts>;

    using hierarchy_t = PHARE::amr::Hierarchy;

    using Splitter_t = PHARE::amr::Splitter<PHARE::core::DimConst<dimension>,
                                            PHARE::core::InterpConst<interp_order>,
                                            PHARE::core::RefinedParticlesConst<nbRefinedPart>>;

    using RefinementParams
        = PHARE::amr::RefinementParams<typename core_types::ParticleArray_t, Splitter_t>;
};

} // namespace PHARE::amr


#endif // PHARE_AMR_INCLUDE_HPP
