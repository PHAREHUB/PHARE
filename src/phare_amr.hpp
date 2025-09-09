#ifndef PHARE_AMR_INCLUDE_HPP
#define PHARE_AMR_INCLUDE_HPP


#include "core/def/phare_mpi.hpp" // IWYU pragma: keep


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



#include "phare_core.hpp"
#include "amr/amr_constants.hpp"
#include "amr/types/amr_types.hpp"
#include "amr/wrappers/hierarchy.hpp"
#include "amr/wrappers/integrator.hpp"
#include "amr/messengers/messenger_factory.hpp"
#include "amr/data/particles/refine/splitter.hpp"
#include "amr/data/particles/refine/particles_data_split.hpp"



namespace PHARE::amr
{
template<std::size_t dimension_, std::size_t interp_order_, std::size_t nbRefinedPart_>
struct PHARE_Types
{
    static auto constexpr dimension     = dimension_;
    static auto constexpr interp_order  = interp_order_;
    static auto constexpr nbRefinedPart = nbRefinedPart_;

    using hierarchy_t = PHARE::amr::Hierarchy;

    using Splitter_t = PHARE::amr::Splitter<PHARE::core::DimConst<dimension>,
                                            PHARE::core::InterpConst<interp_order>,
                                            PHARE::core::RefinedParticlesConst<nbRefinedPart>>;

    using core_types = PHARE::core::PHARE_Types<dimension, interp_order>;

    using RefinementParams
        = PHARE::amr::RefinementParams<typename core_types::ParticleArray_t, Splitter_t>;
};

} // namespace PHARE::amr


#endif // PHARE_AMR_INCLUDE_HPP
