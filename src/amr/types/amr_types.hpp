#ifndef PHARE_AMR_TYPES_HPP
#define PHARE_AMR_TYPES_HPP

#include "core/def/phare_mpi.hpp" // IWYU pragma: keep

#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "core/def.hpp"

namespace PHARE
{
namespace amr
{
    struct SAMRAI_Types
    {
        using patch_t     = SAMRAI::hier::Patch;
        using level_t     = SAMRAI::hier::PatchLevel;
        using hierarchy_t = SAMRAI::hier::PatchHierarchy;

        NO_DISCARD static level_t& getLevel(hierarchy_t const& hierarchy, int levelNumber)
        {
            return *hierarchy.getPatchLevel(levelNumber);
        }
    };

} // namespace amr
} // namespace PHARE

#endif
