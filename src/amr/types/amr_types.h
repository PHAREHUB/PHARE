#ifndef PHARE_AMR_TYPES_H
#define PHARE_AMR_TYPES_H


#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchLevel.h"

namespace PHARE
{
namespace amr
{
    struct SAMRAI_Types
    {
        using patch_t = SAMRAI::hier::Patch;
        using level_t = SAMRAI::hier::PatchLevel;
    };
} // namespace amr
} // namespace PHARE

#endif
