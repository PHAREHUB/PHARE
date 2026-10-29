#ifndef PHARE_AMR_DATA_FIELD_FIELD_LAYOUT_DATA_HPP
#define PHARE_AMR_DATA_FIELD_FIELD_LAYOUT_DATA_HPP

#include <SAMRAI/hier/Patch.h>

namespace PHARE::amr
{
// patch data storing a GridLayout built from its patch at allocation. On restarts, SAMRAI
//  allocates patch data before restored patches have a geometry, so such layouts get the
//  fallback mesh size and origin of layoutFromPatch and must be rebuilt once it is set
//  see https://github.com/PHAREHUB/PHARE/issues/1518
class FieldLayoutData
{
public:
    virtual ~FieldLayoutData() = default;

    virtual void resetLayout(SAMRAI::hier::Patch const& patch) = 0;
};

} // namespace PHARE::amr

#endif // PHARE_AMR_DATA_FIELD_FIELD_LAYOUT_DATA_HPP
