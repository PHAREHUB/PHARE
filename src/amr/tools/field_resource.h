#ifndef PHARE_FIELD_RESOURCE_H
#define PHARE_FIELD_RESOURCE_H

#include "data/field/field_data.h"
#include "data/field/field_variable.h"

namespace PHARE
{
/** @brief tells SAMRAI which kind of variable, patchdata are used for a Field Resource
 * also says the type of the actual data buffer
 */
template<typename GridLayoutT, typename ResourcesUser>
struct UserFieldType
{
    using patch_data_type   = FieldData<GridLayoutT, typename ResourcesUser::field_type>;
    using variable_type     = FieldVariable<GridLayoutT, typename ResourcesUser::field_type>;
    using internal_type_ptr = typename ResourcesUser::field_type*;
};


} // namespace PHARE


#endif // FIELD_RESOURCE_H
