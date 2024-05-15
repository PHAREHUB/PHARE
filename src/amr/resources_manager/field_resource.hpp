#ifndef PHARE_FIELD_RESOURCE_HPP
#define PHARE_FIELD_RESOURCE_HPP

#include "amr/data/field/field_data.hpp"
#include "amr/data/field/field_variable.hpp"

namespace PHARE
{
namespace amr
{
    /** @brief tells SAMRAI which kind of variable, patchdata are used for a Field Resource
     * also says the type of the actual data buffer
     */
    template<typename Grid_t, typename GridLayoutT>
    struct UserFieldType
    {
        using patch_data_type = FieldData<GridLayoutT, Grid_t>;
        using variable_type   = FieldVariable<GridLayoutT, Grid_t>;
    };


} // namespace amr
} // namespace PHARE


#endif // FIELD_RESOURCE_HPP
