#ifndef PHARE_TENSOR_FIELD_RESOURCE_HPP
#define PHARE_TENSOR_FIELD_RESOURCE_HPP

#include "amr/data/tensorfield/tensor_field_data.hpp"
#include "amr/data/tensorfield/tensor_field_variable.hpp"

namespace PHARE
{
namespace amr
{
    /** @brief tells SAMRAI which kind of variable, patchdata are used for a Field Resource
     * also says the type of the actual data buffer
     */
    template<std::size_t rank, typename Grid_t, typename GridLayoutT, typename PhysicalQuantity>
    struct UserTensorFieldType
    {
        using patch_data_type = TensorFieldData<rank, GridLayoutT, Grid_t, PhysicalQuantity>;
        using variable_type   = TensorFieldVariable<rank, GridLayoutT, Grid_t, PhysicalQuantity>;
    };


} // namespace amr
} // namespace PHARE


#endif // PHARE_TENSOR_FIELD_RESOURCE_HPP
