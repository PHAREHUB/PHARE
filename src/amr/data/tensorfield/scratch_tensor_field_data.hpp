#ifndef PHARE_AMR_DATA_TENSORFIELD_SCRATCH_TENSOR_FIELD_DATA_HPP
#define PHARE_AMR_DATA_TENSORFIELD_SCRATCH_TENSOR_FIELD_DATA_HPP

#include "amr/data/field/scratch_field_data.hpp"
#include "amr/resources_manager/tensor_field_resource.hpp"
#include "amr/data/tensorfield/tensor_field_data_factory.hpp"

#include <cstddef>

// scratch tensor fields, see scratch_field_data.hpp
namespace PHARE::amr
{
template<typename View_t, typename Grid_t, typename GridLayoutT, typename PhysicalQuantity>
    requires requires { View_t::rank; }
struct UserScratchFieldType<View_t, Grid_t, GridLayoutT, PhysicalQuantity>
{
    std::size_t static constexpr rank = View_t::rank;
    using user_type       = UserTensorFieldType<rank, Grid_t, GridLayoutT, PhysicalQuantity>;
    using factory_type    = TensorFieldDataFactory<rank, GridLayoutT, Grid_t, PhysicalQuantity>;
    using patch_data_type = user_type::patch_data_type;
    using variable_type   = ScratchVariable<typename user_type::variable_type, factory_type,
                                            GridLayoutT, View_t::scratch_ghosts>;
};

} // namespace PHARE::amr

#endif /* PHARE_AMR_DATA_TENSORFIELD_SCRATCH_TENSOR_FIELD_DATA_HPP */
