#ifndef PHARE_AMR_PHYSICAL_MODEL_MHD_MODEL_STORAGE_HPP
#define PHARE_AMR_PHYSICAL_MODEL_MHD_MODEL_STORAGE_HPP


#include "core/def.hpp"
#include "core/data/tensorfield/tensorfield.hpp"

#include "amr/resources_manager/resources_manager.hpp"


namespace PHARE::solver
{
// MHD has no particle-driven layout mode to dispatch on (unlike hybrid_model_storage's Ions),
// so tiled-ness is read directly off the state field type.
template<bool is_tiled, typename VecFieldT, typename Grid_t, typename GridLayoutT>
struct mhd_model_storage;

template<typename VecFieldT, typename Grid_t, typename GridLayoutT>
struct mhd_model_storage<true, VecFieldT, Grid_t, GridLayoutT>
{
    using gridlayout_type    = GridLayoutT;
    using tiled_grid_type    = Grid_t;
    using tiled_field_type   = VecFieldT::field_type;
    using PhysicalQuantity_t = VecFieldT::physical_quantity_type;
    using grid_type          = tiled_field_type::grid_type;             // normal non-tiled grid
    using field_type         = tiled_field_type::grid_type::field_type; // normal non-tiled field
    using vecfield_type      = core::TensorField<field_type, PhysicalQuantity_t, 1>;

    using UserField_t = amr::UserFieldType<grid_type, gridlayout_type>;
    using UserVecField_t
        = amr::UserTensorFieldType<1, grid_type, gridlayout_type, PhysicalQuantity_t>;

    using resources_manager_type
        = amr::ResourcesManager<gridlayout_type, Grid_t, UserField_t, UserVecField_t>;

    NO_DISCARD auto getCompileTimeResourcesViewList() const { return std::forward_as_tuple(); }
    NO_DISCARD auto getCompileTimeResourcesViewList() { return std::forward_as_tuple(); }
};

template<typename VecFieldT, typename Grid_t, typename GridLayoutT>
struct mhd_model_storage<false, VecFieldT, Grid_t, GridLayoutT>
{
    using gridlayout_type        = GridLayoutT;
    using resources_manager_type = amr::ResourcesManager<gridlayout_type, Grid_t>;
    using field_type             = VecFieldT::field_type;
    using vecfield_type          = VecFieldT;

    NO_DISCARD auto getCompileTimeResourcesViewList() const { return std::forward_as_tuple(); }
    NO_DISCARD auto getCompileTimeResourcesViewList() { return std::forward_as_tuple(); }
};

} // namespace PHARE::solver

#endif // PHARE_AMR_PHYSICAL_MODEL_MHD_MODEL_STORAGE_HPP
