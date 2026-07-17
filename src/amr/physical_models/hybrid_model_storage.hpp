#ifndef PHARE_AMR_PHYSICAL_MODEL_HYBRID_MODEL_STORAGE_HPP
#define PHARE_AMR_PHYSICAL_MODEL_HYBRID_MODEL_STORAGE_HPP


#include "core/def.hpp"
// #include "core/data/field/field.hpp"
#include "core/data/tensorfield/tensorfield.hpp"
#include "core/data/particles/particle_array_def.hpp"
#include "amr/resources_manager/resources_manager.hpp"
// #include "core/hybrid/hybrid_quantities.hpp"



namespace PHARE::solver
{
template<auto layout_mode, typename Ions, typename Grid_t>
struct hybrid_model_storage;

template<typename TiledIons, typename Grid_t>
struct hybrid_model_storage<core::LayoutMode::AoSTS, TiledIons, Grid_t>
{
    auto static constexpr dim = TiledIons::dimension;
    using gridlayout_type     = TiledIons::gridlayout_type;
    using tiled_grid_type     = Grid_t;
    using tiled_field_type    = TiledIons::field_type;
    using PhysicalQuantity_t  = TiledIons::tensorfield_type::physical_quantity_type;
    using grid_type           = tiled_field_type::grid_type;             // normal non-tiled grid
    using field_type          = tiled_field_type::grid_type::field_type; // normal non-tiled field
    using vecfield_type       = core::TensorField<field_type, PhysicalQuantity_t, 1>;
    using tensorfield_type    = core::TensorField<field_type, PhysicalQuantity_t, 2>;


    using UserField_t = amr::UserFieldType<grid_type, gridlayout_type>;
    using UserVecField_t
        = amr::UserTensorFieldType<1, grid_type, gridlayout_type, core::HybridQuantity>;
    using UserTensorField_t
        = amr::UserTensorFieldType<2, grid_type, gridlayout_type, core::HybridQuantity>;

    using resources_manager_type = amr::ResourcesManager<gridlayout_type, Grid_t, UserField_t,
                                                         UserVecField_t, UserTensorField_t>;

    NO_DISCARD auto getCompileTimeResourcesViewList() const
    {
        return std::forward_as_tuple(scratch_field, scratch_vecfield, scratch_tensorfield);
    }
    NO_DISCARD auto getCompileTimeResourcesViewList()
    {
        return std::forward_as_tuple(scratch_field, scratch_vecfield, scratch_tensorfield);
    }

    field_type scratch_field{"phare_scratch_field", PhysicalQuantity_t::Scalar::Vx};
    vecfield_type scratch_vecfield{"phare_scratch_vec_field", PhysicalQuantity_t::Vector::V};
    tensorfield_type scratch_tensorfield{"phare_scratch_tensor_field",
                                         PhysicalQuantity_t::Tensor::M};
};

template<typename TiledIons, typename Grid_t>
struct hybrid_model_storage<core::LayoutMode::AoSPCTS, TiledIons, Grid_t>
    : hybrid_model_storage<core::LayoutMode::AoSTS, TiledIons, Grid_t>
{
};

template<auto layout_mode, typename Ions, typename Grid_t>
struct hybrid_model_storage
{
    using gridlayout_type        = Ions::gridlayout_type;
    using resources_manager_type = amr::ResourcesManager<gridlayout_type, Grid_t>;
    using grid_type              = Grid_t;           // normal non-tiled grid
    using field_type             = Ions::field_type; // normal non-tiled field
    using vecfield_type          = Ions::vecfield_type;
    using tensorfield_type       = Ions::tensorfield_type;


    NO_DISCARD auto getCompileTimeResourcesViewList() const { return std::forward_as_tuple(); }
    NO_DISCARD auto getCompileTimeResourcesViewList() { return std::forward_as_tuple(); }
};

} // namespace PHARE::solver

#endif // PHARE_AMR_PHYSICAL_MODEL_HYBRID_MODEL_STORAGE_HPP
