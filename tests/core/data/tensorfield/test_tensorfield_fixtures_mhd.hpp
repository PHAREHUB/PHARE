#ifndef PHARE_TEST_CORE_DATA_TEST_TENSORFIELD_FIXTURES_HPP
#define PHARE_TEST_CORE_DATA_TEST_TENSORFIELD_FIXTURES_HPP

#include "core/data/grid/grid.hpp"
#include "core/data/field/field.hpp"
#include "core/mhd/mhd_quantities.hpp"
#include "core/data/tensorfield/tensorfield.hpp"

#include "tests/core/data/field/test_field_fixtures_mhd.hpp"

namespace PHARE::core
{

/*
A UsableTensorFieldMHD is an extension of the TensorField view that owns memory for components and
sets the view pointers. It is useful for tests to easily declare usable (== set views) tensors
*/
template<std::size_t dim, std::size_t rank_ = 2>
class UsableTensorFieldMHD : public TensorField<FieldMHD<dim>, MHDQuantity, rank_>
{
    auto constexpr static N_elements = detail::tensor_field_dim_from_rank<rank_>();

public:
    auto static constexpr dimension = dim;
    using Super                     = TensorField<FieldMHD<dim>, MHDQuantity, rank_>;
    using Grid_t                    = Grid<NdArrayVector<dim>, MHDQuantity::Scalar>;
    using tensor_t                  = typename Super::tensor_t;

    template<typename GridLayout>
    UsableTensorFieldMHD(std::string const& name, GridLayout const& layout, tensor_t qty)
        : Super{name, qty}
        , xyz{make_grids(Super::componentNames(), layout, qty)}
    {
        for (std::size_t i = 0; i < N_elements; ++i)
            super()[i].setBuffer(&xyz[i]);
    }

    void set_on(Super& tensorfield)
    {
        // used for setting on normal model tensorfields
        for (std::size_t i = 0; i < N_elements; ++i)
            tensorfield[i].setBuffer(&xyz[i]);
    }

    Super& super() { return *this; }
    Super const& super() const { return *this; }

protected:
    template<typename ComponentNames, typename GridLayout>
    auto static make_grids(ComponentNames const& compNames, GridLayout const& layout, tensor_t qty)
    {
        auto qts = MHDQuantity::componentsQuantities(qty);
        return for_N<N_elements, for_N_R_mode::make_array>(
            [&](auto i) { return Grid_t{compNames[i], qts[i], layout.allocSize(qts[i])}; });
    }

    std::array<Grid_t, N_elements> xyz;
};


} // namespace PHARE::core


#endif /*PHARE_TEST_CORE_DATA_TEST_TENSORFIELD_FIXTURES_HPP*/
