#ifndef PHARE_TEST_CORE_DATA_TEST_TENSORFIELD_FIXTURES_HPP
#define PHARE_TEST_CORE_DATA_TEST_TENSORFIELD_FIXTURES_HPP

#include "core/data/field/field.hpp"
#include "core/data/tensorfield/tensorfield.hpp"

namespace PHARE::core
{
template<std::size_t dim>
using Field_t = Field<dim, HybridQuantity::Scalar>;

template<std::size_t dim, std::size_t rank_ = 2>
class UsableTensorField : public TensorField<Field_t<dim>, HybridQuantity, rank_>
{
    auto constexpr static N_elements = detail::tensor_field_dim_from_rank<rank_>();

    using Grid_t   = Grid<NdArrayVector<dim>, HybridQuantity::Scalar>;
    using Super    = TensorField<Field_t<dim>, HybridQuantity, rank_>;
    using tensor_t = typename Super::tensor_t;

public:
    auto static constexpr dimension = dim;

    template<typename GridLayout>
    UsableTensorField(std::string const& name, GridLayout const& layout, tensor_t qty)
        : Super{name, qty}
        , xyz{make_grids(Super::componentNames(), layout, qty)}
    {
        for (std::size_t i = 0; i < N_elements; ++i)
            Super::setBuffer(Super::componentNames()[i], &xyz[i]);
    }

    void set_on(Super& tensorfield)
    {
        // used for setting on normal model tensorfields
        for (std::size_t i = 0; i < N_elements; ++i)
            tensorfield.setBuffer(Super::componentNames()[i], &xyz[i]);
    }

protected:
    template<typename ComponentNames, typename GridLayout>
    auto static make_grids(ComponentNames const& compNames, GridLayout const& layout, tensor_t qty)
    {
        auto qts = HybridQuantity::componentsQuantities(qty);
        return for_N<N_elements, for_N_R_mode::make_array>([&](auto i) {
            return Grid_t{compNames[i], qts[i], layout.allocSize(qts[i])};
        });
    }

    std::array<Grid_t, N_elements> xyz;
};


} // namespace PHARE::core


#endif /*PHARE_TEST_CORE_DATA_TEST_TENSORFIELD_FIXTURES_HPP*/
