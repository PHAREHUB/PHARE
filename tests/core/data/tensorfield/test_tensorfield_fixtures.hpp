#ifndef PHARE_TEST_CORE_DATA_TEST_TENSORFIELD_FIXTURES_HPP
#define PHARE_TEST_CORE_DATA_TEST_TENSORFIELD_FIXTURES_HPP

#include "core/data/grid/grid.hpp"
#include "core/data/field/field.hpp"
#include "core/utilities/equality.hpp"
#include "core/hybrid/hybrid_quantities.hpp"
#include "core/data/tensorfield/tensorfield.hpp"

#include "tests/core/data/field/test_field_fixtures.hpp"

namespace PHARE::core
{

/*
A UsableTensorField is an extension of the TensorField view that owns memory for components and sets
the view pointers. It is useful for tests to easily declare usable (== set views) tensors

Note: UsableTensorFields hold Grids that are default initialized to zero for convenience rather
than NaN (default grid init value)

*/
template<std::size_t dim, std::size_t rank_ = 2>
class UsableTensorField : public TensorField<Field_t<dim>, HybridQuantity, rank_>
{
    auto constexpr static N_elements = detail::tensor_field_dim_from_rank<rank_>();
    void _set()
    {
        for (std::size_t i = 0; i < N_elements; ++i)
            super()[i].setBuffer(&xyz[i]);
    }

public:
    auto static constexpr dimension = dim;
    using Super                     = TensorField<Field_t<dim>, HybridQuantity, rank_>;
    using Grid_t                    = Grid<NdArrayVector<dim>, HybridQuantity::Scalar>;
    using tensor_t                  = typename Super::tensor_t;

    template<typename GridLayout>
    UsableTensorField(std::string const& name, GridLayout const& layout, tensor_t qty,
                      std::optional<double> v = std::nullopt)
        : Super{name, qty}
        , xyz{make_grids(Super::componentNames(), layout, qty)}
    {
        if (v)
            for (std::size_t i = 0; i < N_elements; ++i)
                xyz[i].fill(*v);
        _set();
    }

    UsableTensorField(UsableTensorField&& that)
        : Super{std::forward<Super>(that)}
        , xyz{std::move(that.xyz)}
    {
        _set();
    }

    UsableTensorField(UsableTensorField const& that)
        : Super{that}
        , xyz{that.xyz}
    {
        _set();
    }

    void set_on(Super& tensorfield)
    {
        // used for setting on normal model tensorfields
        for (std::size_t i = 0; i < N_elements; ++i)
            tensorfield[i].setBuffer(&xyz[i]);
    }

    Super& super() { return *this; }
    Super& super() const { return *this; }

protected:
    template<typename ComponentNames, typename GridLayout>
    auto static make_grids(ComponentNames const& compNames, GridLayout const& layout, tensor_t qty)
    {
        auto qts = HybridQuantity::componentsQuantities(qty);
        return for_N<N_elements, for_N_R_mode::make_array>(
            [&](auto i) { return Grid_t{compNames[i], qts[i], layout.allocSize(qts[i]), 0.}; });
    }

    std::array<Grid_t, N_elements> xyz;
};




template<typename T, typename PhysicalQuantity, std::size_t rank>
EqualityReport compare_tensor_fields(TensorField<T, PhysicalQuantity, rank> const& ref,
                                     TensorField<T, PhysicalQuantity, rank> const& cmp,
                                     double const diff)
{
    auto constexpr static N_elements = detail::tensor_field_dim_from_rank<rank>();

    if (ref.componentNames() != cmp.componentNames())
        return EqualityReport{false, "Tensorfield component mismatch"};

    auto const same_sizes = [&]() {
        return core::for_N_all<N_elements>([&](auto i) { return ref[i].size() == cmp[i].size(); });
    }();

    if (!same_sizes)
        return EqualityReport{false, "Tensorfield shape/size mismatch"};

    std::stringstream log;
    log << std::endl;
    for (std::size_t ci = 0; ci < N_elements; ++ci)
        if (auto eqr = compare_fields(ref[ci], cmp[ci], diff); !eqr)
            return eqr;
        else
            log << eqr.what() << std::endl;

    return EqualityReport{true, log.str()};
}



} // namespace PHARE::core


#endif /*PHARE_TEST_CORE_DATA_TEST_TENSORFIELD_FIXTURES_HPP*/
