#ifndef PHARE_TEST_CORE_DATA_TEST_TENSORFIELD_FIXTURES_HPP
#define PHARE_TEST_CORE_DATA_TEST_TENSORFIELD_FIXTURES_HPP

#include "core/def/phare_config.hpp"

#include "phare_core.hpp"
// #include "core/data/field/field.hpp"
// #include "core/data/grid/gridlayoutdefs.hpp"
// #include "core/data/particles/particle_array_def.hpp"

// #include "core/data/grid/grid.hpp"
#include "core/utilities/types.hpp"
#include "core/utilities/equality.hpp"
#include "core/hybrid/hybrid_quantities.hpp"
// #include "core/data/ndarray/ndarray_vector.hpp"
#include "core/data/tensorfield/tensorfield.hpp"


#include "tests/core/data/grid/test_grid_fixtures.hpp"
#include "tests/core/data/field/test_field_fixtures.hpp"

#include <cstddef>

namespace PHARE::core
{


template<typename GridLayout_, std::size_t rank, auto alloc_mde, auto layout_mde>
struct UsableTensorFieldImpl
{
    using Resolver_t = UsingResolver<GridLayout_, layout_mde, alloc_mde>;
    using Grid_t     = Resolver_t::Grid_t;
    using Field_t    = Resolver_t::Field_t;
    using Super_t    = TensorField<Field_t, HybridQuantity, rank>;
};




/*
A UsableTensorField is an extension of the TensorField view that owns memory for components and sets
the view pointers. It is useful for tests to easily declare usable (== set views) tensors

Note: UsableTensorFields hold Grids that are default initialized to zero for convenience rather
than NaN (default grid init value)

*/
template<typename GridLayout_, std::size_t rank_ = 2, auto alloc_mde = AllocatorMode::CPU,
         auto layout_mde = LayoutMode::AoS>
class UsableTensorField
    : public UsableTensorFieldImpl<GridLayout_, rank_, alloc_mde, layout_mde>::Super_t
{
    // auto constexpr static alloc_mode = alloc_mde;
    static_assert(std::is_same_v<decltype(alloc_mde), AllocatorMode>);
    static_assert(std::is_same_v<decltype(layout_mde), LayoutMode>);

    auto constexpr static N_elements = detail::tensor_field_dim_from_rank<rank_>();
    using Impl_t = UsableTensorFieldImpl<GridLayout_, rank_, alloc_mde, layout_mde>;
    using This   = UsableTensorField;

public:
    auto static constexpr dimension = GridLayout_::dimension;
    using Super                     = Impl_t::Super_t;
    using Grid_t                    = Impl_t::Grid_t;
    using Field_t                   = Impl_t::Field_t;
    auto static constexpr rank      = Super::rank;


protected:
    using tensor_t = Super::tensor_t;

    void _set()
    {
        for (std::size_t i = 0; i < N_elements; ++i)
            super()[i].setBuffer(&xyz[i]);
    }



public:
    UsableTensorField(std::string const& name, GridLayout_ const& layout, tensor_t qty,
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

    Super& view() { return *this; }
    Super const& view() const { return *this; }

    auto& operator*() { return view(); }
    auto& operator*() const { return view(); }

    auto& grids() const { return xyz; }

    auto& operator[](std::size_t const i) _PHARE_ALL_FN_ { return xyz[i]; }
    auto& operator[](std::size_t const i) const _PHARE_ALL_FN_ { return xyz[i]; }


    template<auto AM>
    bool isclose(UsableTensorField<GridLayout_, rank, AM, layout_mde> const& that,
                 double diff = 1e-15) const
    {
        return core::for_N_all<N_elements>(
            [&](auto i) { return (*this)[i].isclose(that[i], diff); });
    }

    template<auto AM>
    bool operator==(UsableTensorField<GridLayout_, rank, AM, layout_mde> const& that) const
    {
        return core::for_N_all<N_elements>([&](auto i) { return (*this)[i] == that[i]; });
    }

    void fill(double const v)
    {
        for (auto& c : xyz)
            c.fill(v);
    }

    template<typename GL, std::size_t rank, auto am, auto lm>
    friend std::ostream& operator<<(std::ostream&, UsableTensorField const&);

protected:
    template<typename ComponentNames, typename GridLayout>
    auto static make_grids(ComponentNames const& compNames, GridLayout const& layout, tensor_t qty)
    {
        auto qts = HybridQuantity::componentsQuantities(qty);
        return for_N<N_elements, for_N_R_mode::make_array>(
            [&](auto i) { return Grid_t{compNames[i], layout, qts[i], 0}; });
    }

    std::array<Grid_t, N_elements> xyz;
};




template<typename GL, std::size_t R, auto am0, auto lm0, auto am1, auto lm1>
EqualityReport compare_reduced_tensor_fields(UsableTensorField<GL, R, am0, lm0> const& ref,
                                             UsableTensorField<GL, R, am1, lm1> const& cmp,
                                             double const diff)
{
    auto constexpr static N_elements = detail::tensor_field_dim_from_rank<R>();

    std::stringstream log;
    log << std::endl;
    for (std::size_t ci = 0; ci < N_elements; ++ci)
        if (auto eqr = compare_reduced_fields(ref[ci], cmp[ci], diff); !eqr)
            return eqr;
        else
            log << eqr.what() << std::endl;

    return EqualityReport{true, log.str()};
}




template<typename F0, typename F1, typename PhysicalQuantity, std::size_t rank>
EqualityReport compare_reduced_tensor_fields(TensorField<F0, PhysicalQuantity, rank> const& ref,
                                             TensorField<F1, PhysicalQuantity, rank> const& cmp,
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
        if (auto eqr = compare_reduced_fields(ref[ci], cmp[ci], diff); !eqr)
            return eqr;
        else
            log << eqr.what() << std::endl;

    return EqualityReport{true, log.str()};
}




template<typename GL, std::size_t R, auto am0, auto lm0, auto am1, auto lm1>
EqualityReport compare_tensor_fields(UsableTensorField<GL, R, am0, lm0> const& ref,
                                     UsableTensorField<GL, R, am1, lm1> const& cmp,
                                     double const diff)
{
    auto constexpr static N_elements = detail::tensor_field_dim_from_rank<R>();

    std::stringstream log;
    log << std::endl;
    for (std::size_t ci = 0; ci < N_elements; ++ci)
        if (auto eqr = compare_fields(ref[ci], cmp[ci], diff); !eqr)
            return eqr;
        else
            log << eqr.what() << std::endl;

    return EqualityReport{true, log.str()};
}




template<typename F0, typename F1, typename PhysicalQuantity, std::size_t rank>
EqualityReport compare_tensor_fields(TensorField<F0, PhysicalQuantity, rank> const& ref,
                                     TensorField<F1, PhysicalQuantity, rank> const& cmp,
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


template<typename GL, std::size_t rank, auto am, auto lm>
std::ostream& operator<<(std::ostream& out, UsableTensorField<GL, rank, am, lm> const& ts)
{
    return (out << *ts);
}


template<typename F0, typename PhysicalQuantity, std::size_t rank>
void zero_ghost_layer(TensorField<F0, PhysicalQuantity, rank>& tf, auto const& layout)
{
    for (auto& f : tf)
        zero_ghost_layer(f, layout);
}

void zero_ghost_layer([[maybe_unused]] auto& layout, auto&... args)
{
    (zero_ghost_layer(args), ...);
}

} // namespace PHARE::core


#endif /*PHARE_TEST_CORE_DATA_TEST_TENSORFIELD_FIXTURES_HPP*/
