#ifndef PHARE_TEST_CORE_DATA_TEST_TENSORFIELD_FIXTURES_HPP
#define PHARE_TEST_CORE_DATA_TEST_TENSORFIELD_FIXTURES_HPP

#include "phare_core.hpp"
#include "core/def/phare_config.hpp"

#include "core/utilities/types.hpp"
#include "core/utilities/equality.hpp"
#include "core/data/tensorfield/tensorfield.hpp"
#include "core/models/quantities/hybrid_quantities.hpp"

#include "tests/core/data/grid/test_grid_fixtures.hpp"
#include "tests/core/data/field/test_field_fixtures.hpp"

#include <cstddef>

namespace PHARE::core
{


template<typename core_types>
struct TensorFieldOptions
{
    using GridLayout_t         = core_types::GridLayout_t;
    using Grid_t               = core_types::Grid_t;
    using Field_t              = core_types::Field_t;
    auto constexpr static opts = core_types::opts;
    std::size_t rank           = 2;

    auto constexpr with_rank(std::size_t const r) const
    {
        auto copy = *this;
        copy.rank = r;
        return copy;
    }
};


template<auto opts>
struct UsableTensorFieldImpl
{
    using Options      = decltype(opts);
    using GridLayout_t = Options::GridLayout_t;
    using Grid_t       = Options::Grid_t;
    using Field_t      = Options::Field_t;
    using Super_t      = TensorField<Field_t, HybridQuantity, opts.rank>;
};




/*
A UsableTensorField is an extension of the TensorField view that owns memory for components and sets
the view pointers. It is useful for tests to easily declare usable (== set views) tensors

Note: UsableTensorFields hold Grids that are default initialized to zero for convenience rather
than NaN (default grid init value)

*/
template<auto opts>
class UsableTensorField : public UsableTensorFieldImpl<opts>::Super_t
{
    using Options = decltype(opts);

    auto constexpr static N_elements = detail::tensor_field_dim_from_rank<opts.rank>();
    using Impl_t                     = UsableTensorFieldImpl<opts>;
    using This                       = UsableTensorField;

public:
    using GridLayout_               = Options::GridLayout_t;
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


    template<auto o>
    bool isclose(UsableTensorField<o> const& that, double diff = 1e-15) const
    {
        return core::for_N_all<N_elements>(
            [&](auto i) { return (*this)[i].isclose(that[i], diff); });
    }

    template<auto o>
    bool operator==(UsableTensorField<o> const& that) const
    {
        return core::for_N_all<N_elements>([&](auto i) { return (*this)[i] == that[i]; });
    }

    void fill(double const v)
    {
        for (auto& c : xyz)
            c.fill(v);
    }

    template<auto o>
    friend std::ostream& operator<<(std::ostream&, UsableTensorField<o> const&);

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




template<auto o0, auto o1>
EqualityReport compare_reduced_tensor_fields(UsableTensorField<o0> const& ref,
                                             UsableTensorField<o1> const& cmp, double const diff)
{
    static_assert(o0.rank == o1.rank);
    auto constexpr static N_elements = detail::tensor_field_dim_from_rank<o0.rank>();

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




template<auto o0, auto o1>
EqualityReport compare_tensor_fields(UsableTensorField<o0> const& ref,
                                     UsableTensorField<o1> const& cmp, double const diff)
{
    static_assert(o0.rank == o1.rank);
    auto constexpr static N_elements = detail::tensor_field_dim_from_rank<o0.rank>();

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


template<auto o>
std::ostream& operator<<(std::ostream& out, UsableTensorField<o> const& ts)
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
