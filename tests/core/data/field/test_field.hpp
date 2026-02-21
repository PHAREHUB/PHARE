#ifndef PHARE_TEST_CORE_FIELD_TEST_HPP
#define PHARE_TEST_CORE_FIELD_TEST_HPP

#include "core/utilities/types.hpp"
#include "core/data/field/field.hpp"
#include "core/hybrid/hybrid_quantities.hpp"
#include "core/data/ndarray/ndarray_vector.hpp"

#include "gtest/gtest.h" // EXPECT_FLOAT_EQ

#include <cassert>




namespace PHARE::core
{
template<std::size_t dim>
struct FieldMock
{
    static auto constexpr dimension   = dim;
    static auto constexpr is_host_mem = true;
    double data;

    FieldMock() = default;

    template<typename... Args>
    auto& operator()(Args...)
    {
        return data;
    }
    template<typename... Args>
    auto& operator()(Args...) const
    {
        return data;
    }
    auto physicalQuantity() const { return qty; }
    std::string name() const { return "FieldMock"; }

    HybridQuantity::Scalar qty = HybridQuantity::Scalar::Ex;
};


template<typename PQ>
bool valid_ghost_box(PQ const physicalQuantity)
{
    if constexpr (std::is_same_v<PQ, HybridQuantity::Scalar>)
    {
        using enum HybridQuantity::Scalar;
        // cause last ghost can have no interpolation
        return not any_in(physicalQuantity, rho, Vx, Vy, Vz);
    }
    throw std::runtime_error("No other impl");
}


template<typename GridLayout, typename Field, typename T1>
void test_fields(GridLayout const& layout, Field const& field0, T1 const& field1)
{
    constexpr auto dim = GridLayout::dimension;

    EXPECT_EQ(field0.shape(), field1.shape());

    auto not_nan = [&](auto const& f0, auto const& f1, auto&... idxs) {
        auto const& v0 = f0(idxs...);
        auto const& v1 = f1(idxs...);

        std::ostringstream oss;
        oss << "(";
        ((oss << idxs << ","), ...);
        std::string idx_str = oss.str();
        if constexpr (sizeof...(idxs) > 0)
            idx_str.back() = ')';
        else
            idx_str += ')';

        if (std::isnan(v0) || std::isnan(v1))
            throw std::runtime_error("This 1dfield should not be NaN at index " + idx_str);

        EXPECT_FLOAT_EQ(v0, v1) << " at index " << idx_str;
    };

    if (valid_ghost_box(field0.physicalQuantity()))
        layout.evalOnGhostBox(field0, [&](auto&... idxs) { not_nan(field0, field1, idxs...); });
    else
        layout.evalOnBox(field0, [&](auto&... idxs) { not_nan(field0, field1, idxs...); });
}


template<typename GridLayout, typename NdArrayImpl>
void test(GridLayout const& layout,
          Field<GridLayout::dimension, HybridQuantity::Scalar> const& field0,
          Field<GridLayout::dimension, HybridQuantity::Scalar> const& field1)
{
    test_fields(layout, field0, field1);
}



template<typename GridLayout, typename Field, typename T>
void test(GridLayout const& layout, Field const& field0, std::vector<T> const& fieldV)
{
    EXPECT_EQ(field0.size(), fieldV.size());
    core::NdArrayView<GridLayout::dimension, T const> const field1{fieldV.data(), field0.shape()};
    test_fields(layout, field0, field1);
}



template<typename GridLayout, typename Field0, typename Field1>
void test(GridLayout const& layout, Field0 const& field0, Field1 const& field1)
{
    test_fields(layout, field0, field1);
}


} // namespace PHARE::core


#endif /* PHARE_TEST_CORE_FIELD_TEST_H */
