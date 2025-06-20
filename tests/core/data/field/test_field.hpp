#ifndef PHARE_TEST_CORE_FIELD_TEST_HPP
#define PHARE_TEST_CORE_FIELD_TEST_HPP

#include <cassert>
#include <functional>

#include "core/data/ndarray/ndarray_vector.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"


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




template<typename GridLayout, typename Field, typename T1>
void test_fields(GridLayout const& layout, Field const& field0, T1 const& field1)
{
    constexpr auto dim = GridLayout::dimension;

    EXPECT_EQ(field0.shape(), field1.shape());


    for (std::size_t i = 0; i < field0.size(); ++i)
    {
        auto const& v0 = field0.data()[i];
        auto const& v1 = field1.data()[i];
        if (std::isnan(v0) || std::isnan(v1))
            throw std::runtime_error("This 1dfield should not be NaN index " + std::to_string(i));
        EXPECT_FLOAT_EQ(v0, v1);
    }
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
