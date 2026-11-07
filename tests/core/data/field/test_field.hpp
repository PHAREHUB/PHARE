#ifndef PHARE_TEST_CORE_FIELD_TEST_HPP
#define PHARE_TEST_CORE_FIELD_TEST_HPP

#include "core/data/field/field.hpp"
#include "core/models/quantities/hybrid_quantities.hpp"

#include "tests/core/data/field/test_field_fixtures.hpp"

#include "gtest/gtest.h" // EXPECT_FLOAT_EQ

#include <cassert>

namespace PHARE::core
{

template<std::size_t dim>
struct FieldMock
{
    static auto constexpr dimension  = dim;
    auto constexpr static alloc_mode = AllocatorMode::CPU;

    FieldMock() = default;

    auto& operator()(auto&&...) { return data; }
    auto& operator()(auto&&...) const { return data; }
    auto physicalQuantity() const { return qty; }
    std::string name() const { return "FieldMock"; }

    double data;
    HybridQuantity::Scalar qty = HybridQuantity::Scalar::Ex;
};


template<typename GridLayout, typename Field0, typename Field1>
void test_fields(GridLayout const& layout, Field0 const& field0, Field1 const& field1)
{
    EXPECT_EQ(field0.shape(), field1.shape());

    auto const eq = valid_ghost_box(field0.physicalQuantity())
                        ? compare_fields(field0, field1)
                        : compare_field_domains(layout, field0, field1);
    EXPECT_TRUE(eq) << eq.why();
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
    test_fields(layout, field0, make_field_from(field0, fieldV.data()));
}


template<typename GridLayout, typename Field0, typename Field1>
void test(GridLayout const& layout, Field0 const& field0, Field1 const& field1)
{
    test_fields(layout, field0, field1);
}


} // namespace PHARE::core


#endif /* PHARE_TEST_CORE_FIELD_TEST_HPP */
