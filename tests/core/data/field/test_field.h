#ifndef PHARE_TEST_CORE_FIELD_TEST_H
#define PHARE_TEST_CORE_FIELD_TEST_H

#include <cassert>
#include <functional>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace PHARE
{
template<typename FieldFilter, typename Func, typename GridLayout, typename Field>
std::array<std::size_t, GridLayout::dimension> fieldIndices(FieldFilter ff, Func&& func,
                                                            GridLayout& layout, Field& field)
{
    constexpr auto dim = GridLayout::dimension;
    static_assert(dim >= 1 and dim <= 3, "Invalid dimension.");

    auto get = [&](auto dir) { return func(ff, layout, field, dir); };

    if constexpr (dim == 1)
        return {get(core::Direction::X)};
    if constexpr (dim == 2)
        return {get(core::Direction::X), get(core::Direction::Y)};
    if constexpr (dim == 3)
        return {get(core::Direction::X), get(core::Direction::Y), get(core::Direction::Z)};
}


class FieldNullFilter
{
public:
    template<typename Field, typename GridLayout>
    std::size_t start(GridLayout const& layout, Field const& field, core::Direction const direction)
    {
        return layout.ghostStartIndex(field.physicalQuantity(), direction);
    }

    template<typename Field, typename GridLayout>
    std::size_t end(GridLayout const& layout, Field const& field, core::Direction const direction)
    {
        return layout.ghostEndIndex(field.physicalQuantity(), direction);
    }

    template<typename Field, typename GridLayout>
    std::size_t size(GridLayout const& layout, Field const& field, core::Direction const direction)
    {
        return end(layout, field, direction) - start(layout, field, direction) + 1;
    }
};

class FieldDomainPlusNFilter
{
public:
    FieldDomainPlusNFilter(std::size_t n = 0)
        : n_{n}
    {
    }

    template<typename Field, typename GridLayout>
    std::size_t start(GridLayout const& layout, Field const& field, core::Direction const direction)
    {
        return layout.physicalStartIndex(field.physicalQuantity(), direction) - n_;
    }

    template<typename Field, typename GridLayout>
    std::size_t end(GridLayout const& layout, Field const& field, core::Direction const direction)
    {
        return layout.physicalEndIndex(field.physicalQuantity(), direction) + n_;
    }

private:
    std::size_t n_;
};


struct FieldDomainFilter : public FieldDomainPlusNFilter
{
};
} // end namespace PHARE

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

    PHARE::core::HybridQuantity::Scalar qty = PHARE::core::HybridQuantity::Scalar::Ex;
};




template<typename GridLayout, typename Field, typename T1, typename FF = PHARE::FieldNullFilter>
void test_fields(GridLayout const& layout, Field const& field0, T1 const& field1,
                 FF const ff = FF{})
{
    constexpr auto dim = GridLayout::dimension;

    EXPECT_EQ(field0.shape(), field1.shape());

    auto beg
        = fieldIndices(ff, std::mem_fn(&FF::template start<Field, GridLayout>), layout, field0);
    auto end = fieldIndices(ff, std::mem_fn(&FF::template end<Field, GridLayout>), layout, field0);

    if constexpr (dim == 1)
    {
        for (std::size_t i = beg[0]; i < end[0]; ++i)
        {
            if (std::isnan(field0(i)) || std::isnan(field1(i)))
                throw std::runtime_error("This 1dfield should not be NaN index "
                                         + std::to_string(i));
            EXPECT_FLOAT_EQ(field0(i), field1(i));
        }
    }
    else if constexpr (dim == 2)
    {
        for (std::size_t i = beg[0]; i < end[0]; ++i)
        {
            for (std::size_t j = beg[1]; j < end[1]; ++j)
            {
                if (std::isnan(field0(i, j)) || std::isnan(field1(i, j)))
                    throw std::runtime_error("This 2dfield should not be NaN");
                EXPECT_FLOAT_EQ(field0(i, j), field1(i, j));
            }
        }
    }
    else if constexpr (dim == 3)
    {
        for (std::size_t i = beg[0]; i < end[0]; ++i)
        {
            for (std::size_t j = beg[1]; j < end[1]; ++j)
            {
                for (std::size_t k = beg[2]; k < end[2]; ++k)
                {
                    if (std::isnan(field0(i, j, k)) || std::isnan(field1(i, j, k)))
                        throw std::runtime_error("This 3dfield should not be NaN");
                    EXPECT_FLOAT_EQ(field0(i, j, k), field1(i, j, k));
                }
            }
        }
    }
}


template<typename GridLayout, typename NdArrayImpl, typename FF = PHARE::FieldNullFilter>
void test(GridLayout const& layout,
          PHARE::core::Field<NdArrayImpl, PHARE::core::HybridQuantity::Scalar> const& field0,
          PHARE::core::Field<NdArrayImpl, PHARE::core::HybridQuantity::Scalar> const& field1,
          FF const ff = FF{})
{
    test_fields(layout, field0, field1, ff);
}


template<typename GridLayout, typename NdArrayImpl, std::size_t dim, typename T,
         typename FF = PHARE::FieldNullFilter>
void test(GridLayout const& layout,
          PHARE::core::Field<NdArrayImpl, PHARE::core::HybridQuantity::Scalar> const& field0,
          PHARE::core::NdArrayView<dim, T> const& field1, FF const ff = FF{})
{
    static_assert(NdArrayImpl::dimension == dim);
    static_assert(std::is_same_v<typename NdArrayImpl::type, T>);
    test_fields(layout, field0, field1, ff);
}


template<typename GridLayout, typename NdArrayImpl, typename T,
         typename FF = PHARE::FieldNullFilter>
void test(GridLayout const& layout,
          PHARE::core::Field<NdArrayImpl, PHARE::core::HybridQuantity::Scalar> const& field0,
          std::vector<T> const& fieldV, FF const ff = FF{})
{
    static_assert(std::is_same_v<typename NdArrayImpl::type, T>);
    EXPECT_EQ(field0.size(), fieldV.size());
    core::NdArrayView<GridLayout::dimension, T, T const* const> field1{fieldV, field0.shape()};
    test_fields(layout, field0, field1, ff);
}

template<typename GridLayout, typename NdArrayView, typename Qty, typename T,
         typename FF = PHARE::FieldNullFilter>
void test(GridLayout const& layout, PHARE::core::FieldView<NdArrayView, Qty> field0,
          std::vector<T>&& fieldV, FF const ff = FF{})
{
    EXPECT_EQ(field0.size(), fieldV.size());
    core::NdArrayView<GridLayout::dimension, T> field1{fieldV, field0.shape()};
    test_fields(layout, field0, field1, ff);
}


template<typename GridLayout, typename NdArrayView, typename Qty,
         typename FF = PHARE::FieldNullFilter>
void test(GridLayout const& layout, PHARE::core::FieldView<NdArrayView, Qty> field0,
          PHARE::core::FieldView<NdArrayView, Qty> field1, FF const ff = FF{})
{
    test_fields(layout, field0, field1, ff);
}


} // namespace PHARE::core




#endif /* PHARE_TEST_CORE_FIELD_TEST_H */
