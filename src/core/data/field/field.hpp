
#ifndef PHARE_CORE_DATA_FIELD_FIELD_BASE_HPP
#define PHARE_CORE_DATA_FIELD_FIELD_BASE_HPP

#include <array>
#include <string>
#include <vector>
#include <cassert>
#include <cstddef>
#include <utility>
#include <algorithm>

#include "core/data/ndarray/ndarray_vector.hpp"

namespace PHARE::core
{
template<typename NdArray, typename PhysicalQuantity>
class FieldView : public NdArray
{
public:
    static constexpr bool is_contiguous = true;
    static const std::size_t dimension  = NdArray::dimension;
    using type                          = typename NdArray::type;
    using Super                         = NdArray;
    using pointer_type                  = type*;
    using data_view_t                   = typename NdArray::view_t;
    using view_t                        = FieldView<data_view_t, PhysicalQuantity>;
    using Super::data;
    using Super::shape;
    using Super::size;

    template<typename Type>
    FieldView(std::array<Type, dimension> shape, PhysicalQuantity qty)
        : Super{shape}
        , qty_{qty}
    {
    }

    FieldView(pointer_type ptr, std::array<std::uint32_t, dimension> shape, PhysicalQuantity qty)
        : Super{ptr, shape}
        , qty_{qty}
    {
    }

    FieldView(NdArray&& ndArray, PhysicalQuantity qty)
        : Super{std::move(ndArray)}
        , qty_{qty}
    {
    }

    FieldView(NdArray& ndArray, PhysicalQuantity qty)
        : Super{ndArray}
        , qty_{qty}
    {
    }

    constexpr PhysicalQuantity physicalQuantity() const { return qty_; }

private:
    PhysicalQuantity qty_;
};



//! Class Field represents a multidimensional (1,2 or 3D) scalar field
/** Users of Field objects needing to know which physical quantity a specific
 *  Field instance represents can get this info by calling physicalQuantity().
 *  Users may also give a string name to a field object and get a name by calling
 *  name().
 */
template<typename NdArrayImpl, typename PhysicalQuantity>
class Field : public FieldView<NdArrayImpl, PhysicalQuantity>
{
public:
    using Super                  = FieldView<NdArrayImpl, PhysicalQuantity>;
    using impl_type              = NdArrayImpl;
    using type                   = typename NdArrayImpl::type;
    using physical_quantity_type = PhysicalQuantity;
    using data_view_t            = typename NdArrayImpl::view_t;
    using view_t                 = FieldView<data_view_t, PhysicalQuantity>;
    using Super::data;
    using Super::physicalQuantity;
    using Super::shape;
    using Super::size;

    Field()                    = delete;
    Field(Field const& source) = delete;
    Field(Field&& source)      = default;
    Field& operator=(Field&& source) = delete;
    Field& operator=(Field const& source) = delete;

    template<typename... Dims>
    Field(std::string const& name, PhysicalQuantity qty, Dims... dims)
        : Super{std::array{dims...}, qty}
        , name_{name}
    {
        static_assert(sizeof...(Dims) == NdArrayImpl::dimension, "Invalid dimension");
    }

    template<std::size_t dim>
    Field(std::string const& name, PhysicalQuantity qty, std::array<std::uint32_t, dim> const& dims)
        : Super{dims, qty}
        , name_{name}
    {
    }

    Field(Super ndArray, std::string const& name, PhysicalQuantity qty)
        : Super{ndArray, qty}
        , name_{name}
    {
    }

    void copyData(Field const& source)
    {
        static_cast<NdArrayImpl&>(*this) = static_cast<NdArrayImpl const&>(source);
    }

    std::string name() const { return name_; }


    auto view() const { return view_t{data_view_t{data(), shape()}, physicalQuantity()}; }
    auto view() { return view_t{data_view_t{data(), shape()}, physicalQuantity()}; }

private:
    std::string name_{"No Name"};
};




template<typename NdArrayImpl, typename PhysicalQuantity>
void average(Field<NdArrayImpl, PhysicalQuantity> const& f1,
             Field<NdArrayImpl, PhysicalQuantity> const& f2,
             Field<NdArrayImpl, PhysicalQuantity>& avg)
{
    assert(f1.shape() == f2.shape() and f1.shape() == avg.shape());


    std::transform(std::begin(f1), std::end(f1), std::begin(f2), std::begin(avg),
                   std::plus<double>());

    std::transform(std::begin(avg), std::end(avg), std::begin(avg),
                   [](double x) { return x * 0.5; });
}


} // namespace PHARE::core

#endif
