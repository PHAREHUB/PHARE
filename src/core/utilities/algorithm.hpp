#ifndef PHARE_CORE_UTILITIES_ALGORITHM_HPP
#define PHARE_CORE_UTILITIES_ALGORITHM_HPP


#include "core/def.hpp"

#include "core/utilities/span.hpp"
#include "core/utilities/types.hpp"
#include "core/data/field/field.hpp"
#include "core/data/field/recentered_field.hpp"
#include "core/utilities/box/box.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/data/ndarray/ndarray_vector.hpp"
#include "core/data/tensorfield/tensorfield.hpp"

#include <string>
#include <algorithm>
#include <stdexcept>


namespace PHARE
{
namespace core
{
    template<std::uint32_t lhs, std::uint32_t rhs>
    NO_DISCARD constexpr std::uint32_t max()
    {
        if constexpr (lhs < rhs)
        {
            return rhs;
        }
        else if constexpr (lhs >= rhs)
        {
            return lhs;
        }
    }


    template<typename T>
    NO_DISCARD std::string to_str(T&& t)
    {
        return t.to_str();
    }



    template<typename Container, typename ContainedT = typename Container::value_type>
    NO_DISCARD bool notIn(ContainedT& obj, Container& list)
    {
        auto sameItem = std::find_if(std::begin(list), std::end(list), [&obj](auto& currentItem) {
            return obj->name() == currentItem->name();
        });

        return sameItem == std::end(list);
    }



} // namespace core
} // namespace PHARE

namespace PHARE::core
{

template<Spannable Span>
void average(Span const& f1, Span const& f2, Span& avg)
{
    auto const size = f1.size();
    auto const d1   = f1.data();
    auto const d2   = f2.data();
    auto av         = avg.data();
    for (std::size_t i = 0; i < size; ++i)
        av[i] = (d1[i] + d2[i]) * .5;
}

// TODO benchmark several ways of converting to primal
// https://github.com/PHAREHUB/PHARE/pull/1486#issuecomment-6401620199
template<std::size_t dim, typename... Ts>
auto& convert_to_fortran_primal(  // DOES NOT WORK ON GHOST BOX!
    Field<dim, Ts...>& dst,       //
    Field<dim, Ts...> const& src, //
    auto const& layout            //
)
{
    bool static constexpr c_ordering = false;

    if (not all(layout.centering(dst), [](auto const c) { return c == QtyCentering::primal; }))
        throw std::runtime_error("Invalid operation, all directions must be primal");

    auto const all_primal
        = all(layout.centering(src), [](auto const c) { return c == QtyCentering::primal; });

    auto const lcl_box = layout.AMRToLocal(layout.AMRBoxFor(dst));
    auto dst_view      = make_array_view<c_ordering>(dst.data(), *lcl_box.shape());

    auto const lcl_zero_box = [&]() {
        auto box = lcl_box;
        box.lower -= lcl_box.lower; // 0
        box.upper -= lcl_box.lower; // shift
        return box;
    }();

    if (all_primal)
        for (auto const [lix, lix0] : boxes_iterator(lcl_box, lcl_zero_box))
            dst_view(lix0) = src(lix);

    else
    {
        RecenteredField const primal_src{src, ConstArray<QtyCentering, dim>(QtyCentering::primal),
                                         layout};
        for (auto const [lix, lix0] : boxes_iterator(lcl_box, lcl_zero_box))
            dst_view(lix0) = primal_src(*lix);
    }

    return dst;
}

template<typename Field_t, typename PQ, std::size_t rank>
auto& convert_to_fortran_primal(               //
    TensorField<Field_t, PQ, rank>& dst,       //
    TensorField<Field_t, PQ, rank> const& src, //
    auto const& layout                         //
)
{
    for (std::size_t ci = 0; ci < src.size(); ++ci)
        convert_to_fortran_primal(dst[ci], src[ci], layout);
    return dst;
}


template<template<typename> typename Op, std::size_t dim, typename PhysicalQuantity,
         typename Data_t>
void operate(Field<dim, PhysicalQuantity, Data_t>& dst,
             Field<dim, PhysicalQuantity, Data_t> const& src, auto&&... args)
{
    using Operator = Op<Data_t>;
    for (std::size_t i = 0; i < dst.size(); ++i)
        Operator{dst.data()[i]}(src.data()[i], args...);
}

template<template<typename> typename Op, typename Field_t, typename PQ, std::size_t rank>
void operate(TensorField<Field_t, PQ, rank>& dst, TensorField<Field_t, PQ, rank> const& src,
             auto&&... args)
{
    using Operator = Op<typename Field_t::type>;
    for (std::size_t ci = 0; ci < dst.size(); ++ci)
        for (std::size_t i = 0; i < dst[ci].size(); ++i)
            Operator{dst[ci].data()[i]}(src[ci].data()[i], args...);
}


} // namespace PHARE::core

#endif // PHARE_CORE_UTILITIES_ALGORITHM_HPP
