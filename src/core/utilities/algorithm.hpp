#ifndef PHARE_CORE_UTILITIES_ALGORITHM_HPP
#define PHARE_CORE_UTILITIES_ALGORITHM_HPP


#include "core/def.hpp"

#include "core/utilities/span.hpp"
#include "core/utilities/types.hpp"
#include "core/data/field/field.hpp"
#include "core/utilities/box/box.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/data/ndarray/ndarray_vector.hpp"
#include "core/data/tensorfield/tensorfield.hpp"


#include "core/data/grid/grid.hpp"
#include "core/data/grid/grid_tiles.hpp"



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


template<typename PhysicalQuantity>
auto convert_to_primal(        //
    auto const& src,           //
    auto const& layout,        //
    auto const lix,            //
    PhysicalQuantity const qty //
)
{
    using PQ         = PhysicalQuantity;
    using GridLayout = std::decay_t<decltype(layout)>;

    if (qty == PQ::Bx)
        return GridLayout::template project<GridLayout::BxToMoments>(src, lix);
    else if (qty == PQ::By)
        return GridLayout::template project<GridLayout::ByToMoments>(src, lix);
    else if (qty == PQ::Bz)
        return GridLayout::template project<GridLayout::BzToMoments>(src, lix);


    // ONLY EVER Scalar type!
    if constexpr (std::is_same_v<PhysicalQuantity, HybridQuantity::Scalar>)
    {
        if (qty == PQ::Ex)
            return GridLayout::template project<GridLayout::ExToMoments>(src, lix);
        else if (qty == PQ::Ey)
            return GridLayout::template project<GridLayout::EyToMoments>(src, lix);
        else if (qty == PQ::Ez)
            return GridLayout::template project<GridLayout::EzToMoments>(src, lix);
    }
    if constexpr (std::is_same_v<PhysicalQuantity, MHDQuantity::Scalar>)
    {
        // if we are not the magnetic field, then all scalars and vectors are cell-centered in MHD
        return GridLayout::template project<GridLayout::cellCenterToFullPrimal>(src, lix);
    }

    throw std::runtime_error("Quantity not supported for conversion to primal.");
}

template<std::size_t dim, typename pq, typename data_t, auto alloc_mode>
auto& convert_to_fortran_primal(                   // DOES NOT WORK ON GHOST BOX!
    Field<dim, pq, data_t, alloc_mode>& dst,       //
    Field<dim, pq, data_t, alloc_mode> const& src, //
    auto const& layout                             //
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
        for (auto const [lix, lix0] : boxes_iterator(lcl_box, lcl_zero_box))
            dst_view(lix0) = convert_to_primal(src, layout, lix, src.physicalQuantity());

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

template<typename GL, typename ND, typename PQ>
void transform(FieldTileSet<GL, ND, PQ> const& in, FieldTileSet<GL, ND, PQ>& out, auto const fn)
{
    std::size_t const size = in().size();
    for (std::size_t i = 0; i < size; ++i)
        std::transform(in[i]().begin(), in[i]().end(), out[i]().begin(), fn);
}


template<typename GL, typename ND, typename PQ>
void transform(FieldTileSet<GL, ND, PQ> const& in0, FieldTileSet<GL, ND, PQ> const& in1,
               FieldTileSet<GL, ND, PQ>& out, auto const fn)
{
    assert(in0.isUsable());
    assert(in1.isUsable());
    assert(out.isUsable());
    auto const size = in0().size();

    for (std::size_t i = 0; i < size; ++i)
    {
        auto& i0 = in0[i]();
        auto& i1 = in1[i]();
        auto& o0 = out[i]();
        std::transform(in0[i]().begin(), in0[i]().end(), in1[i]().begin(), out[i]().begin(), fn);
    }
}


template<typename GL, typename ND, typename PQ>
void average(FieldTileSet<GL, ND, PQ> const& f1, FieldTileSet<GL, ND, PQ> const& f2,
             FieldTileSet<GL, ND, PQ>& avg)
{
    core::transform(f1, f2, avg, std::plus<double>());
    core::transform(avg, avg, [](double x) { return x * 0.5; });
}


template<auto opts>
void transform(basic::Field<opts> const& in, basic::Field<opts>& out, auto const fn)
{
    std::transform(in.begin(), in.end(), out.begin(), fn);
}

template<auto opts>
void transform(basic::Field<opts> const& in0, basic::Field<opts> const& in1,
               basic::Field<opts>& out, auto const fn)
{
    std::transform(in0.begin(), in0.end(), in1.begin(), out.begin(), fn);
}


template<auto opts>
void average(basic::Field<opts> const& f1, basic::Field<opts> const& f2, basic::Field<opts>& avg)
{
    core::transform(f1, f2, avg, std::plus<double>());
    core::transform(avg, avg, [](double x) { return x * 0.5; });
}

template<typename Field_t, std::size_t rank>
void average(basic::TensorField<Field_t, rank> const& vf1,
             basic::TensorField<Field_t, rank> const& vf2, basic::TensorField<Field_t, rank>& Vavg)
{
    auto constexpr static N = detail::tensor_field_dim_from_rank<rank>();

    for (std::size_t i = 0; i < N; ++i)
        average(vf1[i], vf2[i], Vavg[i]);
}


template<typename Field_t>
void accumulate_field(Field_t& dst, Field_t const& src, auto const coef)
    requires(not is_field_tile_set_v<Field_t>)
{
    if (dst.size() != src.size())
        throw std::runtime_error("Cannot accumulate fields with different sizes");

    for (std::size_t i = 0; i < dst.size(); ++i)
        dst.data()[i] = src.data()[i] * coef;
}


template<typename Field_t>
void accumulate_field(Field_t& dst, Field_t const& src, auto const coef)
    requires(is_field_tile_set_v<Field_t>)
{
    auto& dst_tiles = dst();
    auto& src_tiles = src();

    if (dst_tiles.size() != src_tiles.size())
        throw std::runtime_error("Cannot accumulate fields with different number of tiles");

    for (std::size_t i = 0; i < dst_tiles.size(); ++i)
        accumulate_field(dst_tiles[i](), src_tiles[i](), coef);
}



template<typename Field_t>
void accumulate(basic::TensorField<Field_t, 1>& dst, basic::TensorField<Field_t, 1> const& src,
                auto const coeff)
{
    for (std::size_t c = 0; c < dst.size(); ++c)
        accumulate_field(dst[c], src[c], coeff);
}

template<template<typename> typename Op, typename Field_t>
void operate_field(Field_t& dst, Field_t const& src, auto&&... args)
    requires(not is_field_tile_set_v<Field_t>)
{
    using Operator = Op<typename Field_t::type>;
    for (std::size_t i = 0; i < dst.size(); ++i)
        Operator{dst.data()[i]}(src.data()[i], args...);
}

template<template<typename> typename Op, typename Field_t>
void operate_field(Field_t& dst, Field_t const& src, auto&&... args)
    requires(is_field_tile_set_v<Field_t>)
{
    auto& dst_tiles = dst();
    auto& src_tiles = src();
    if (dst_tiles.size() != src_tiles.size())
        throw std::runtime_error("operate_field: tile count mismatch");
    for (std::size_t i = 0; i < dst_tiles.size(); ++i)
        operate_field<Op>(dst_tiles[i](), src_tiles[i](), args...);
}

template<template<typename> typename Op, std::size_t dim, typename PhysicalQuantity, typename Data_t>
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
    for (std::size_t ci = 0; ci < dst.size(); ++ci)
        operate_field<Op>(dst[ci], src[ci], args...);
}


} // namespace PHARE::core

#endif // PHARE_CORE_UTILITIES_ALGORITHM_HPP
