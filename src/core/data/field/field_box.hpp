#ifndef PHARE_CORE_DATA_FIELD_FIELD_BOX_HPP
#define PHARE_CORE_DATA_FIELD_FIELD_BOX_HPP

#include "core/utilities/types.hpp"
#include "core/utilities/box/box.hpp"
#include "core/data/field/field_box_span.hpp"
#include "core/data/ndarray/ndarray_vector.hpp"


#include <vector>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <type_traits>

namespace PHARE::core
{

template<typename Field_t_>
class FieldBox
{
    using value_type = std::decay_t<typename Field_t_::value_type>;

public:
    using Field_t = Field_t_;

    auto constexpr static dimension = Field_t::dimension;

    Field_t& field;
    Box<int, dimension> amr_box;
    Box<std::uint32_t, dimension> lcl_box;

    template<typename GridLayout_t>
    FieldBox(Field_t& field_, GridLayout_t const& layout)
        : field{field_}
        , amr_box{layout.AMRGhostBoxFor(field.physicalQuantity())}
        , lcl_box{layout.ghostBoxFor(field)}
    {
    }

    template<typename GridLayout_t>
    FieldBox(Field_t& field_, GridLayout_t const& layout,
             Box<std::uint32_t, dimension> const& selection)
        : field{field_}
        , amr_box{layout.localToAMR(selection)}
        , lcl_box{selection}
    {
    }

    template<typename GridLayout_t>
    FieldBox(Field_t& field_, GridLayout_t const& layout, Box<int, dimension> const& selection)
        : field{field_}
        , amr_box{selection}
        , lcl_box{layout.AMRToLocal(selection)}
    {
    }


    template<typename Operator>
    void set_from(std::vector<value_type> const& vec, std::size_t seek = 0);

    void append_to(std::vector<value_type>& vec);
};


template<typename Operator>
void operate_on_field_row(auto& dst, auto const& src)
    requires(std::is_same_v<Operator, Equals<double>>)
{
    std::memcpy(dst.data(), src.data(), src.size() * sizeof(double));
}

template<typename Operator>
void operate_on_field_row(auto& dst, auto const& src)
{
    for (std::size_t i = 0; i < src.size(); ++i)
        Operator{dst[i]}(src[i]);
}


template<typename Operator>
void operate_on_fields(auto& dst, auto const& src)
{
    assert(dst.lcl_box.size() == src.lcl_box.size());
    auto d_span       = make_field_box_span(dst.lcl_box, dst.field);
    auto const s_span = make_field_box_span(src.lcl_box, src.field);

    auto d_slabs = d_span.begin();
    auto s_slabs = s_span.begin();
    for (; s_slabs != s_span.end(); ++s_slabs, ++d_slabs)
    {
        auto d_rows = d_slabs.begin();
        auto s_rows = s_slabs.begin();
        for (; s_rows != s_slabs.end(); ++s_rows, ++d_rows)
            operate_on_field_row<Operator>(*d_rows, *s_rows);
    }
}

void max_of_fields(auto& dst, auto const& src)
{
    assert(dst.lcl_box.size() == src.lcl_box.size());
    auto src_it = src.lcl_box.begin();
    auto dst_it = dst.lcl_box.begin();
    for (; dst_it != dst.lcl_box.end(); ++src_it, ++dst_it)
    {
        auto& dst_val = dst.field(*dst_it);
        auto& src_val = src.field(*src_it);
        dst_val       = std::max(dst_val, src_val);
    }
}


template<typename Operator>
void set_field_row_from(auto& dst, auto const* src_data)
    requires(std::is_same_v<Operator, Equals<double>>)
{
    std::memcpy(dst.data(), src_data, dst.size() * sizeof(double));
}

template<typename Operator>
void set_field_row_from(auto& dst, auto const* src_data)
{
    for (std::size_t i = 0; i < dst.size(); ++i, ++src_data)
        Operator{dst[i]}(*src_data);
}


template<typename Field_t>
template<typename Operator>
void FieldBox<Field_t>::set_from(std::vector<value_type> const& vec, std::size_t seek)
{
    for (auto& slab : make_field_box_span(lcl_box, field))
        for (auto& row : slab)
        {
            set_field_row_from<Operator>(row, vec.data() + seek);
            seek += row.size();
        }
}


template<typename Field_t>
void FieldBox<Field_t>::append_to(std::vector<value_type>& vec)
{
    using value_type = Field_t::value_type;
    // reserve vec before use!
    std::size_t seek = vec.size();
    vec.resize(vec.size() + lcl_box.size());

    for (auto const& slab : make_field_box_span(lcl_box, field))
        for (auto const& row : slab)
        {
            std::memcpy(vec.data() + seek, row.data(), row.size() * sizeof(value_type));
            seek += row.size();
        }
}



template<typename Field_t>
void fill_ghost(Field_t& field, auto const& layout, auto const v)
{
    std::size_t const nghosts  = layout.nbrGhosts();
    std::size_t const max      = nghosts - 1;
    field[NdArrayMask{0, max}] = v;
}


} // namespace PHARE::core


#endif
