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

template<typename Field_t>
class FieldBox
{
    using value_type = std::decay_t<typename Field_t::value_type>;

public:
    auto constexpr static dimension = Field_t::dimension;

    Field_t& field;
    Box<int, dimension> amr_ghost_box;
    Box<std::uint32_t, dimension> lcl_box;

    template<typename GridLayout_t>
    FieldBox(Field_t& field_, GridLayout_t const& layout)
        : field{field_}
        , amr_ghost_box{layout.AMRGhostBoxFor(field.physicalQuantity())}
        , lcl_box{layout.ghostBoxFor(field)}
    {
    }

    template<typename GridLayout_t>
    FieldBox(Field_t& field_, GridLayout_t const& layout,
             Box<std::uint32_t, dimension> const& selection)
        : field{field_}
        , amr_ghost_box{layout.AMRGhostBoxFor(field.physicalQuantity())}
        , lcl_box{selection}
    {
    }

    template<typename GridLayout_t>
    FieldBox(Field_t& field_, GridLayout_t const& layout, Box<int, dimension> const& selection)
        : field{field_}
        , amr_ghost_box{layout.AMRGhostBoxFor(field.physicalQuantity())}
        , lcl_box{layout.AMRToLocal(selection)}
    {
    }


    template<typename Operator>
    void set_from(std::vector<value_type> const& vec, std::size_t seek = 0);

    void append_to(std::vector<value_type>& vec);
};



template<typename Operator, typename T>
void operate_on_span(auto& dst, T const* src_data)
    requires(std::is_same_v<Operator, Equals<T>>)
{
    std::memcpy(dst.data(), src_data, dst.size() * sizeof(T));
}

template<typename Operator, typename T>
void operate_on_span(auto& dst, T const* src_data)
{
    for (std::size_t i = 0; i < dst.size(); ++i, ++src_data)
        Operator{dst[i]}(*src_data);
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
        auto d_spans = d_slabs.begin();
        auto s_spans = s_slabs.begin();
        for (; s_spans != s_slabs.end(); ++s_spans, ++d_spans)
            operate_on_span<Operator>(*d_spans, (*s_spans).data());
    }
}




template<typename Field_t>
template<typename Operator>
void FieldBox<Field_t>::set_from(std::vector<value_type> const& vec, std::size_t seek)
{
    for (auto& slab : make_field_box_span(lcl_box, field))
        for (auto& span : slab)
        {
            operate_on_span<Operator>(span, vec.data() + seek);
            seek += span.size();
        }
}


template<typename Field_t>
void FieldBox<Field_t>::append_to(std::vector<value_type>& vec)
{
    // reserve vec before use!
    std::size_t seek = vec.size();
    vec.resize(vec.size() + lcl_box.size());

    for (auto const& slab : make_field_box_span(lcl_box, field))
        for (auto const& span : slab)
        {
            std::memcpy(vec.data() + seek, span.data(), span.size() * sizeof(value_type));
            seek += span.size();
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
