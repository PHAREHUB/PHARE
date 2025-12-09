#ifndef PHARE_CORE_DATA_FIELD_FIELD_BOX_HPP
#define PHARE_CORE_DATA_FIELD_FIELD_BOX_HPP

#include "core/def.hpp"
#include "core/logger.hpp"
#include "core/utilities/types.hpp"
#include "core/utilities/box/box.hpp"

#include <vector>
#include <cstddef>
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


template<typename Operator>
void operate_on_fields(auto& dst, auto const& src)
{
    assert(dst.lcl_box.size() == src.lcl_box.size());
    auto src_it = src.lcl_box.begin();
    auto dst_it = dst.lcl_box.begin();
    for (; dst_it != dst.lcl_box.end(); ++src_it, ++dst_it)
        Operator{dst.field(*dst_it)}(src.field(*src_it));
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



template<typename Field_t>
template<typename Operator>
void FieldBox<Field_t>::set_from(std::vector<value_type> const& vec, std::size_t seek)
{
    auto dst_it = lcl_box.begin();
    for (; dst_it != lcl_box.end(); ++seek, ++dst_it)
        Operator{field(*dst_it)}(vec[seek]);
}

template<typename Field_t>
void FieldBox<Field_t>::append_to(std::vector<value_type>& vec)
{
    // reserve vec before use!
    auto src_it = lcl_box.begin();
    for (; src_it != lcl_box.end(); ++src_it)
        vec.push_back(field(*src_it));
}

} // namespace PHARE::core


#endif
