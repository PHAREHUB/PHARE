#ifndef PHARE_CORE_DATA_FIELD_FIELD_BOX_HPP
#define PHARE_CORE_DATA_FIELD_FIELD_BOX_HPP

#include "core/utilities/types.hpp"
#include "core/utilities/box/box.hpp"
#include "core/utilities/box/box_spans.hpp"

#include <vector>
#include <cassert>
#include <cstddef>
#include <type_traits>

namespace PHARE::core
{

template<typename Field_t>
class FieldBox
{
    using value_type = std::decay_t<typename Field_t::value_type>;

    // set_from/append_to/operate_on_fields walk lcl_box span by span (box_spans.hpp)
    // and a span is contiguous in memory only for C ordering
    static_assert(Field_t::is_c_ordered);

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


// Operator{dst[i]}(src[i]) for i in [0, size) (Operator: SetEqual, PlusEquals, SetMax, see
//  types.hpp). A plain loop on purpose, which the compiler vectorizes: memcpy is an
//  out-of-line call that costs more than the copy on the 2-cell spans of thin ghost layers.
template<typename Operator>
void operate_on_span(auto* dst, auto const* src, std::size_t const size)
{
    for (std::size_t i = 0; i < size; ++i)
        Operator{dst[i]}(src[i]);
}


// dst.field = Operator(dst.field, src.field) over two same-shaped boxes, e.g. to copy or
//  add one patch's region into another's
template<typename Operator>
void operate_on_fields(auto& dst, auto const& src)
{
    assert(dst.lcl_box.shape() == src.lcl_box.shape());

    for_each_span(dst.lcl_box, src.lcl_box,
                  [&](auto const& dst_start, auto const& src_start, auto const size) {
                      operate_on_span<Operator>(&dst.field(dst_start), &src.field(src_start), size);
                  });
}


template<typename Field_t>
template<typename Operator>
void FieldBox<Field_t>::set_from(std::vector<value_type> const& vec, std::size_t seek)
{
    for_each_span(lcl_box, [&](auto const& start, auto const size) {
        operate_on_span<Operator>(&field(start), vec.data() + seek, size);
        seek += size;
    });
}


// reserve vec beforehand, or the resize below may reallocate it on every call
template<typename Field_t>
void FieldBox<Field_t>::append_to(std::vector<value_type>& vec)
{
    std::size_t seek = vec.size();
    vec.resize(seek + lcl_box.size());

    for_each_span(lcl_box, [&](auto const& start, auto const size) {
        operate_on_span<SetEqual<value_type>>(vec.data() + seek, &field(start), size);
        seek += size;
    });
}

} // namespace PHARE::core


#endif
