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
    Box<int, dimension> amr_box;
    Box<std::uint32_t, dimension> lcl_box;

    template<typename GridLayout_t>
    FieldBox(Field_t& field_, GridLayout_t const& layout)
        : field{field_}
        , amr_box{layout.AMRBox()}
        , lcl_box{layout.ghostBoxFor(field)}
    {
    }

    template<typename GridLayout_t>
    FieldBox(Field_t& field_, GridLayout_t const& layout,
             Box<std::uint32_t, dimension> const& selection)
        : field{field_}
        , amr_box{layout.AMRBox()}
        , lcl_box{selection}
    {
        // assert(amr_box * selection);
    }

    template<typename GridLayout_t>
    FieldBox(Field_t& field_, GridLayout_t const& layout, Box<int, dimension> const& selection)
        : field{field_}
        , amr_box{layout.AMRBox()}
        , lcl_box{layout.AMRToLocal(selection)}
    {
        // assert(amr_box * selection);
    }


    template<typename Operator, typename Field_t0>
    void op(FieldBox<Field_t0> const& that);

    template<typename Operator>
    void op(std::vector<value_type> const& vec, std::size_t seek = 0);

    void append_to(std::vector<value_type>& vec);
};


template<typename Field_t>
template<typename Operator, typename Field_t0>
void FieldBox<Field_t>::op(FieldBox<Field_t0> const& that)
{
    auto src_it = that.lcl_box.begin();
    auto dst_it = lcl_box.begin();
    for (; dst_it != lcl_box.end(); ++src_it, ++dst_it)
        Operator{field(*dst_it)}(that.field(*src_it));
}



template<typename Field_t>
template<typename Operator>
void FieldBox<Field_t>::op(std::vector<value_type> const& vec, std::size_t seek)
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
