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

/** @brief Pairs a field with the same region of interest expressed in both index
 * spaces: amr_box (global AMR indices) and lcl_box (local, ghost-box-relative
 * storage indices), so code that needs to walk the field data while also
 * knowing where each cell sits on the AMR grid doesn't have to convert back
 * and forth (see operator_across_adjacent_levels() below for such a use).
 *
 * With only a field and layout, the region defaults to the field's whole ghost
 * box; passing a local or AMR box restricts it to that sub-region instead (the
 * other box is then derived from it via the layout).
 */
template<typename Field_t_>
class FieldBox
{
public:
    using Field_t    = Field_t_;
    using value_type = std::decay_t<typename Field_t::value_type>;

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

    // selection given in local (ghost-box-relative) indices
    template<typename GridLayout_t>
    FieldBox(Field_t& field_, GridLayout_t const& layout,
             Box<std::uint32_t, dimension> const& selection)
        : field{field_}
        , amr_box{layout.localToAMR(selection)}
        , lcl_box{selection}
    {
    }

    // selection given in AMR (global) indices
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



/** @brief Bundles, for a coarse/fine FieldBox pair, the two kinds of span iteration
 * operator_across_adjacent_levels() needs to walk in lockstep:
 *  - the "field" spans, which walk local field storage and yield (data span, local point)
 *  - the "amr" spans, which walk the same box in AMR (global) index space and yield
 *    (amr point, span size), used only to know where each span sits on the AMR grid
 *    (so operator_across_adjacent_levels can tell, via modulo ratio, which fine cells
 *    sit on top of a coarse one and which fall in between).
 */
template<typename Field_t>
struct CrossLevelIndices
{
    auto constexpr static dim = Field_t::dimension;
    using FieldPointSpans_t   = FieldBoxPointSpans<Field_t>;
    using FieldPointSpans_ct  = FieldBoxPointSpans<Field_t const>;

    FieldBox<Field_t>& fine;
    FieldBox<Field_t const> const& coarse;

    FieldBoxSpan<Field_t, FieldPointSpans_t> fine_field_span
        = make_field_box_point_span(fine.lcl_box, fine.field);
    FieldBoxSpan<Field_t const, FieldPointSpans_ct> coarse_field_span
        = make_field_box_point_span(coarse.lcl_box, coarse.field);

    BoxSpan<int, dim> fine_amr_span   = make_box_span(fine.amr_box);
    BoxSpan<int, dim> coarse_amr_span = make_box_span(coarse.amr_box);
};



/** @brief Applies fn to every fine cell of `fine` overlapping `coarse`, walking both
 * fields span-by-span (see core/utilities/box/box_span.hpp) instead of point-by-point.
 *
 * `fine` and `coarse` are iterated in lockstep on the fine grid; `coarse`'s iterators
 * only advance every ratio-th fine step (skipped via the two "% ratio" checks below),
 * since a single coarse cell/span underlies `ratio` fine ones (ratio is 2 for the usual
 * AMR refinement ratio, see amr::refinementRatio at call sites).
 *
 * Assumes `coarse` is exactly the box `fine` coarsens to (callers build both FieldBoxes
 * from the same overlap, one coarsened, so this isn't re-checked here).
 */
template<typename Field_t>
void operator_across_adjacent_levels(FieldBox<Field_t const>& coarse, FieldBox<Field_t>& fine,
                                     std::size_t ratio, auto& fn)
{
    auto constexpr static FASTEST_DIR = Field_t::dimension - 1;

    CrossLevelIndices<Field_t> indices{fine, coarse};

    auto fine_field_slabs   = indices.fine_field_span.begin();
    auto fine_amr_slabs     = indices.fine_amr_span.begin();
    auto coarse_field_slabs = indices.coarse_field_span.begin();
    auto coarse_amr_slabs   = indices.coarse_amr_span.begin();

    auto slab_idx = indices.fine_amr_span.slab_begin();

    for (; fine_field_slabs != indices.fine_field_span.end();
         ++fine_field_slabs, ++fine_amr_slabs, ++slab_idx)
    {
        auto fine_field_spans   = fine_field_slabs.begin();
        auto coarse_field_spans = coarse_field_slabs.begin();
        auto fine_amr_spans     = fine_amr_slabs.begin();
        auto coarse_amr_spans   = coarse_amr_slabs.begin();

        auto span_idx = fine_amr_slabs.span_begin();
        for (; fine_field_spans != fine_field_slabs.end();
             ++fine_field_spans, ++fine_amr_spans, ++span_idx)
        {
            auto&& [fine_amr_point, fine_amr_span_size]     = *fine_amr_spans;
            auto&& [coarse_amr_point, coarse_amr_span_size] = *coarse_amr_spans;

            auto&& [fine_span, fine_lcl_point]     = *fine_field_spans;
            auto&& [coarse_span, coarse_lcl_point] = *coarse_field_spans;

            std::size_t fine_idx   = 0;
            std::size_t coarse_idx = 0;
            for (; fine_idx < fine_span.size(); ++fine_idx)
            {
                fn(coarse.field, fine.field, fine_amr_point, coarse_amr_point, fine_lcl_point,
                   coarse_lcl_point, fine_span[fine_idx], coarse_span[coarse_idx]);

                if (fine_amr_point[FASTEST_DIR] % ratio != 0)
                {
                    ++coarse_idx;
                    ++coarse_lcl_point[FASTEST_DIR];
                    ++coarse_amr_point[FASTEST_DIR];
                }

                ++fine_amr_point[FASTEST_DIR];
                ++fine_lcl_point[FASTEST_DIR];
            }

            if (span_idx % ratio != 0)
            {
                ++coarse_field_spans;
                ++coarse_amr_spans;
            }
        }

        if (slab_idx % ratio != 0)
        {
            ++coarse_field_slabs;
            ++coarse_amr_slabs;
        }
    }
}



// plain assignment (Operator == SetEqual<T>) is just a bulk copy, so it can
// skip the element-by-element loop below and memcpy the whole span instead
template<typename Operator, typename T>
void operate_on_span(auto& dst, T const* src_data)
    requires(std::is_same_v<Operator, SetEqual<T>>)
{
    std::memcpy(dst.data(), src_data, dst.size() * sizeof(T));
}

// generic per-element path for any other Operator (e.g. PlusEquals, SetMax, see types.hpp)
template<typename Operator, typename T>
void operate_on_span(auto& dst, T const* src_data)
{
    for (std::size_t i = 0; i < dst.size(); ++i, ++src_data)
        Operator{dst[i]}(*src_data);
}



// applies Operator span-by-span between two same-shaped FieldBoxes (see
// core/utilities/box/box_span.hpp); used to add/copy/max one field's local
// region into another's, e.g. periodic boundary data exchange
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



} // namespace PHARE::core


#endif
