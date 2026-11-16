#ifndef PHARE_CORE_NUMERICS_FILTER_BINOMIAL_FILTER_HPP
#define PHARE_CORE_NUMERICS_FILTER_BINOMIAL_FILTER_HPP

#include "core/data/field/field_traits.hpp"
#include "core/data/field/field_scratch.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"

#include <array>
#include <cassert>
#include <cstddef>
#include <algorithm>


namespace PHARE::core
{
/**
 * @brief BinomialFilter applies one pass of the binomial filter (1/4, 1/2, 1/4) in every
 * direction of the grid to the physical nodes of a field, on its own centering:
 *
 *     f(i, j, k) <- sum_{a,b,c = -1,0,1} w_a w_b w_c f(i + a, j + b, k + c),  w = (1/4, 1/2, 1/4)
 *
 * The transfer function of a pass is prod_d cos^2(k_d dl_d / 2): it keeps constants and removes
 * the Nyquist mode. A pass reads one ghost layer in every direction, edges and corners included,
 * and writes the physical nodes only: the ghosts must be filled again before the next pass.
 *
 * The tensor product is done as one sweep per direction. A sweep along d is applied to the
 * physical nodes along d and the directions already swept, and to one more node on each side
 * along the directions still to sweep, so that the last sweep reads filtered ghost rows. The
 * intermediate results live in the field scratch buffer, which grows once and is then reused.
 */
template<typename GridLayout>
class BinomialFilter
{
    static constexpr auto dimension = GridLayout::dimension;
    using Nodes_t                   = std::array<std::size_t, 3>;

public:
    explicit BinomialFilter(GridLayout const& layout)
        : layout_{layout}
    {
    }


    template<IsField Field>
    void operator()(Field& field) const
    {
        // indexes padded to 3 directions at the front, the last one being contiguous
        constexpr auto pad = 3 - dimension;
        Nodes_t nodes{1, 1, 1}, lo{0, 0, 0}, hi{0, 0, 0};
        for (std::size_t d = 0; d < dimension; ++d)
        {
            auto const direction = static_cast<Direction>(d);
            nodes[pad + d]       = field.shape()[d];
            lo[pad + d]          = layout_.physicalStartIndex(field, direction);
            hi[pad + d]          = layout_.physicalEndIndex(field, direction);
            assert(lo[pad + d] >= 1 and hi[pad + d] + 1 < nodes[pad + d]);
        }
        Nodes_t const strides{nodes[1] * nodes[2], nodes[2], 1};

        auto const size = field.size();
        ScratchWrite scratch;
        auto* const buffers = scratch.get().get_no_copy((dimension == 3 ? 2 : 1) * size).data();

        // the last sweep writes into the field: in 1D it is also the first, it reads a copy
        double const* src = field.data();
        if constexpr (dimension == 1)
        {
            std::copy(field.data(), field.data() + size, buffers);
            src = buffers;
        }

        for (std::size_t d = 0; d < dimension; ++d)
        {
            double* const dst = d == dimension - 1 ? field.data() : buffers + (d % 2) * size;

            Nodes_t from = lo, to = hi;
            for (std::size_t e = d + 1; e < dimension; ++e)
            {
                --from[pad + e];
                ++to[pad + e];
            }
            sweep_(src, dst, from, to, strides, strides[pad + d]);
            src = dst;
        }
    }


    // each component of a vector field
    template<typename VecField>
        requires(!IsField<VecField>)
    void operator()(VecField& vecfield) const
    {
        for (std::size_t c = 0; c < 3; ++c)
            (*this)(vecfield[c]);
    }


private:
    // dst = (1/4, 1/2, 1/4) stencil of src along the direction of the given stride, on the nodes
    // from ... to (inclusive)
    static void sweep_(double const* const src, double* const dst, Nodes_t const& from,
                       Nodes_t const& to, Nodes_t const& strides, std::size_t const stride)
    {
        for (auto i = from[0]; i <= to[0]; ++i)
            for (auto j = from[1]; j <= to[1]; ++j)
            {
                auto const row = i * strides[0] + j * strides[1];
                for (auto k = from[2]; k <= to[2]; ++k)
                {
                    auto const n = row + k;
                    dst[n]       = 0.25 * src[n - stride] + 0.5 * src[n] + 0.25 * src[n + stride];
                }
            }
    }


    GridLayout const& layout_;
};

} // namespace PHARE::core

#endif
