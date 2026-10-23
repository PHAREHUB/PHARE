#ifndef CORE_NUMERICS_RECONSTRUCTIONS_STENCIL_HPP
#define CORE_NUMERICS_RECONSTRUCTIONS_STENCIL_HPP

#include "core/utilities/index/index.hpp"

#include <array>
#include <cstddef>
#include <type_traits>

namespace PHARE::core
{
// the 2w+1 values of F along `direction` centred on `index`
template<typename GridLayout, int w, auto direction, typename Field>
auto centred_stencil(Field const& F, MeshIndex<Field::dimension> index)
{
    std::array<std::decay_t<decltype(F(index))>, 2 * w + 1> u;
    auto at = index;
    for (int m = 0; m < w; ++m)
        at = GridLayout::template previous<direction>(at);
    for (int m = 0; m < 2 * w + 1; ++m)
    {
        u[m] = F(at);
        at   = GridLayout::template next<direction>(at);
    }
    return u;
}

// Reconstruction::faces_of for n consecutive cells of a contiguous line, the stencil running
// with `stride` through `u`: same values as cell_faces, in a loop the compiler can vectorise
template<typename Reconstruction>
void cell_faces_line(double const* u, std::ptrdiff_t const stride, int const n, double* lower,
                     double* upper)
{
    constexpr int w = Reconstruction::nghosts - 1;
    for (int k = 0; k < n; ++k)
    {
        std::array<double, 2 * w + 1> s;
        for (int m = 0; m < 2 * w + 1; ++m)
            s[m] = u[k + (m - w) * stride];
        auto const [l, up] = Reconstruction::faces_of(s);
        lower[k]           = l;
        upper[k]           = up;
    }
}
} // namespace PHARE::core

#endif
