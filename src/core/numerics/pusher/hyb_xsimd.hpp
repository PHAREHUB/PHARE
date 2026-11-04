#ifndef PHARE_CORE_PUSHER_HYB_XSIMD_HPP
#define PHARE_CORE_PUSHER_HYB_XSIMD_HPP

// PHARE_HYB_XSIMD (1 or 0) selects the backend of the hybrid lane kernels at compile time:
// explicit xsimd batches of hyb::laneCount = 8 doubles, which needs xsimd and AVX-512, or the
// lane loops left to the compiler (hyb::AutoVec). PHARE_WITH_XSIMD is set by CMake (option
// withXsimd); without it the xsimd header is used if it is on the include path.
#ifndef PHARE_WITH_XSIMD
#if __has_include(<xsimd/xsimd.hpp>)
#define PHARE_WITH_XSIMD 1
#else
#define PHARE_WITH_XSIMD 0
#endif
#endif

#if PHARE_WITH_XSIMD && defined(__AVX512F__)
#define PHARE_HYB_XSIMD 1
#else
#define PHARE_HYB_XSIMD 0
#endif

#include "core/numerics/pusher/particle_lanes.hpp"

#include <type_traits>

namespace PHARE::core::hyb
{
using Backend = std::conditional_t<PHARE_HYB_XSIMD, Xsimd, AutoVec>;
}

#if PHARE_HYB_XSIMD

// xsimd 13 lane helpers of the hyb::Xsimd backend. One block of hyb::laneCount = 8 particles
// is one batch:
//   doubles: xsimd::batch<double> (default arch, avx512)
//   ints   : xsimd::batch<int32_t, avx2> (8 x int32, the lane storage type)
// int <-> double and the field row loads use AVX-512F intrinsics on the batch registers.
// Only exact operations (abs, compare + select, int <-> double of integral values) replace
// scalar library calls; no xsimd min/max (they differ from std on ties/NaN).

#include <xsimd/xsimd.hpp>
#include <immintrin.h>
#include <cstdint>

namespace PHARE::core::hyb::xs
{

using bd = xsimd::batch<double>;
using bi = xsimd::batch<std::int32_t, xsimd::avx2>;
using bm = xsimd::batch_bool<double>;
static_assert(bd::size == laneCount, "hyb::Xsimd needs 8 double lanes (avx512)");
static_assert(bi::size == laneCount);

[[gnu::always_inline]] inline bd ld(double const* p)
{
    return bd::load_unaligned(p);
}
[[gnu::always_inline]] inline void st(double* p, bd const v)
{
    v.store_unaligned(p);
}
[[gnu::always_inline]] inline bi ldi(int const* p)
{
    return bi::load_unaligned(p);
}
[[gnu::always_inline]] inline void sti(int* p, bi const v)
{
    v.store_unaligned(p);
}
// exact: int32 -> double
[[gnu::always_inline]] inline bd toD(bi const v)
{
    return bd(_mm512_cvtepi32_pd(static_cast<__m256i>(v)));
}
// exact for integral values
[[gnu::always_inline]] inline bi toI(bd const v)
{
    return bi(_mm512_cvttpd_epi32(static_cast<__m512d>(v)));
}
/** cols[k][l] = base[row[l] + k], k < K <= 4 (the rows4 field gather): lane l loads one
 * contiguous row of W = 4 doubles (ymm) at base + row[l], the 8 x W tile is transposed
 * in registers (unpacklo/hi + shuf_f64x2; xsimd 13 has no transpose) and only its first K columns
 * are kept. The values are the ones gather() reads, so the interpolation sums are unchanged.
 *
 * A row reads W - K doubles past the stencil. Ghost nodes do not bound this (the stencil may end
 * at the last node of the field), so a lane whose row passes the end of the buffer (row + W >
 * size) loads only the elements below size, by an AVX-512 masked load (masked-off elements are
 * not accessed, no fault). `edge` (any lane of the stencil past size - W, from its last row)
 * selects that path for the whole call; the usual call takes the unmasked path.
 */
template<std::size_t K>
[[gnu::always_inline]] inline void loadRows(double const* base, bi const row, int const size,
                                            bool const edge, bd (&cols)[K])
{
    constexpr int W = 4;
    static_assert(K >= 1 and K <= 4);
    alignas(32) int r[laneCount];
    sti(r, row);
    __m256d y[laneCount];
    if (!edge) [[likely]]
    {
        for (std::size_t l = 0; l < laneCount; ++l)
            y[l] = _mm256_loadu_pd(base + r[l]);
    }
    else
    {
        for (std::size_t l = 0; l < laneCount; ++l)
        {
            int const n      = size - r[l];
            __mmask8 const m = n >= W ? 0xF : (n <= 0 ? 0 : (1u << n) - 1u);
            y[l]             = _mm256_maskz_loadu_pd(m, base + r[l]);
        }
    }
    // z0 = rows 0|2, z1 = 1|3, z2 = 4|6, z3 = 5|7 (one row per 256-bit half)
    __m512d const z0 = _mm512_insertf64x4(_mm512_castpd256_pd512(y[0]), y[2], 1);
    __m512d const z1 = _mm512_insertf64x4(_mm512_castpd256_pd512(y[1]), y[3], 1);
    __m512d const z2 = _mm512_insertf64x4(_mm512_castpd256_pd512(y[4]), y[6], 1);
    __m512d const z3 = _mm512_insertf64x4(_mm512_castpd256_pd512(y[5]), y[7], 1);
    // 128-bit chunks: t0 = (c0|c2 of lanes 01, c0|c2 of 23), t1 the same for c1|c3,
    // t2/t3 for lanes 45, 67
    __m512d const t0 = _mm512_unpacklo_pd(z0, z1);
    __m512d const t1 = _mm512_unpackhi_pd(z0, z1);
    __m512d const t2 = _mm512_unpacklo_pd(z2, z3);
    __m512d const t3 = _mm512_unpackhi_pd(z2, z3);
    cols[0]          = bd(_mm512_shuffle_f64x2(t0, t2, 0x88));
    if constexpr (K > 1)
        cols[1] = bd(_mm512_shuffle_f64x2(t1, t3, 0x88));
    if constexpr (K > 2)
        cols[2] = bd(_mm512_shuffle_f64x2(t0, t2, 0xDD));
    if constexpr (K > 3)
        cols[3] = bd(_mm512_shuffle_f64x2(t1, t3, 0xDD));
}

// 1. where m, else 0.
[[gnu::always_inline]] inline bd ones(bm const m)
{
    return xsimd::select(m, bd(1.), bd(0.));
}


/** start index and weights of 8 lanes, as Interpolator::localIndexAndWeights_ /
 * computeStartLeftShift + Weighter<order>::computeWeight, same expression order.
 * lc: local cell, d: delta.
 */
// the multiply-adds are explicit fma where GCC contracts the AutoVec lane loops
template<int order, bool dual, std::size_t nbrPts>
[[gnu::always_inline]] inline void startAndWeights(bi const lc, bd const d, bi& s, bd (&w)[nbrPts])
{
    static_assert(order >= 1 and order <= 3, "hyb::Xsimd: interp order 1, 2 or 3");
    if constexpr (order == 1)
    {
        if constexpr (dual)
            s = lc - toI(ones(d < bd(.5)));
        else
            s = lc;
    }
    else if constexpr (order == 3)
    {
        if constexpr (dual)
            s = lc - bi(1) - toI(ones(d < bd(.5)));
        else
            s = lc - bi(1);
    }
    else
    {
        if constexpr (dual)
            s = lc - bi(1);
        else
            s = lc - toI(ones(d < bd(.5)));
    }
    bd normalizedPos = toD(lc) + d;
    if constexpr (dual)
        normalizedPos = normalizedPos - bd(.5);
    if constexpr (order == 1)
    {
        w[1] = normalizedPos - toD(s);
        w[0] = bd(1.) - w[1];
    }
    else if constexpr (order == 3)
    {
        // Weighter<3>; a product with another use (coef_sq) is not contracted
        bd const index    = toD(s) - normalizedPos;
        bd const coef1    = xsimd::fma(bd(0.5), index, bd(1.));
        bd const coef2    = index + bd(1.);
        bd const coef3    = index + bd(2.);
        bd const coef4    = xsimd::fnma(bd(0.5), index + bd(3.), bd(1.));
        bd const coef2_sq = coef2 * coef2;
        bd const coef2_cu = coef2_sq * coef2;
        bd const coef3_sq = coef3 * coef3;
        bd const coef3_cu = coef3_sq * coef3;
        w[0]              = bd(4. / 3.) * coef1 * coef1 * coef1;
        w[1]              = xsimd::fnma(bd(0.5), coef2_cu, bd(2. / 3.) - coef2_sq);
        w[2]              = xsimd::fma(bd(0.5), coef3_cu, bd(2. / 3.) - coef3_sq);
        w[3]              = bd(4. / 3.) * coef4 * coef4 * coef4;
    }
    else
    {
        bd const delta = toD(s + bi(1)) - normalizedPos;
        bd const coef1 = bd(0.5) + delta;
        bd const coef2 = delta;
        bd const coef3 = bd(0.5) - delta;
        w[0]           = bd(0.5) * coef1 * coef1;
        w[1]           = xsimd::fnma(coef2, coef2, bd(0.75)); // contracted in AutoVec too
        w[2]           = bd(0.5) * coef3 * coef3;
    }
}

} // namespace PHARE::core::hyb::xs

#endif // PHARE_HYB_XSIMD

#endif
