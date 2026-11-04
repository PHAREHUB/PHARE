#ifndef PHARE_CORE_NUMERICS_INTERPOLATOR_CELL_DEPOSIT_HPP
#define PHARE_CORE_NUMERICS_INTERPOLATOR_CELL_DEPOSIT_HPP

// Cell deposit of sorted domain arrays (xsimd backend only): density and flux deposit of blocks of
// laneCount particles into a lane-private window of the four moments, reduced over the lanes and
// added to the grid once per window flush.
//
// Window: per dim W nodes from a base node, buffer acc[4 moments][W^dim nodes][laneCount lanes].
// A block of lanes (one cell, or neighbouring cells) has the bounding box of its stencils; when
// that box lies in the window the block is accumulated there (per window node and lane, the
// lane's weight product or 0: per-dim weight select on the lane's stencil offset). Otherwise the
// window is flushed (lanes summed, added to the grid over the touched nodes, zeroed) and
// re-anchored on the block's key cell: base = key - maxShift - margin, clamped so that the block
// fits. A block wider than the window is deposited lane by lane on the grid (fallback).
// On cell-sorted particles a window holds one cell (corrector, W = U = union of the stencils of
// one cell's particles) or one pre-push cell with a margin of one cell (predictor, W = U + 2):
// one flush per cell. Unsorted particles flush more, then fall back.
//
// Per lane the products are the scalar deposit's (dN * xw * yw then multiply-add with the last
// weight); only the summation order differs (lane-private partial sums, then the lane sum).

#include "core/numerics/interpolator/interpolator.hpp"
#include "core/numerics/pusher/particle_lanes.hpp"
#include "core/numerics/pusher/hyb_xsimd.hpp"

#if PHARE_HYB_XSIMD

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <vector>

namespace PHARE::core
{

template<std::size_t dim, std::size_t interpOrder>
class CellDeposit
{
    static constexpr std::size_t B = hyb::laneCount;
    static constexpr int nbrPts    = nbrPointsSupport(interpOrder);
    // primal start = cell - leftShift, leftShift in [minShift, maxShift]
    static constexpr int maxShift = interpOrder == 1 ? 0 : 1;
    static constexpr int minShift = interpOrder == 3 ? 1 : 0;
    // union of the stencils of the particles of one cell, per dim
    static constexpr int U    = nbrPts + maxShift - minShift;
    static constexpr int Wmax = 16;
    // cells a particle moves in a push: the predictor window grows by this margin on each side
    static constexpr int predMargin = 1;
    static_assert(U + 2 * predMargin <= Wmax);

    using Interp = Interpolator<dim, interpOrder>;

    struct alignas(64) L8
    {
        double x[B];
    };

    // Lane-axis arithmetic as one xsimd batch of the B lanes (same per-lane operations as the
    // lane loops: product, exact fma). With plain lane loops nested in the 4-moment loop, clang
    // unrolls the lane loop and vectorizes the moment loop instead: gather + scatter on the
    // accumulator at a stride of nodes x lanes (3D order-2 corrector deposit ~7x slower than GCC).
    using bd = hyb::xs::bd;
    static_assert(bd::size == B);
    [[gnu::always_inline]] static inline bd ld_(double const* p) { return bd::load_unaligned(p); }
    // a = fma(x, y, a) on the B lanes
    [[gnu::always_inline]] static inline void fmaLanes_(double* a, bd const x, bd const y)
    {
        xsimd::fma(x, y, bd::load_unaligned(a)).store_unaligned(a);
    }

public:
    using Block = hyb::LaneBlock<dim>;

    /** start a deposit pass on targets (Interpolator::depositTargets); withMargin: the window
     * has a margin of predMargin cells around the key cell (pushed particles keyed on their
     * pre-push cell, or an unsorted array)
     */
    template<typename Targets>
    void begin(Targets const& t, std::array<int, dim> const& amrToLocal, bool const withMargin)
    {
        int const margin = withMargin ? predMargin : 0;
        nd_              = t.density;
        fx_              = t.xFlux;
        fy_              = t.yFlux;
        fz_              = t.zFlux;
        a2l_             = amrToLocal;
        for (std::size_t d = 0; d < dim; ++d)
        {
            shape_[d] = static_cast<int>(t.shape[d]);
            W_[d]     = U + 2 * margin;
            shift_[d] = maxShift + margin;
        }
        nNodes_ = 1;
        for (std::size_t d = 0; d < dim; ++d)
            nNodes_ *= W_[d];
        if (acc_.size() < static_cast<std::size_t>(4 * nNodes_))
            acc_.assign(4 * nNodes_, L8{});
        active_ = false;
    }

    void end() { flush_(); }

    /** deposit the lanes of `in` (members weight, cell (AMR), delta, v as [component][lane],
     * every lane holding valid data) with mask[l] != 0. key: AMR cell the block is sorted on
     */
    template<typename Lanes>
    [[gnu::always_inline]] inline void block(Lanes const& in, int const* mask,
                                             std::array<int, dim> const& key)
    {
        int starts[dim][B];
        double w[dim][nbrPts][B];
        for (std::size_t d = 0; d < dim; ++d)
        {
            int const a2l = a2l_[d];
            for (std::size_t l = 0; l < B; ++l)
            {
                double const dl = in.delta[d][l];
                int const lc    = in.cell[d][l] + a2l;
                int const s
                    = lc
                      - Interp::template computeStartLeftShift<QtyCentering, QtyCentering::primal>(
                          dl);
                double const normalizedPos = lc + dl;
                std::array<double, nbrPts> wl;
                Weighter<interpOrder>{}.computeWeight(normalizedPos, s, wl);
                starts[d][l] = s;
                for (int k = 0; k < nbrPts; ++k)
                    w[d][k][l] = wl[k];
            }
        }

        // bounding box of the stencils of the masked lanes (lane loops: min/max reductions)
        std::array<int, dim> lo, hi;
        int nOn = 0;
        for (std::size_t l = 0; l < B; ++l)
            nOn += mask[l] != 0;
        if (nOn == 0)
            return;
        for (std::size_t d = 0; d < dim; ++d)
        {
            int mn = 1 << 30, mx = -(1 << 30);
            for (std::size_t l = 0; l < B; ++l)
            {
                mn = std::min(mn, mask[l] ? starts[d][l] : (1 << 30));
                mx = std::max(mx, mask[l] ? starts[d][l] : -(1 << 30));
            }
            lo[d] = mn;
            hi[d] = mx + nbrPts;
        }

        if (!(active_ and fits_(base_, lo, hi)))
        {
            std::array<int, dim> nb;
            bool wide = false;
            for (std::size_t d = 0; d < dim; ++d)
            {
                if (hi[d] - lo[d] > W_[d])
                    wide = true;
                int const kb = key[d] + a2l_[d] - shift_[d];
                nb[d]        = std::clamp(kb, hi[d] - W_[d], lo[d]);
            }
            if (wide)
            {
                scalarLanes_(in, mask, starts, w);
                return;
            }
            flush_();
            base_   = nb;
            active_ = true;
            for (std::size_t d = 0; d < dim; ++d)
            {
                tlo_[d] = W_[d];
                thi_[d] = 0;
            }
        }

        // the deposit coefficients, 0 on masked lanes
        double dM[4][B];
        for (std::size_t l = 0; l < B; ++l)
        {
            double const wgt = mask[l] ? in.weight[l] : 0.;
            dM[0][l]         = 1. * wgt;
            dM[1][l]         = in.v[0][l] * wgt;
            dM[2][l]         = in.v[1][l] * wgt;
            dM[3][l]         = in.v[2][l] * wgt;
        }

        // per dim and window row r of the block box: lane weight w[r - offset] or 0
        std::array<int, dim> ext, o;
        double ws[dim][U + 2 * 2 + Wmax][B];
        for (std::size_t d = 0; d < dim; ++d)
        {
            ext[d]  = hi[d] - lo[d];
            o[d]    = lo[d] - base_[d];
            tlo_[d] = std::min(tlo_[d], o[d]);
            thi_[d] = std::max(thi_[d], o[d] + ext[d]);
            // lane offsets as doubles (small integers, exact) so the select is a batch compare:
            // as a lane loop GCC emits one data-dependent branch per lane, row and weight
            // (~1.9 branch mispredictions per particle in the 3D order-2 predictor)
            alignas(64) double off[B];
            for (std::size_t l = 0; l < B; ++l)
                off[l] = mask[l] ? static_cast<double>(starts[d][l] - lo[d]) : 0.;
            bd const offB = ld_(off);
            for (int r = 0; r < ext[d]; ++r)
            {
                bd const k = bd(static_cast<double>(r)) - offB;
                bd v(0.);
                for (int kk = 0; kk < nbrPts; ++kk)
                    v = xsimd::select(k == bd(static_cast<double>(kk)), ld_(w[d][kk]), v);
                v.store_unaligned(ws[d][r]);
            }
        }

        L8* __restrict__ acc = acc_.data();
        auto const nN        = static_cast<std::size_t>(nNodes_);
        if constexpr (dim == 1)
        {
            for (int rx = 0; rx < ext[0]; ++rx)
            {
                auto const node = static_cast<std::size_t>(o[0] + rx);
                bd const wx     = ld_(ws[0][rx]);
                for (std::size_t m = 0; m < 4; ++m)
                    fmaLanes_(acc[m * nN + node].x, ld_(dM[m]), wx);
            }
        }
        else if constexpr (dim == 2)
        {
            for (int rx = 0; rx < ext[0]; ++rx)
            {
                bd const wx = ld_(ws[0][rx]);
                bd px[4];
                for (std::size_t m = 0; m < 4; ++m)
                    px[m] = ld_(dM[m]) * wx;
                auto const row = static_cast<std::size_t>((o[0] + rx) * W_[1] + o[1]);
                for (int ry = 0; ry < ext[1]; ++ry)
                {
                    bd const wy = ld_(ws[1][ry]);
                    for (std::size_t m = 0; m < 4; ++m)
                        fmaLanes_(acc[m * nN + row + ry].x, px[m], wy);
                }
            }
        }
        else
        {
            for (int rx = 0; rx < ext[0]; ++rx)
            {
                bd const wx = ld_(ws[0][rx]);
                bd px[4];
                for (std::size_t m = 0; m < 4; ++m)
                    px[m] = ld_(dM[m]) * wx;
                for (int ry = 0; ry < ext[1]; ++ry)
                {
                    bd const wy = ld_(ws[1][ry]);
                    bd pxy[4];
                    for (std::size_t m = 0; m < 4; ++m)
                        pxy[m] = px[m] * wy;
                    auto const row = static_cast<std::size_t>(
                        ((o[0] + rx) * W_[1] + (o[1] + ry)) * W_[2] + o[2]);
                    for (int rz = 0; rz < ext[2]; ++rz)
                    {
                        bd const wz = ld_(ws[2][rz]);
                        for (std::size_t m = 0; m < 4; ++m)
                            fmaLanes_(acc[m * nN + row + rz].x, pxy[m], wz);
                    }
                }
            }
        }
    }

    /** corrector: all particles of a domain array, by CSR rows when the array is sorted (one
     * row = one cell, blocks of laneCount inside it, last block padded and masked), else by
     * blocks of laneCount in array order
     */
    template<typename Particles, typename Targets>
    [[gnu::noinline]] void depositArray(Particles const& particles, Targets const& t,
                                        std::array<int, dim> const& amrToLocal)
    {
        // unsorted (decayed order): the window gets the predictor's margin of cells
        begin(t, amrToLocal, !particles.csrValid());
        auto const n     = particles.size();
        std::size_t done = 0;
        if (particles.csrValid())
        {
            auto const& csr       = particles.csr();
            int const* const ends = particles.rowEnds();
            int const K           = particles.nbrCells();
            for (int k = 0; k < K; ++k)
            {
                auto const s = static_cast<std::size_t>(csr[k]);
                auto const e = static_cast<std::size_t>(ends[k]);
                if (s == e)
                    continue;
                depositRun_(particles, s, e);
            }
            done = particles.csrPrefix();
        }
        range_(particles, done, n);
        end();
    }

private:
    /** one CSR row [s, e): every particle in the same cell, so the block box is the cell's
     * stencil union (U nodes per dim from cell - maxShift, compile-time), lane-private sums in
     * a local buffer, added to the grid once at the end of the row
     */
    template<typename Particles>
    void depositRun_(Particles const& particles, std::size_t const s, std::size_t const e)
    {
        runCore_(
            particles[s].iCell, e - s,
            [&](std::size_t const i0, std::size_t const n, double (&wgt)[B], double (&dl)[dim][B],
                double (&vl)[3][B]) __attribute__((always_inline)) {
                // lanes: full block straight, last block padded with its last particle, weight 0
                if (n == B)
                    for (std::size_t l = 0; l < B; ++l)
                    {
                        auto const& part = particles[i0 + l];
                        wgt[l]           = part.weight;
                        for (std::size_t d = 0; d < dim; ++d)
                            dl[d][l] = part.delta[d];
                        for (std::size_t c = 0; c < 3; ++c)
                            vl[c][l] = part.v[c];
                    }
                else
                    for (std::size_t l = 0; l < B; ++l)
                    {
                        auto const& part = particles[i0 + std::min(l, n - 1)];
                        wgt[l]           = l < n ? part.weight : 0.;
                        for (std::size_t d = 0; d < dim; ++d)
                            dl[d][l] = part.delta[d];
                        for (std::size_t c = 0; c < 3; ++c)
                            vl[c][l] = part.v[c];
                    }
            },
            s);
    }

    /** one cell's row of `count` particles: loadLanes(i0, n, wgt, dl, vl) fills the lanes of the
     * block of particles [i0, i0 + n) (i0 from first, step B), padded lanes weight 0
     */
    template<typename Cell, typename LoadLanes>
    [[gnu::always_inline]] inline void runCore_(Cell const& cell0, std::size_t const count,
                                                LoadLanes&& loadLanes, std::size_t const first)
    {
        auto const s = first, e = first + count;
        constexpr std::size_t NN            = dim == 1 ? U : dim == 2 ? U * U : U * U * U;
        alignas(64) double accMem[4][NN][B] = {};
#if defined(__clang__)
        // clang promotes the fully unrolled accumulator to registers and, when it exceeds them
        // (2D orders 2 and 3: 64 zmm for 32 registers), shuffles the spilled values between stack slots
        // every block (~3x GCC's zmm stores); an opaque pointer keeps it in memory, as GCC does
        double (*acc)[NN][B] = accMem;
        if constexpr (4 * NN > 16)
            asm("" : "+r"(acc));
#else
        auto& acc = accMem;
#endif
        std::array<int, dim> lo;
        for (std::size_t d = 0; d < dim; ++d)
            lo[d] = cell0[d] + a2l_[d] - maxShift;

        std::array<int, dim> lc; // local cell of the row, shared by all its particles
        for (std::size_t d = 0; d < dim; ++d)
            lc[d] = cell0[d] + a2l_[d];
        for (std::size_t i0 = s; i0 < e; i0 += B)
        {
            std::size_t const n = std::min(B, e - i0);
            alignas(64) double wgt[B], dl[dim][B], vl[3][B];
            loadLanes(i0, n, wgt, dl, vl);
            // per dim and union row r: lane weight w[r - (start - lo)] or 0
            alignas(64) double ws[dim][U][B];
            for (std::size_t d = 0; d < dim; ++d)
            {
                int const c = lc[d];
                for (std::size_t l = 0; l < B; ++l)
                {
                    double const dd = dl[d][l];
                    int const sh    = Interp::template computeStartLeftShift<QtyCentering,
                                                                             QtyCentering::primal>(dd);
                    int const start = c - sh;
                    double const normalizedPos = c + dd;
                    std::array<double, nbrPts> wl;
                    Weighter<interpOrder>{}.computeWeight(normalizedPos, start, wl);
                    int const off = maxShift - sh;
                    for (int r = 0; r < U; ++r)
                    {
                        double v = 0.;
                        for (int k = 0; k < nbrPts; ++k)
                            v = r - off == k ? wl[k] : v;
                        ws[d][r][l] = v;
                    }
                }
            }
            alignas(64) double dM[4][B];
            for (std::size_t l = 0; l < B; ++l)
            {
                dM[0][l] = 1. * wgt[l];
                dM[1][l] = vl[0][l] * wgt[l];
                dM[2][l] = vl[1][l] * wgt[l];
                dM[3][l] = vl[2][l] * wgt[l];
            }
            if constexpr (dim == 1)
            {
                for (int rx = 0; rx < U; ++rx)
                    for (std::size_t m = 0; m < 4; ++m)
                        fmaLanes_(acc[m][rx], ld_(dM[m]), ld_(ws[0][rx]));
            }
            else if constexpr (dim == 2)
            {
                for (int rx = 0; rx < U; ++rx)
                    for (std::size_t m = 0; m < 4; ++m)
                    {
                        bd const px = ld_(dM[m]) * ld_(ws[0][rx]);
                        for (int ry = 0; ry < U; ++ry)
                            fmaLanes_(acc[m][rx * U + ry], px, ld_(ws[1][ry]));
                    }
            }
            else
            {
                for (int rx = 0; rx < U; ++rx)
                    for (int ry = 0; ry < U; ++ry)
                        for (std::size_t m = 0; m < 4; ++m)
                        {
                            bd const pxy = ld_(dM[m]) * ld_(ws[0][rx]) * ld_(ws[1][ry]);
                            for (int rz = 0; rz < U; ++rz)
                                fmaLanes_(acc[m][(rx * U + ry) * U + rz], pxy, ld_(ws[2][rz]));
                        }
            }
        }

        double* const fields[4] = {nd_, fx_, fy_, fz_};
        auto lane_sum           = [](double const* a) {
            return ((a[0] + a[1]) + (a[2] + a[3])) + ((a[4] + a[5]) + (a[6] + a[7]));
        };
        for (std::size_t node = 0; node < NN; ++node)
        {
            std::size_t gi;
            if constexpr (dim == 1)
                gi = static_cast<std::size_t>(lo[0]) + node;
            else if constexpr (dim == 2)
                gi = static_cast<std::size_t>(lo[0] + static_cast<int>(node / U))
                         * static_cast<std::size_t>(shape_[1])
                     + static_cast<std::size_t>(lo[1] + static_cast<int>(node % U));
            else
                gi = (static_cast<std::size_t>(lo[0] + static_cast<int>(node / (U * U)))
                          * static_cast<std::size_t>(shape_[1])
                      + static_cast<std::size_t>(lo[1] + static_cast<int>((node / U) % U)))
                         * static_cast<std::size_t>(shape_[2])
                     + static_cast<std::size_t>(lo[2] + static_cast<int>(node % U));
            for (std::size_t m = 0; m < 4; ++m)
                fields[m][gi] += lane_sum(acc[m][node]);
        }
    }

    template<typename Particles>
    void range_(Particles const& particles, std::size_t const s, std::size_t const e)
    {
        for (std::size_t i0 = s; i0 < e; i0 += B)
        {
            std::size_t const n = std::min(B, e - i0);
            int mask[B];
            for (std::size_t l = 0; l < B; ++l)
            {
                auto const& part = particles[i0 + std::min(l, n - 1)];
                mask[l]          = l < n;
                blk_.weight[l]   = part.weight;
                for (std::size_t d = 0; d < dim; ++d)
                {
                    blk_.cell[d][l]  = part.iCell[d];
                    blk_.delta[d][l] = part.delta[d];
                }
                for (std::size_t c = 0; c < 3; ++c)
                    blk_.v[c][l] = part.v[c];
            }
            std::array<int, dim> key;
            for (std::size_t d = 0; d < dim; ++d)
                key[d] = blk_.cell[d][0];
            block(blk_, mask, key);
        }
    }

    bool fits_(std::array<int, dim> const& base, std::array<int, dim> const& lo,
               std::array<int, dim> const& hi) const
    {
        for (std::size_t d = 0; d < dim; ++d)
            if (lo[d] < base[d] or hi[d] > base[d] + W_[d])
                return false;
        return true;
    }

    // lanes summed (fixed pairwise order), added to the grid on the touched nodes, zeroed
    void flush_()
    {
        if (!active_)
            return;
        active_ = false;
        double* const fields[4] = {nd_, fx_, fy_, fz_};
        auto const nN           = static_cast<std::size_t>(nNodes_);
        auto lane_sum           = [](double* __restrict__ a) {
            double const s = ((a[0] + a[1]) + (a[2] + a[3])) + ((a[4] + a[5]) + (a[6] + a[7]));
            for (std::size_t l = 0; l < B; ++l)
                a[l] = 0.;
            return s;
        };
        auto addNode = [&](std::size_t const node, std::size_t const gi) {
            for (std::size_t m = 0; m < 4; ++m)
                fields[m][gi] += lane_sum(acc_[m * nN + node].x);
        };
        if constexpr (dim == 1)
        {
            for (int x = tlo_[0]; x < thi_[0]; ++x)
                addNode(static_cast<std::size_t>(x), static_cast<std::size_t>(base_[0] + x));
        }
        else if constexpr (dim == 2)
        {
            auto const ny = static_cast<std::size_t>(shape_[1]);
            for (int x = tlo_[0]; x < thi_[0]; ++x)
                for (int y = tlo_[1]; y < thi_[1]; ++y)
                    addNode(static_cast<std::size_t>(x * W_[1] + y),
                            static_cast<std::size_t>(base_[0] + x) * ny
                                + static_cast<std::size_t>(base_[1] + y));
        }
        else
        {
            auto const ny = static_cast<std::size_t>(shape_[1]);
            auto const nz = static_cast<std::size_t>(shape_[2]);
            for (int x = tlo_[0]; x < thi_[0]; ++x)
                for (int y = tlo_[1]; y < thi_[1]; ++y)
                    for (int z = tlo_[2]; z < thi_[2]; ++z)
                        addNode(static_cast<std::size_t>((x * W_[1] + y) * W_[2] + z),
                                (static_cast<std::size_t>(base_[0] + x) * ny
                                 + static_cast<std::size_t>(base_[1] + y))
                                        * nz
                                    + static_cast<std::size_t>(base_[2] + z));
        }
    }

    // fallback: the masked lanes deposited one by one on the grid, depositOne_'s expressions
    template<typename Lanes>
    [[gnu::noinline]] void scalarLanes_(Lanes const& in, int const* mask,
                                        int const (&starts)[dim][B],
                                        double const (&w)[dim][nbrPts][B])
    {
        for (std::size_t l = 0; l < B; ++l)
        {
            if (!mask[l])
                continue;
            double const wgt = in.weight[l];
            double const dN = 1. * wgt, dX = in.v[0][l] * wgt, dY = in.v[1][l] * wgt,
                         dZ = in.v[2][l] * wgt;
            if constexpr (dim == 1)
            {
                for (int ix = 0; ix < nbrPts; ++ix)
                {
                    auto const i = static_cast<std::size_t>(starts[0][l] + ix);
                    nd_[i]       = std::fma(dN, w[0][ix][l], nd_[i]);
                    fx_[i]       = std::fma(dX, w[0][ix][l], fx_[i]);
                    fy_[i]       = std::fma(dY, w[0][ix][l], fy_[i]);
                    fz_[i]       = std::fma(dZ, w[0][ix][l], fz_[i]);
                }
            }
            else if constexpr (dim == 2)
            {
                auto const ny = static_cast<std::size_t>(shape_[1]);
                for (int ix = 0; ix < nbrPts; ++ix)
                    for (int iy = 0; iy < nbrPts; ++iy)
                    {
                        auto const i = static_cast<std::size_t>(starts[0][l] + ix) * ny
                                       + static_cast<std::size_t>(starts[1][l] + iy);
                        double const xw = w[0][ix][l], yw = w[1][iy][l];
                        nd_[i] += dN * xw * yw;
                        fx_[i] += dX * xw * yw;
                        fy_[i] += dY * xw * yw;
                        fz_[i] += dZ * xw * yw;
                    }
            }
            else
            {
                auto const ny = static_cast<std::size_t>(shape_[1]);
                auto const nz = static_cast<std::size_t>(shape_[2]);
                for (int ix = 0; ix < nbrPts; ++ix)
                    for (int iy = 0; iy < nbrPts; ++iy)
                        for (int iz = 0; iz < nbrPts; ++iz)
                        {
                            auto const i = (static_cast<std::size_t>(starts[0][l] + ix) * ny
                                            + static_cast<std::size_t>(starts[1][l] + iy))
                                               * nz
                                           + static_cast<std::size_t>(starts[2][l] + iz);
                            double const xw = w[0][ix][l], yw = w[1][iy][l], zw = w[2][iz][l];
                            nd_[i] += dN * xw * yw * zw;
                            fx_[i] += dX * xw * yw * zw;
                            fy_[i] += dY * xw * yw * zw;
                            fz_[i] += dZ * xw * yw * zw;
                        }
            }
        }
    }

    std::vector<L8> acc_;
    Block blk_;
    double *nd_ = nullptr, *fx_ = nullptr, *fy_ = nullptr, *fz_ = nullptr;
    std::array<int, dim> a2l_{}, shape_{}, W_{}, shift_{}, base_{}, tlo_{}, thi_{};
    int nNodes_ = 0;
    bool active_ = false;
};

} // namespace PHARE::core

#endif // PHARE_HYB_XSIMD

#endif
