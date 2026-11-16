#ifndef PHARE_CORE_NUMERICS_ESIRKEPOV_ESIRKEPOV_HPP
#define PHARE_CORE_NUMERICS_ESIRKEPOV_ESIRKEPOV_HPP

#include "core/def.hpp"
#include "core/errors.hpp"
#include "core/utilities/point/point.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/data/vecfield/vecfield_component.hpp"
#include "core/numerics/interpolator/interpolator.hpp"

#include <array>
#include <cmath>
#include <string>
#include <cassert>
#include <cstddef>
#include <algorithm>


namespace PHARE::core
{
/**
 * @brief Esirkepov deposits a charge conserving current (Esirkepov 2001, Comput. Phys. Commun.
 * 135, 144).
 *
 * The current of a particle is the one of its straight motion from x - back*(dt/2)*v to
 * x + front*(dt/2)*v, during (back + front)*dt/2. With back = front = 1 and particles at time n,
 * this is the current at time n of the leapfrog, between the half time positions.
 * The charge density deposited with depositCharge at the two ends of the segment then verifies
 * the discrete continuity equation to round-off, with the same shape functions as the moments.
 *
 * Along each direction d, S0 and S1 are the shapes at the two ends of the segment and
 * dS = S1 - S0. Along the segment, the shapes are taken linear in time, so that their averages
 * are <S> = S0 + dS/2 and <Sa Sb> = S0a S0b + (dSa S0b + S0a dSb)/2 + dSa dSb/3.
 * - the current along a direction d of the grid is dual in d and primal in the others,
 *   J_d(i+1/2) - J_d(i-1/2) = -q w dl_d/duration dS_d <product of the shapes of the others>,
 *   the sum of these differences over the directions is S1x S1y S1z - S0x S0y S0z: continuity.
 * - the current along a direction the grid does not resolve is primal,
 *   J_c = q w v_c <product of the shapes of all directions>,
 *   which is the current of the next dimension summed along that direction.
 * This is the Yee centering of J.
 */
template<std::size_t dim, std::size_t interpOrder>
class Esirkepov
{
    static_assert(dim >= 1 and dim <= 3);

    static constexpr auto support = nbrPointsSupport(interpOrder);
    // a displacement of less than a cell moves the start of the shape by one node at most
    static constexpr auto maxNodes = support + 1;
    using Weights_t                = std::array<double, support>;

    struct Shape
    {
        int start;
        Weights_t weights;
    };

    // shapes at both ends of the segment along one direction, on nodes lo ... lo + nbrNodes - 1
    struct Segment
    {
        int lo;
        int nbrNodes;
        std::array<double, maxNodes> S0{}; // shape at the start
        std::array<double, maxNodes> dS{}; // shape at the end minus shape at the start
    };

public:
    static constexpr auto dimension    = dim;
    static constexpr auto interp_order = interpOrder;


    template<typename Particles, typename VecField, typename GridLayout>
    void operator()(Particles const& particles, VecField& J, GridLayout const& layout,
                    double const charge, double const dt, double const back = 1.,
                    double const front = 1.)
    {
        // the displacement guard below assumes the segment lasts at most dt
        assert(back + front <= 2.);

        auto const& dl      = layout.meshSize();
        auto const halfDt   = 0.5 * dt;
        auto const duration = (back + front) * halfDt;

        for (auto const& particle : particles)
        {
            for (std::size_t d = 0; d < dim; ++d)
                checkDisplacement_(particle.v[d], dt, dl[d], d);

            auto const local = layout.AMRToLocal(Point{particle.iCell});
            std::array<Segment, dim> segments;
            for (std::size_t d = 0; d < dim; ++d)
            {
                auto const pos = static_cast<double>(local[d]) + particle.delta[d];
                auto const dx  = halfDt * particle.v[d] / dl[d];
                segments[d]    = segment_(pos - back * dx, pos + front * dx);

                assert(segments[d].lo >= 0
                       and static_cast<std::size_t>(segments[d].lo + segments[d].nbrNodes)
                               <= J[d].shape()[d] + 1);
            }

            auto const qw = charge * particle.weight;
            std::array<double, dim> coefs; // J_d jump per unit dS_d <...>
            for (std::size_t d = 0; d < dim; ++d)
                coefs[d] = -qw * dl[d] / duration;

            deposit_(segments, coefs, qw, particle.v, J);
        }
    }



    /**
     * @brief depositCharge deposits charge*weight at the position x + timeShift*v of each
     * particle, the particle itself is not modified.
     */
    template<typename Particles, typename Field, typename GridLayout>
    void depositCharge(Particles const& particles, Field& rho, GridLayout const& layout,
                       double const charge, double const timeShift)
    {
        auto const& dl = layout.meshSize();

        for (auto const& particle : particles)
        {
            auto const local = layout.AMRToLocal(Point{particle.iCell});
            std::array<Shape, dim> s;
            for (std::size_t d = 0; d < dim; ++d)
            {
                s[d] = shape_(static_cast<double>(local[d]) + particle.delta[d]
                              + timeShift * particle.v[d] / dl[d]);
                assert(s[d].start >= 0
                       and static_cast<std::size_t>(s[d].start) + support <= rho.shape()[d]);
            }

            // same order of products as the Interpolator density deposit
            auto const qw = charge * particle.weight;
            for (std::size_t ix = 0; ix < support; ++ix)
            {
                auto const xw = qw * s[0].weights[ix];
                if constexpr (dim == 1)
                    rho(s[0].start + ix) += xw;
                else
                    for (std::size_t iy = 0; iy < support; ++iy)
                    {
                        auto const xyw = xw * s[1].weights[iy];
                        if constexpr (dim == 2)
                            rho(s[0].start + ix, s[1].start + iy) += xyw;
                        else
                            for (std::size_t iz = 0; iz < support; ++iz)
                                rho(s[0].start + ix, s[1].start + iy, s[2].start + iz)
                                    += xyw * s[2].weights[iz];
                    }
            }
        }
    }



private:
    template<typename VecField>
    static void deposit_(std::array<Segment, 1> const& segments, std::array<double, 1> const& coefs,
                         double const qw, std::array<double, 3> const& v, VecField& J)
    {
        auto& Jx        = J(Component::X);
        auto& Jy        = J(Component::Y);
        auto& Jz        = J(Component::Z);
        auto const& sx  = segments[0];
        auto const qwvy = qw * v[1];
        auto const qwvz = qw * v[2];

        double runningJx = 0.;
        for (int i = 0; i < sx.nbrNodes; ++i)
        {
            auto const x    = sx.lo + i;
            auto const avgX = average_(sx, i);

            Jy(x) += qwvy * avgX;
            Jz(x) += qwvz * avgX;

            // the last dual node of the range carries no current, the shapes have equal sums
            if (i < sx.nbrNodes - 1)
            {
                runningJx += coefs[0] * sx.dS[i];
                Jx(x) += runningJx;
            }
        }
    }


    template<typename VecField>
    static void deposit_(std::array<Segment, 2> const& segments, std::array<double, 2> const& coefs,
                         double const qw, std::array<double, 3> const& v, VecField& J)
    {
        auto& Jx             = J(Component::X);
        auto& Jy             = J(Component::Y);
        auto& Jz             = J(Component::Z);
        auto const& [sx, sy] = segments;
        auto const& [cx, cy] = coefs;
        auto const qwvz      = qw * v[2];

        // running sums of Jx along x, one per y node
        std::array<double, maxNodes> runningJx{};

        for (int i = 0; i < sx.nbrNodes; ++i)
        {
            auto const x    = sx.lo + i;
            auto const avgX = average_(sx, i);
            // Jx has no node at the last x of the range
            auto* const jx = i < sx.nbrNodes - 1 ? &Jx(x, sy.lo) : nullptr;
            auto* const jy = &Jy(x, sy.lo);
            auto* const jz = &Jz(x, sy.lo);

            double runningJy = 0.;
            for (int j = 0; j < sy.nbrNodes; ++j)
            {
                jz[j] += qwvz * average_(sx, i, sy, j);

                if (i < sx.nbrNodes - 1)
                {
                    runningJx[j] += cx * sx.dS[i] * average_(sy, j);
                    jx[j] += runningJx[j];
                }
                if (j < sy.nbrNodes - 1)
                {
                    runningJy += cy * sy.dS[j] * avgX;
                    jy[j] += runningJy;
                }
            }
        }
    }


    template<typename VecField>
    static void deposit_(std::array<Segment, 3> const& segments, std::array<double, 3> const& coefs,
                         double const /*qw*/, std::array<double, 3> const& /*v*/, VecField& J)
    {
        auto& Jx                 = J(Component::X);
        auto& Jy                 = J(Component::Y);
        auto& Jz                 = J(Component::Z);
        auto const& [sx, sy, sz] = segments;
        auto const& [cx, cy, cz] = coefs;

        // running sums of Jx along x, one per (y, z) node, and of Jy along y, one per z node
        std::array<std::array<double, maxNodes>, maxNodes> runningJx{};

        for (int i = 0; i < sx.nbrNodes; ++i)
        {
            auto const x = sx.lo + i;
            std::array<double, maxNodes> runningJy{};

            for (int j = 0; j < sy.nbrNodes; ++j)
            {
                auto const y  = sy.lo + j;
                auto const wz = cz * average_(sx, i, sy, j);
                // Jx (Jy) has no node at the last x (y) of the range
                auto* const jx = i < sx.nbrNodes - 1 ? &Jx(x, y, sz.lo) : nullptr;
                auto* const jy = j < sy.nbrNodes - 1 ? &Jy(x, y, sz.lo) : nullptr;
                auto* const jz = &Jz(x, y, sz.lo);

                double runningJz = 0.;
                for (int k = 0; k < sz.nbrNodes; ++k)
                {
                    if (i < sx.nbrNodes - 1)
                    {
                        runningJx[j][k] += cx * sx.dS[i] * average_(sy, j, sz, k);
                        jx[k] += runningJx[j][k];
                    }
                    if (j < sy.nbrNodes - 1)
                    {
                        runningJy[k] += cy * sy.dS[j] * average_(sx, i, sz, k);
                        jy[k] += runningJy[k];
                    }
                    if (k < sz.nbrNodes - 1)
                    {
                        runningJz += wz * sz.dS[k];
                        jz[k] += runningJz;
                    }
                }
            }
        }
    }


    // <S> over the segment at node i
    static double average_(Segment const& s, int const i) { return s.S0[i] + 0.5 * s.dS[i]; }

    // <Sa Sb> over the segment at node (i, j)
    static double average_(Segment const& a, int const i, Segment const& b, int const j)
    {
        constexpr double third = 1. / 3.;
        return a.S0[i] * average_(b, j) + a.dS[i] * (0.5 * b.S0[j] + third * b.dS[j]);
    }


    Segment segment_(double const pos0, double const pos1)
    {
        auto const s0 = shape_(pos0);
        auto const s1 = shape_(pos1);

        Segment seg;
        seg.lo       = std::min(s0.start, s1.start);
        seg.nbrNodes = std::max(s0.start, s1.start) - seg.lo + static_cast<int>(support);

        auto const o0 = s0.start - seg.lo;
        auto const o1 = s1.start - seg.lo;
        for (std::size_t ik = 0; ik < support; ++ik)
        {
            seg.S0[o0 + ik] = s0.weights[ik];
            seg.dS[o0 + ik] -= s0.weights[ik];
            seg.dS[o1 + ik] += s1.weights[ik];
        }
        return seg;
    }


    Shape shape_(double const pos)
    {
        auto const cell  = std::floor(pos);
        auto const delta = pos - cell;

        Shape s;
        s.start = static_cast<int>(cell)
                  - Interpolator<dim, interpOrder>::template computeStartLeftShift<
                      QtyCentering, QtyCentering::primal>(delta);
        weighter_.computeWeight(pos, s.start, s.weights);
        return s;
    }


    static void checkDisplacement_(double const v, double const dt, double const dl,
                                   std::size_t const direction)
    {
        if (std::abs(v) * dt >= dl)
            throw DictionaryException{}("ID", "Esirkepov::displacement")("direction",
                                                                         std::to_string(direction))(
                "v", std::to_string(v))("dt", std::to_string(dt))("dx", std::to_string(dl));
    }


    Weighter<interpOrder> weighter_;
};

} // namespace PHARE::core

#endif
