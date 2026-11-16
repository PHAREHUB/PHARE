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
 * In 1D, the longitudinal current Jx is dual and the transverse ones Jy and Jz are primal, as
 * on the Yee layout. Dual node i sits between primal nodes i and i+1.
 *
 * Only 1D is implemented for now.
 */
template<std::size_t dim, std::size_t interpOrder>
class Esirkepov
{
    static_assert(dim == 1, "Esirkepov is only implemented in 1D");

    static constexpr auto support = nbrPointsSupport(interpOrder);
    using Weights_t               = std::array<double, support>;

    struct Shape
    {
        int start;
        Weights_t weights;
    };

public:
    static constexpr auto dimension    = dim;
    static constexpr auto interp_order = interpOrder;


    template<typename Particles, typename VecField, typename GridLayout>
    void operator()(Particles const& particles, VecField& J, GridLayout const& layout,
                    double const charge, double const dt, double const back = 1.,
                    double const front = 1.)
    {
        auto& Jx = J(Component::X);
        auto& Jy = J(Component::Y);
        auto& Jz = J(Component::Z);

        auto const dx       = layout.meshSize()[0];
        auto const halfDt   = 0.5 * dt;
        auto const duration = (back + front) * halfDt;

        for (auto const& particle : particles)
        {
            auto const vx = particle.v[0];
            checkDisplacement_(vx, dt, dx);

            auto const pos = localPosition_(particle, layout);
            auto const s0  = shape_(pos - back * halfDt * vx / dx);
            auto const s1  = shape_(pos + front * halfDt * vx / dx);

            auto const lo = std::min(s0.start, s1.start);
            auto const hi = std::max(s0.start, s1.start) + static_cast<int>(support) - 1;

            assert(lo >= 0 and static_cast<std::size_t>(hi) < Jy.size());

            auto const qw      = charge * particle.weight;
            auto const coefJx  = -qw * dx / duration;
            auto const coefJyz = 0.5 * qw;

            double runningJx = 0.;
            for (int i = lo; i <= hi; ++i)
            {
                auto const w0 = weight_(s0, i);
                auto const w1 = weight_(s1, i);

                Jy(i) += coefJyz * particle.v[1] * (w0 + w1);
                Jz(i) += coefJyz * particle.v[2] * (w0 + w1);

                // the last dual node of the range carries no current, the shapes have equal sums
                if (i < hi)
                {
                    runningJx += coefJx * (w1 - w0);
                    Jx(i) += runningJx;
                }
            }
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
        auto const dx = layout.meshSize()[0];

        for (auto const& particle : particles)
        {
            auto const s
                = shape_(localPosition_(particle, layout) + timeShift * particle.v[0] / dx);

            assert(s.start >= 0 and static_cast<std::size_t>(s.start) + support <= rho.size());

            auto const qw = charge * particle.weight;
            for (std::size_t ik = 0; ik < support; ++ik)
                rho(s.start + ik) += qw * s.weights[ik];
        }
    }



private:
    // position in units of the mesh size, local primal node i being at position i
    template<typename Particle, typename GridLayout>
    static double localPosition_(Particle const& particle, GridLayout const& layout)
    {
        auto const iCell = layout.AMRToLocal(Point{particle.iCell});
        return static_cast<double>(iCell[0]) + particle.delta[0];
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


    static double weight_(Shape const& s, int const i)
    {
        auto const ik = i - s.start;
        return (ik >= 0 and ik < static_cast<int>(support)) ? s.weights[ik] : 0.;
    }


    static void checkDisplacement_(double const vx, double const dt, double const dx)
    {
        if (std::abs(vx) * dt >= dx)
            throw DictionaryException{}("ID", "Esirkepov::displacement")("v", std::to_string(vx))(
                "dt", std::to_string(dt))("dx", std::to_string(dx));
    }


    Weighter<interpOrder> weighter_;
};

} // namespace PHARE::core

#endif
