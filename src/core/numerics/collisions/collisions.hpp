#ifndef PHARE_CORE_NUMERICS_COLLISIONS_COLLISIONS_HPP
#define PHARE_CORE_NUMERICS_COLLISIONS_COLLISIONS_HPP


#include "core/def.hpp"
#include "core/utilities/box/box.hpp"
#include "core/utilities/meta/meta_utilities.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numbers>
#include <tuple>
#include <numeric>
#include <utility>
#include <vector>


/* Binary Monte-Carlo Coulomb collisions between the macro-ions of one cell.
 *
 * In each cell, particles are pooled in species groups (populations with the same mass and
 * charge), then paired within each group and between each pair of groups (Takizuka & Abe 1977).
 * The deflection of a pair is drawn from the scattering parameter
 *
 *     s = nu0 Za^2 Zb^2 n_eff dt / (m_ab^2 u^3)
 *
 * in PHARE units, nu0 being the Coulomb logarithm times e^4 n0 / (4 pi eps0^2 mp^2 VA^3 Omega_p).
 * Unequal weights follow Higginson et al. 2020 as written by Angus et al. 2025: n_eff uses the
 * larger weight of the pair and each partner is deflected with probability w_other / w_max.
 *
 * Random draws are pure functions of (seed, level, step, cell, block, draw index), so results do
 * not depend on how the level is split into patches and ranks.
 *
 * Accuracy: relaxation rates converge to the Landau-Fokker-Planck ones as dt -> 0, the error
 * growing like sqrt(s_th) (Wang et al. 2008), s_th being s at the thermal relative speed:
 * roughly 0.23 sqrt(s_th) for Nanbu's law, 1.5 times that for Takizuka & Abe's, half that for
 * Bobylev's, twice that with weights spread over two decades (measured on temperature
 * isotropization within about 30 %, see the tests).
 */
namespace PHARE::core
{
enum class ScatteringLaw { nanbu, takizuka_abe, bobylev, count };


struct CollisionInfo
{
    double const nu0;
    std::int64_t const every;
    ScatteringLaw const law;
    bool const exactConservation;
    std::uint64_t const seed;

    NO_DISCARD bool isActive() const { return nu0 > 0.; }
    NO_DISCARD bool collidesAt(std::int64_t const step) const
    {
        assert(every > 0);
        return isActive() and step % every == 0;
    }
};


struct CollisionCounters
{
    std::size_t pairs   = 0;
    std::size_t skipped = 0; // pairs with no relative velocity, e.g. split siblings
    std::size_t large_s = 0; // s > 1 for a pair moving faster than 0.1 of the block thermal speed
    std::size_t unconverged = 0; // blocks whose energy could not be fully restored
    double residualEnergy   = 0; // energy those blocks gained, in units of n0 mp VA^2
};


namespace collisions
{
    // stateless 64-bit mixer, the output function of SplitMix64 (Steele et al. 2014)
    NO_DISCARD constexpr std::uint64_t splitmix64(std::uint64_t x)
    {
        x += 0x9e3779b97f4a7c15ULL;
        x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
        x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
        return x ^ (x >> 31);
    }

    NO_DISCARD constexpr std::uint64_t combine(std::uint64_t const key, std::uint64_t const value)
    {
        return splitmix64(key ^ splitmix64(value));
    }

    // uniform in (0, 1], never 0 since the scattering laws take its logarithm
    NO_DISCARD constexpr double toUnit(std::uint64_t const x)
    {
        return static_cast<double>((x >> 11) + 1) * 0x1.0p-53;
    }

    // draw number i of the stream `key` (a Weyl sequence through the mixer)
    struct Stream
    {
        std::uint64_t const key;
        std::uint64_t index = 0;

        double operator()() { return toUnit(splitmix64(key + index++ * 0x9e3779b97f4a7c15ULL)); }
    };


    template<ScatteringLaw law>
    NO_DISCARD double cosChi(double const s, Stream& draw);

    // Nanbu 1997 cumulative angle, A(s) fitted by Perez et al. 2012 (as in WarpX and PICNIC)
    template<>
    inline double cosChi<ScatteringLaw::nanbu>(double const s, Stream& draw)
    {
        auto const fromA = [&](double const A) {
            return std::log(std::exp(-A) + 2. * draw() * std::sinh(A)) / A;
        };

        if (s < 0.1)
            while (true) // the law is truncated at -1
                if (auto const c = 1. + s * std::log(draw()); c >= -1.)
                    return c;
        if (s < 3.)
            return fromA(1.
                         / (0.0056958 + 0.9560202 * s - 0.508139 * s * s + 0.47913906 * s * s * s
                            - 0.12788975 * s * s * s * s + 0.02389567 * s * s * s * s * s));
        if (s < 6.)
            return fromA(3. * std::exp(-s));
        return 2. * draw() - 1.;
    }

    // Takizuka & Abe 1977: tan(chi/2) is gaussian with variance s/2
    template<>
    inline double cosChi<ScatteringLaw::takizuka_abe>(double const s, Stream& draw)
    {
        if (s >= 6.)
            return 2. * draw() - 1.;
        auto const r     = std::sqrt(-2. * std::log(draw()));
        auto const delta = std::sqrt(0.5 * s) * r * std::cos(2. * std::numbers::pi * draw());
        return (1. - delta * delta) / (1. + delta * delta);
    }

    // Bobylev & Nanbu 2000: one fixed angle with the mean deflection of the cumulative law
    template<>
    inline double cosChi<ScatteringLaw::bobylev>(double const s, Stream& /*draw*/)
    {
        return 1. - std::min(s, 2.);
    }


    // rotates the relative velocity u = va - vb by the polar angle chi and the azimuth phi
    // (Takizuka & Abe 1977) and shares the change between the partners in proportion to
    // the other's mass, so that the momentum and energy of the pair are kept
    inline void scatter(std::array<double, 3>& va, std::array<double, 3>& vb, double const ma,
                        double const mb, double const cosChi, double const phi, bool const updateA,
                        bool const updateB)
    {
        std::array<double, 3> const u{va[0] - vb[0], va[1] - vb[1], va[2] - vb[2]};
        auto const uPerp  = std::sqrt(u[0] * u[0] + u[1] * u[1]);
        auto const uNorm  = std::sqrt(uPerp * uPerp + u[2] * u[2]);
        auto const sinChi = std::sqrt(std::max(0., 1. - cosChi * cosChi));
        auto const cosPhi = std::cos(phi); // the compiler fuses these into one sincos
        auto const sinPhi = std::sin(phi);

        std::array<double, 3> du;
        if (uPerp > 1e-12 * uNorm)
        {
            du[0] = (u[0] * u[2] * sinChi * cosPhi - u[1] * uNorm * sinChi * sinPhi) / uPerp
                    - u[0] * (1. - cosChi);
            du[1] = (u[1] * u[2] * sinChi * cosPhi + u[0] * uNorm * sinChi * sinPhi) / uPerp
                    - u[1] * (1. - cosChi);
            du[2] = -uPerp * sinChi * cosPhi - u[2] * (1. - cosChi);
        }
        else // u along z, any perpendicular axes do
        {
            du[0] = uNorm * sinChi * cosPhi;
            du[1] = uNorm * sinChi * sinPhi;
            du[2] = -u[2] * (1. - cosChi);
        }

        auto const mTotal = ma + mb;
        for (std::size_t i = 0; i < 3; ++i)
        {
            if (updateA)
                va[i] += mb / mTotal * du[i];
            if (updateB)
                vb[i] -= ma / mTotal * du[i];
        }
    }


    // pairs of indexes in 0..n-1, consecutive in a shuffled group. An odd group starts with the
    // triangle (0,1), (1,2), (2,0) at half strength (Takizuka & Abe 1977).
    void forEachIntraPair(std::size_t const n, auto&& fn)
    {
        auto first = std::size_t{0};
        if (n % 2 == 1 and n >= 3)
        {
            fn(0, 1, 0.5);
            fn(1, 2, 0.5);
            fn(2, 0, 0.5);
            first = 3;
        }
        for (auto i = first; i + 1 < n; i += 2)
            fn(i, i + 1, 1.);
    }

    // every particle of the larger group meets one of the smaller group, used in turn
    void forEachInterPair(std::size_t const nLarge, std::size_t const nSmall, auto&& fn)
    {
        for (std::size_t i = 0; i < nLarge; ++i)
            fn(i, i % nSmall);
    }

    // Fisher-Yates, drawing from `draw`
    void shuffle(auto& range, Stream& draw)
    {
        for (auto i = range.size(); i > 1; --i)
            std::swap(range[i - 1], range[std::min(static_cast<std::size_t>(draw() * i), i - 1)]);
    }

    // n_eff / w_max for a pair, s being nu0 (Za Zb)^2 n_eff dt / (m_ab^2 u^3): within a group
    // of n, the n - 1 others (halved in the odd triangle); between groups, the size of the
    // smaller one (Angus et al. 2025, Eq. 12)
    NO_DISCARD inline double intraStrength(std::size_t const n, double const factor)
    {
        return factor * static_cast<double>(n - 1);
    }

    NO_DISCARD inline double interStrength(std::size_t const nSmall)
    {
        return static_cast<double>(nSmall);
    }

    // probabilities that each partner of a pair is deflected, s being built on the larger
    // weight: the lighter partner always, the heavier one with the ratio of the weights
    struct Acceptance
    {
        double a, b;
    };

    NO_DISCARD inline Acceptance acceptance(double const wa, double const wb)
    {
        auto const wMax = std::max(wa, wb);
        return {wb / wMax, wa / wMax};
    }
} // namespace collisions



template<typename Ions>
class Collisions : public CollisionInfo
{
    using ParticleArray_t           = Ions::particle_array_type;
    static constexpr auto dimension = ParticleArray_t::dimension;
    using Box_t                     = Box<int, dimension>;

    struct Slot
    {
        std::array<double, 3> v;
        std::array<double, dimension> delta;
        double weight;
        std::uint64_t key;
        ParticleArray_t* particles;
        std::size_t index;
        std::uint64_t population;
    };

    struct Group
    {
        double mass;
        double charge;
        std::vector<Slot> slots; // the group's particles in the current cell, capacity kept
    };

public:
    explicit Collisions(CollisionInfo const& info)
        : CollisionInfo{info}
    {
    }

    // collides the domain particles of all populations in every cell of `cells`; one Collisions
    // object per thread, its buffers are reused from cell to cell
    void operator()(Ions& ions, Box_t const& cells, std::uint64_t const stepKey, double const dt)
    {
        if (groups_.empty())
            makeGroups_(ions);

        Constexprifier{law, exactConservation}([&]<ScatteringLaw lawV, bool exactV>() {
            for (auto const& cell : cells)
                collideCell_<lawV, exactV>(ions, cell, cellKey_(stepKey, cell), dt);
        });
    }

    NO_DISCARD auto const& counters() const { return counters_; }


private:
    void makeGroups_(Ions const& ions)
    {
        for (auto const& pop : ions)
        {
            auto const mass   = pop.mass();
            auto const charge = pop.particleInitializerInfo()["charge"].template to<double>();
            auto const same   = std::find_if(groups_.begin(), groups_.end(), [&](auto const& g) {
                return g.mass == mass and g.charge == charge;
            });
            groupOf_.push_back(static_cast<std::size_t>(same - groups_.begin()));
            if (same == groups_.end())
                groups_.push_back({mass, charge, {}});
        }
    }

    static std::uint64_t cellKey_(std::uint64_t key, auto const& cell)
    {
        for (auto const c : cell)
            key = collisions::combine(key, static_cast<std::uint64_t>(c));
        return key;
    }

    // sorts a group by a hash of the key, of the population and of every bit of the particle:
    // a random order that does not depend on how the particles are stored. Only bitwise
    // identical particles of one population tie, and their order does not matter.
    void sortRandomly_(std::vector<Slot>& slots, std::uint64_t const blockKey)
    {
        // words folded by an odd multiplication, a bijection, then one full mix
        auto const hash = [&](Slot const& slot) {
            auto h         = blockKey ^ slot.population;
            auto const add = [&](double const x) {
                h = (h ^ std::bit_cast<std::uint64_t>(x)) * 0x9e3779b97f4a7c15ULL;
            };
            add(slot.weight);
            for (auto const x : slot.delta)
                add(x);
            for (auto const x : slot.v)
                add(x);
            return collisions::splitmix64(h);
        };

        for (auto& slot : slots)
            slot.key = hash(slot);
        sortBy_(slots, [](auto const& a, auto const& b) { return a.key < b.key; });
    }

    // sorts through an index, cheaper than moving slots around
    void sortBy_(std::vector<Slot>& slots, auto&& less)
    {
        index_.resize(slots.size());
        std::iota(index_.begin(), index_.end(), std::size_t{0});
        std::sort(index_.begin(), index_.end(),
                  [&](auto const i, auto const j) { return less(slots[i], slots[j]); });

        sorted_.clear();
        for (auto const i : index_)
            sorted_.push_back(slots[i]);
        slots.swap(sorted_);
    }

    static double thermalSpeed_(std::vector<Slot> const& a) { return thermalSpeed_(a, {}); }

    static double thermalSpeed_(std::vector<Slot> const& a, std::vector<Slot> const& b)
    {
        double w = 0, sq = 0;
        std::array<double, 3> mean{};
        for (auto const* slots : {&a, &b})
            for (auto const& slot : *slots)
            {
                w += slot.weight;
                for (std::size_t i = 0; i < 3; ++i)
                {
                    mean[i] += slot.weight * slot.v[i];
                    sq += slot.weight * slot.v[i] * slot.v[i];
                }
            }
        auto const meanSq = (mean[0] * mean[0] + mean[1] * mean[1] + mean[2] * mean[2]) / (w * w);
        return std::sqrt(std::max(0., sq / w - meanSq) / 3.);
    }

    // weighted sums of a block, m_p = m w_p being the mass of a macro-particle
    struct Sums
    {
        std::array<double, 3> momentum{}; // sum m_p v
        double energy = 0;                // sum m_p v^2 / 2
        double mass   = 0;                // sum m_p
        double massW  = 0;                // sum m_p w
    };

    static Sums sums_(auto const&... groups)
    {
        Sums sums;
        auto const add = [&](Group const& g) {
            for (auto const& slot : g.slots)
            {
                auto const mp = g.mass * slot.weight;
                for (std::size_t i = 0; i < 3; ++i)
                {
                    sums.momentum[i] += mp * slot.v[i];
                    sums.energy += 0.5 * mp * slot.v[i] * slot.v[i];
                }
                sums.mass += mp;
                sums.massW += mp * slot.weight;
            }
        };
        (add(groups), ...);
        return sums;
    }

    // Angus et al. 2025, Eqs. 17-24: the weighted acceptance conserves momentum and energy only
    // on average. A velocity shift proportional to the weight restores the momentum the block had
    // before its pairs scattered, then zero-angle exchanges of relative energy between neighbours
    // of each group restore its energy, shared between two groups in proportion to their mean
    // weight times their energy in the block's centre of mass frame.
    // equal weights: every pair conserves exactly, nothing to restore
    static bool equalWeights_(auto const&... groups)
    {
        auto const w   = std::get<0>(std::forward_as_tuple(groups...)).slots.front().weight;
        auto const all = [&](Group const& g) {
            return std::all_of(g.slots.begin(), g.slots.end(),
                               [&](auto const& slot) { return slot.weight == w; });
        };
        return (all(groups) and ...);
    }

    void conserve_(Sums const& before, Group& g)
    {
        if (equalWeights_(g))
            return;
        shiftMomentum_(before, g);
        record_(absorb_(g, sums_(g).energy - before.energy));
    }

    void conserve_(Sums const& before, Group& g, Group& h)
    {
        if (equalWeights_(g, h))
            return;
        shiftMomentum_(before, g, h);
        auto const dE = sums_(g, h).energy - before.energy;

        auto const centreOfMass = [&](std::size_t i) { return before.momentum[i] / before.mass; };
        auto const share        = [&](Group const& group) {
            double energy = 0, weight = 0;
            for (auto const& slot : group.slots)
            {
                weight += slot.weight;
                for (std::size_t i = 0; i < 3; ++i)
                    energy += 0.5 * group.mass * slot.weight
                              * std::pow(slot.v[i] - centreOfMass(i), 2);
            }
            return weight / group.slots.size() * energy;
        };

        auto const shareG = share(g), shareH = share(h);
        auto const dEg = g.slots.size() == 1   ? 0.
                         : h.slots.size() == 1 ? dE
                                               : dE * shareG / (shareG + shareH);
        record_(absorb_(g, dEg) + absorb_(h, dE - dEg));
    }

    static void shiftMomentum_(Sums const& before, auto&... groups)
    {
        auto const after = sums_(groups...);
        auto const shift = [&](Group& g) {
            for (auto& slot : g.slots)
                for (std::size_t i = 0; i < 3; ++i)
                    slot.v[i]
                        -= slot.weight * (after.momentum[i] - before.momentum[i]) / after.massW;
        };
        (shift(groups), ...);
    }

    // removes dE from the group through zero-angle exchanges between neighbours, heaviest first,
    // and returns what could not be removed
    double absorb_(Group& g, double dE)
    {
        auto constexpr maxSweeps = 10;
        auto& slots              = g.slots;
        auto const n             = slots.size();
        if (dE == 0. or n < 2)
            return dE;

        // heaviest first, ties in the shuffled order: a total order, the same on every platform
        sortBy_(slots, [](auto const& a, auto const& b) {
            return std::tuple{-a.weight, a.key} < std::tuple{-b.weight, b.key};
        });
        for (auto sweep = 0; sweep < maxSweeps and dE != 0.; ++sweep)
            for (auto p = static_cast<std::size_t>(sweep % 2); p < n and dE != 0.; p += 2)
                exchange_(slots[p], slots[(p + 1) % n], g.mass, dE);
        return dE;
    }

    // inelastic event at zero angle: the pair gives up U of its relative energy, at most 5 %,
    // and keeps its momentum
    static void exchange_(Slot& a, Slot& b, double const mass, double& dE)
    {
        auto constexpr energyFraction = 0.05;
        auto const ma                 = mass * a.weight;
        auto const mb                 = mass * b.weight;
        auto const mu                 = ma * mb / (ma + mb);
        std::array const u{a.v[0] - b.v[0], a.v[1] - b.v[1], a.v[2] - b.v[2]};
        auto const eRel = 0.5 * mu * (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
        if (eRel <= 0.)
            return;

        auto U = std::copysign(energyFraction * eRel, dE);
        if (std::abs(U) >= std::abs(dE))
            U = std::exchange(dE, 0.);
        else
            dE -= U;

        auto const factor = std::sqrt(1. - U / eRel) - 1.;
        for (std::size_t i = 0; i < 3; ++i)
        {
            a.v[i] += mu / ma * factor * u[i];
            b.v[i] -= mu / mb * factor * u[i];
        }
    }

    void record_(double const residual)
    {
        if (residual != 0.)
        {
            ++counters_.unconverged;
            counters_.residualEnergy += residual;
        }
    }

    template<ScatteringLaw law>
    void collidePair_(Slot& a, Slot& b, Group const& ga, Group const& gb, double const strength,
                      double const vth, collisions::Stream& draw)
    {
        auto const squared = [](auto const& x) { return x[0] * x[0] + x[1] * x[1] + x[2] * x[2]; };
        auto const uSq     = squared(std::array{a.v[0] - b.v[0], a.v[1] - b.v[1], a.v[2] - b.v[2]});

        ++counters_.pairs;
        if (uSq <= 1e-28 * (squared(a.v) + squared(b.v))) // split siblings
        {
            ++counters_.skipped;
            return;
        }
        auto const uNorm = std::sqrt(uSq);

        auto const wMax = std::max(a.weight, b.weight);
        auto const mab  = ga.mass * gb.mass / (ga.mass + gb.mass);
        auto const qq   = ga.charge * gb.charge;
        auto const s    = nu0 * qq * qq * wMax * strength / (mab * mab * uNorm * uNorm * uNorm);
        if (s > 1. and uNorm > 0.1 * vth)
            ++counters_.large_s;

        auto const cosChi   = collisions::cosChi<law>(s, draw);
        auto const phi      = 2. * std::numbers::pi * draw();
        auto const [pA, pB] = collisions::acceptance(a.weight, b.weight);
        auto const updateA  = pA >= 1. or draw() < pA;
        auto const updateB  = pB >= 1. or draw() < pB;
        collisions::scatter(a.v, b.v, ga.mass, gb.mass, cosChi, phi, updateA, updateB);
    }

    template<ScatteringLaw law, bool exact>
    void collideCell_(Ions& ions, auto const& cell, std::uint64_t const key, double const dt)
    {
        for (auto& group : groups_)
            group.slots.clear();

        std::uint64_t ipop = 0;
        for (auto& pop : ions)
        {
            auto& particles = pop.domainParticles();
            auto& slots     = groups_[groupOf_[ipop]].slots;
            for (auto const index : particles.cellIndexes(cell))
            {
                auto const& particle = particles[index];
                slots.push_back(
                    {particle.v, particle.delta, particle.weight, 0, &particles, index, ipop});
            }
            ++ipop;
        }

        // one sort per group, then each block between groups shuffles that order further
        std::uint64_t block = 0;
        for (auto& g : groups_)
        {
            auto& slots = g.slots;
            ++block;
            if (slots.size() < 2)
                continue;

            auto const blockKey = collisions::combine(key, block);
            auto draw           = collisions::Stream{blockKey};
            auto const n        = slots.size();
            auto const vth      = thermalSpeed_(slots);
            auto const before   = exact ? sums_(g) : Sums{};
            sortRandomly_(slots, blockKey);
            collisions::forEachIntraPair(n, [&](auto const i, auto const j, double const factor) {
                collidePair_<law>(slots[i], slots[j], g, g,
                                  collisions::intraStrength(n, factor) * dt, vth, draw);
            });
            if constexpr (exact)
                conserve_(before, g);
        }

        for (std::size_t ig = 0; ig < groups_.size(); ++ig)
            for (std::size_t ih = ig + 1; ih < groups_.size(); ++ih)
            {
                ++block;
                auto& g = groups_[ig];
                auto& h = groups_[ih];
                if (g.slots.empty() or h.slots.empty())
                    continue;

                auto draw         = collisions::Stream{collisions::combine(key, block)};
                auto const vth    = thermalSpeed_(g.slots, h.slots);
                auto const before = exact ? sums_(g, h) : Sums{};
                collisions::shuffle(g.slots, draw);
                collisions::shuffle(h.slots, draw);

                auto const gIsLarge = g.slots.size() >= h.slots.size();
                auto& large         = gIsLarge ? g : h;
                auto& small         = gIsLarge ? h : g;
                auto const strength = collisions::interStrength(small.slots.size()) * dt;
                collisions::forEachInterPair(
                    large.slots.size(), small.slots.size(), [&](auto const i, auto const j) {
                        collidePair_<law>(large.slots[i], small.slots[j], large, small, strength,
                                          vth, draw);
                    });
                if constexpr (exact)
                    conserve_(before, g, h);
            }

        for (auto const& g : groups_)
            for (auto const& slot : g.slots)
                (*slot.particles)[slot.index].v = slot.v;
    }

    std::vector<Group> groups_;
    std::vector<std::size_t> groupOf_; // population index -> group index
    std::vector<std::size_t> index_;   // sort scratch
    std::vector<Slot> sorted_;         // sort scratch
    CollisionCounters counters_;
};

} // namespace PHARE::core

#endif
