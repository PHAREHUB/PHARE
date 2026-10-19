#ifndef PHARE_CORE_NUMERICS_COLLISIONS_COLLISIONS_HPP
#define PHARE_CORE_NUMERICS_COLLISIONS_COLLISIONS_HPP


#include "core/def.hpp"
#include "core/utilities/box/box.hpp"
#include "core/utilities/meta/meta_utilities.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numbers>
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
 */
namespace PHARE::core
{
enum class ScatteringLaw { nanbu, takizuka_abe, bobylev, count };


struct CollisionInfo
{
    double const nu0;
    std::int64_t const every;
    ScatteringLaw const law;
    std::uint64_t const seed;

    NO_DISCARD bool isActive() const { return nu0 > 0.; }
    NO_DISCARD bool collidesAt(std::int64_t const step) const
    {
        return isActive() and step % every == 0;
    }
};


struct CollisionCounters
{
    std::size_t pairs   = 0;
    std::size_t skipped = 0; // pairs with no relative velocity, e.g. split siblings
    std::size_t large_s = 0; // s > 1 for a pair moving faster than 0.1 of the block thermal speed
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
        auto const uPerp  = std::hypot(u[0], u[1]);
        auto const uNorm  = std::hypot(uPerp, u[2]);
        auto const sinChi = std::sqrt(std::max(0., 1. - cosChi * cosChi));
        auto const cosPhi = std::cos(phi);
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
        double weight;
        std::uint64_t key;
        ParticleArray_t* particles;
        std::size_t index;
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

        Constexprifier{law}([&]<ScatteringLaw lawV>() {
            for (auto const& cell : cells)
                collideCell_<lawV>(ions, cell, cellKey_(stepKey, cell), dt);
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

    // shuffles a group: sorted by a hash of the block key and of every bit of the particle,
    // a random order that does not depend on how the particles are stored
    static void shuffle_(std::vector<Slot>& slots, std::uint64_t const blockKey)
    {
        for (auto& slot : slots)
        {
            auto const& particle = (*slot.particles)[slot.index];
            auto key             = blockKey;
            for (auto const x : particle.delta)
                key = collisions::combine(key, std::bit_cast<std::uint64_t>(x));
            for (auto const x : slot.v)
                key = collisions::combine(key, std::bit_cast<std::uint64_t>(x));
            slot.key = key;
        }
        std::sort(slots.begin(), slots.end(),
                  [](auto const& a, auto const& b) { return a.key < b.key; });
    }

    static double thermalSpeed_(std::vector<Slot> const& a, std::vector<Slot> const& b)
    {
        double w = 0, mean[3] = {0, 0, 0}, sq = 0;
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

    template<ScatteringLaw law>
    void collidePair_(Slot& a, Slot& b, Group const& ga, Group const& gb, double const strength,
                      double const vth, collisions::Stream& draw)
    {
        std::array<double, 3> const u{a.v[0] - b.v[0], a.v[1] - b.v[1], a.v[2] - b.v[2]};
        auto const uNorm = std::sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
        auto const vNorm = std::hypot(a.v[0], a.v[1], a.v[2]) + std::hypot(b.v[0], b.v[1], b.v[2]);

        ++counters_.pairs;
        if (uNorm <= 1e-14 * vNorm)
        {
            ++counters_.skipped;
            return;
        }

        auto const wMax = std::max(a.weight, b.weight);
        auto const mab  = ga.mass * gb.mass / (ga.mass + gb.mass);
        auto const qq   = ga.charge * gb.charge;
        auto const s    = nu0 * qq * qq * wMax * strength / (mab * mab * uNorm * uNorm * uNorm);
        if (s > 1. and uNorm > 0.1 * vth)
            ++counters_.large_s;

        auto const cosChi  = collisions::cosChi<law>(s, draw);
        auto const phi     = 2. * std::numbers::pi * draw();
        auto const updateA = b.weight >= wMax or draw() * wMax < b.weight;
        auto const updateB = a.weight >= wMax or draw() * wMax < a.weight;
        collisions::scatter(a.v, b.v, ga.mass, gb.mass, cosChi, phi, updateA, updateB);
    }

    template<ScatteringLaw law>
    void collideCell_(Ions& ions, auto const& cell, std::uint64_t const key, double const dt)
    {
        for (auto& group : groups_)
            group.slots.clear();

        std::size_t ipop = 0;
        for (auto& pop : ions)
        {
            auto& particles = pop.domainParticles();
            auto& slots     = groups_[groupOf_[ipop++]].slots;
            for (auto const index : particles.cellIndexes(cell))
                slots.push_back(
                    {particles[index].v, particles[index].weight, 0, &particles, index});
        }

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
            auto const vth      = thermalSpeed_(slots, {});
            shuffle_(slots, blockKey);
            collisions::forEachIntraPair(n, [&](auto const i, auto const j, double const factor) {
                collidePair_<law>(slots[i], slots[j], g, g, factor * (n - 1) * dt, vth, draw);
            });
        }

        for (std::size_t ig = 0; ig < groups_.size(); ++ig)
            for (std::size_t ih = ig + 1; ih < groups_.size(); ++ih)
            {
                ++block;
                auto& g = groups_[ig];
                auto& h = groups_[ih];
                if (g.slots.empty() or h.slots.empty())
                    continue;

                auto const blockKey = collisions::combine(key, block);
                auto draw           = collisions::Stream{blockKey};
                auto const vth      = thermalSpeed_(g.slots, h.slots);
                shuffle_(g.slots, blockKey);
                shuffle_(h.slots, collisions::combine(blockKey, block));

                auto const gIsLarge = g.slots.size() >= h.slots.size();
                auto& large         = gIsLarge ? g : h;
                auto& small         = gIsLarge ? h : g;
                auto const strength = small.slots.size() * dt;
                collisions::forEachInterPair(
                    large.slots.size(), small.slots.size(), [&](auto const i, auto const j) {
                        collidePair_<law>(large.slots[i], small.slots[j], large, small, strength,
                                          vth, draw);
                    });
            }

        for (auto const& g : groups_)
            for (auto const& slot : g.slots)
                (*slot.particles)[slot.index].v = slot.v;
    }

    std::vector<Group> groups_;
    std::vector<std::size_t> groupOf_; // population index -> group index
    CollisionCounters counters_;
};

} // namespace PHARE::core

#endif
