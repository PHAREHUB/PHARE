#ifndef PHARE_TEST_CORE_DATA_PARTICLES_TEST_PARTICLES_FIXTURES_HPP
#define PHARE_TEST_CORE_DATA_PARTICLES_TEST_PARTICLES_FIXTURES_HPP

#include "core/utilities/types.hpp"
#include "core/utilities/equality.hpp"
#include "core/data/particles/particle.hpp"
#include "core/data/ions/ion_population/particle_pack.hpp"

#include <random>
#include <string>

namespace PHARE::core
{

template<typename ParticleArray_t>
struct UsableParticlesPopulation
{
    template<typename... Args>
    UsableParticlesPopulation(std::string const& _name, Args&&... args)
        : name{_name}
        , domain_particles{args...}
    {
    }

    auto& pack() { return particles_pack; }
    auto& pack() const { return particles_pack; }

    std::string name;
    ParticleArray_t domain_particles;
    ParticleArray_t patch_ghost_particles = domain_particles;
    ParticleArray_t level_ghost_particles = domain_particles;
    ParticlesPack<ParticleArray_t> particles_pack{name,
                                                  &domain_particles, //
                                                  &patch_ghost_particles,
                                                  &level_ghost_particles,
                                                  /*levelGhostParticlesOld=*/nullptr,
                                                  /*levelGhostParticlesNew=*/nullptr};
};

template<std::size_t dim, typename Particle_t = Particle<dim>>
Particle_t particle(std::array<int, dim> const& icell)
{
    return {
        /*.weight = */ .001,
        /*.charge = */ 1,
        /*.iCell  = */ icell,
        /*.delta  = */ ConstArray<double, dim>(.5),
        /*.v      = */ {{.002002002002, .003003003003, .004004004004}} //
    };
}

template<typename Particles, typename Box>
void add_particles_in(Particles& particles, Box const& box, std::size_t const ppc)
{
    for (auto const& bix : box)
        for (std::size_t i = 0; i < ppc; ++i)
            particles.emplace_back(particle(*bix));
}


auto inline rando(std::optional<int> seed = std::nullopt)
{
    if (seed.has_value())
        return std::mt19937_64(*seed);
    std::random_device rd;
    std::seed_seq seed_seq{rd(), rd(), rd(), rd(), rd()};
    return std::mt19937_64(seed_seq);
}

template<typename Particles>
void delta_disperse(Particles& particles, std::optional<int> seed = std::nullopt)
{
    auto gen = rando(seed);
    ParticleDeltaDistribution<double> deltaDistrib;
    for (auto& p : particles)
        p.delta = for_N_make_array<Particles::dimension>([&](auto) { return deltaDistrib(gen); });
}

template<typename Particles>
void vary_velocity(Particles& particles, double const min, double const max,
                   std::optional<int> seed = std::nullopt)
{
    auto gen = rando(seed);
    std::uniform_real_distribution<double> dist{min, max};
    for (auto& p : particles)
        p.v = for_N_make_array<3>([&](auto) { return dist(gen); });
}

template<typename Box_t, typename RValue = std::size_t>
struct CellFlattener
{
    template<typename Icell>
    NO_DISCARD RValue operator()(Icell const& icell) const
    {
        if constexpr (Box_t::dimension == 2)
            return icell[1] + icell[0] * shape[1] * shape[0];
        if constexpr (Box_t::dimension == 3)
            return icell[2] + icell[1] * shape[2] + icell[0] * shape[1] * shape[2];
        return icell[0];
    }

    Box_t const box;
    std::array<int, Box_t::dimension> shape = box.shape().toArray();
};

template<typename I0, typename I1> // support const vs non-const iterators
auto it_dist(I0&& i0, I1&& i1)
{
    return std::distance(i0, i1);
}


template<typename GridLayout>
auto& sort_particles(GridLayout const& layout, auto& particles)
{
    auto constexpr cmp_deltas = [](auto const& a, auto const& b) -> bool {
        return as_tuple(a.delta) < as_tuple(b.delta);
    };
    auto const by_deltas = [&](std::uint64_t const& l, std::uint64_t const& r) {
        std::sort(particles.begin() + l, particles.begin() + r, cmp_deltas);
    };

    auto const ghostBox = grow(layout.AMRBox(), GridLayout::nbrParticleGhosts());
    CellFlattener cell_flattener{ghostBox};

    std::sort(
        particles.begin(), particles.end(),
        [cf = cell_flattener](auto const& a, auto const& b) { return cf(a.iCell) < cf(b.iCell); });

    auto const end = particles.end();
    auto beg       = particles.begin();
    auto lst       = particles.begin();

    auto const check = [&]() { return lst != end and lst->iCell == beg->iCell; };

    while (lst != end)
    {
        lst = beg + 1;
        while (check())
            ++lst;
        auto const s = it_dist(particles.begin(), beg);
        auto const e = it_dist(particles.begin(), lst);
        by_deltas(s, e);
        beg = lst;
    }

    return particles;
}

EqualityReport compare_particles(auto const& ps0, auto const& ps1, double const atol)
{
    if (ps0.size() != ps1.size())
        return EqualityReport{false, "different sizes: " + std::to_string(ps0.size()) + " vs "
                                         + std::to_string(ps1.size())};

    for (std::size_t i = 0; i < ps0.size(); ++i)
    {
        std::string const idx = std::to_string(i);
        if (ps0[i].iCell != ps1[i].iCell)
            return EqualityReport{false, "icell mismatch at index: " + idx, i};

        if (!float_equals(ps0[i].v, ps1[i].v, atol))
            return EqualityReport{false, "v mismatch at index: " + idx, i};

        if (!float_equals(ps0[i].delta, ps1[i].delta, atol))
            return EqualityReport{false, "delta mismatch at index: " + idx, i};
    }

    return EqualityReport{true};
}


} // namespace PHARE::core


#endif /* PHARE_TEST_CORE_DATA_PARTICLES_TEST_PARTICLES_FIXTURES_HPP */
