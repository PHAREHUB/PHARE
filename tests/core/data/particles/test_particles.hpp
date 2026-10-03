#ifndef PHARE_CORE_DATA_TEST_PARTICLES_HPP
#define PHARE_CORE_DATA_TEST_PARTICLES_HPP

#include "core/def.hpp"
#include "phare_core.hpp"
#include "core/utilities/types.hpp"
#include "core/data/particles/particle_array_def.hpp"

#include "core/data/particles/particle_array_sorter.hpp"
#include "core/data/particles/particle_array_appender.hpp"
#include "core/data/particles/particle_array_converter.hpp"
#include "core/data/particles/particle_array_comparator.hpp"



#include <fstream>

namespace PHARE::core
{
template<std::size_t dim, typename Particle_t = Particle<dim>>
Particle_t particle(std::array<int, dim> const& icell, [[maybe_unused]] std::size_t const id = 0)
{
    using Tup = std::tuple<double, double, std::array<int, dim>, std::array<double, dim>,
                           std::array<double, 3>>;
    Tup params{
        /*.weight = */ .001,
        /*.charge = */ 1,
        /*.iCell  = */ icell,
        /*.delta  = */ ConstArray<double, dim>(.51),
        /*.v      = */ {{.002002002002, .003003003003, .004004004004}} //
    };

    return std::make_from_tuple<Particle_t>(params);
}

template<std::size_t dim>
Particle<dim> particle(int const icell = 15)
{
    return particle(ConstArray<int, dim>(icell));
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
void shuffle(Particles& particles, std::optional<int> seed = std::nullopt)
{
    auto gen = rando(seed);

    std::shuffle(particles.begin(), particles.end(), gen);

    if constexpr (any_in(Particles::layout_mode, LayoutMode::AoSMapped))
    {
        particles.empty_map();
        particles.map_particles();
    }
}

template<typename Particles>
void delta_disperse(Particles& particles, std::optional<int> seed = std::nullopt)
{
    auto gen = rando(seed);
    ParticleDeltaDistribution<double> deltaDistrib;
    for (auto& p : particles)
        p.delta() = core::ConstArrayFrom<Particles::dimension>([&] { return deltaDistrib(gen); });
}

template<typename Particles>
void vary_velocity(Particles& particles, double const min, double const max,
                   std::optional<int> seed = std::nullopt)
{
    auto gen = rando(seed);
    std::uniform_real_distribution<double> dist{min, max};
    for (auto& p : particles)
        p.v() = core::ConstArrayFrom<3>([&] { return dist(gen); });
}


template<typename Particles, typename Box>
void add_ghost_particles(Particles& particles, Box const& box, std::size_t const ppc)
{
    add_particles<ParticleType::Ghost>(particles, box, ppc);
}

template<typename Particles, typename Box>
void add_ghost_particles(Particles& particles, Box const& box, std::size_t const ppc,
                         std::size_t const ghosts)
{
    // add_particles deals in one plain box at a time, so the ghost ring (which excludes
    // the interior) is decomposed into disjoint fragments first, one call per fragment.
    for (auto const& fragment : grow(box, ghosts).remove(box))
        add_particles<ParticleType::Ghost>(particles, fragment, ppc);
}


template<auto type = ParticleType::Domain, typename Particles, typename Box>
void add_particles(Particles& particles, Box const& box, std::size_t const ppc = 100)
{
    static_assert(std::is_same_v<decltype(type), ParticleType>);
    core::add_particles<type>(particles, box, ppc,
                              [](auto const& icell, auto const id) { return particle(icell, id); });
}


template<auto type, typename Src, typename Dst>
void add_particles_from(Src const& src, Dst& dst)
{
    static_assert(std::is_same_v<decltype(type), ParticleType>);
    // SLOW!
    for (auto const& p : src)
        dst.emplace_back(p);

    dst.template on_appended<type>();
}


template<typename Particles>
void write_raw_to_file(Particles const& particles, std::string const& filename)
{
    using Particle_t = typename Particles::value_type;
    std::ofstream f{filename, std::ios::binary};
    f.write(reinterpret_cast<char const*>(particles.vector().data()),
            particles.size() * sizeof(Particle_t));
}

template<typename Particles>
Particles& read_raw_from_file(Particles& particles, std::string const& filename)
{
    using Particle_t = typename Particles::value_type;

    std::ifstream f{filename, std::ios::binary};

    // Stop eating new lines in binary mode
    f.unsetf(std::ios::skipws);

    // get its size:
    std::streampos fileSize;
    f.seekg(0, std::ios::end);
    fileSize = f.tellg();
    f.seekg(0, std::ios::beg);
    particles.resize(fileSize / sizeof(Particle_t));

    // read the data:
    f.read(reinterpret_cast<char*>(particles.vector().data()),
           particles.size() * sizeof(Particle_t));

    return particles;
}


template<typename Particles>
auto read_raw_from_file(std::string const& filename)
{
    Particles particles;
    return read_raw_from_file(particles, filename);
}




template<typename ParticleArray>
std::size_t memory_for_particles(ParticleArray const& ps)
{
    return ParticleArray::size_of_particle() * ps.size();
}

template<typename ParticleArray>
std::size_t ram_in_mbs(ParticleArray const& ps)
{
    return memory_for_particles(ps) / 1e6;
}

template<typename ParticleArray>
std::size_t ram_in_gbs(ParticleArray const& ps)
{
    return memory_for_particles(ps) / 1e9;
}

template<typename P0>
struct ParticleComparable
{
    template<typename P1>
    bool operator<(P1 const& p1)
    {
        auto const ic0 = p0.iCell();
        auto const ic1 = p1.iCell();

        if (ic0 > ic1)
            return false;

        if (ic0 == ic1)
            return as_tuple(p0.delta()) < as_tuple(p1.delta());

        return true;
    }

    P0 const& p0;
};

template<typename ParticleArray0, typename ParticleArray1>
std::size_t count_equal(ParticleArray0 const& p0, ParticleArray1 const& p1)
{
    PHARE_LOG_LINE_SS("");
    auto inc = [](auto&... is) { (++is, ...); };

    std::size_t i0 = 0, i1 = 0, eq = 0;

    while (i0 < p0.size() and i1 < p1.size())
    {
        PHARE_LOG_LINE_SS(i0 << " " << i1 << " " << eq);
        PHARE_LOG_LINE_SS(Point{p0[i0].iCell()} << " " << Point{p1[i1].iCell()});

        auto const eqr = particle_compare(p0[i0], p1[i1]);
        PHARE_LOG_LINE_SS(eqr.why());

        if (eqr)
        {
            PHARE_LOG_LINE_SS(i0 << " " << i1 << " " << eq);
            inc(i0, i1, eq);
        }
        else
        {
            PHARE_LOG_LINE_SS(i0 << " " << i1 << " " << eq);
            inc(ParticleComparable{p0[i0]} < p1[i1] ? ++i0 : ++i1);
        }
    }

    PHARE_LOG_LINE_SS(eq);

    return eq;
}

template<typename ParticleArray0, typename ParticleArray1, typename Shift>
std::size_t count_equal(ParticleArray0 const& p0, ParticleArray1 const& p1, Shift shift)
{
    auto inc = [](auto&... is) { (++is, ...); };

    std::size_t i0 = 0, i1 = 0, eq = 0;

    while (i0 < p0.size() and i1 < p1.size())
    {
        auto const p = shift(p0[i0]);

        if (particle_compare(p, p1[i1]))
            inc(i0, i1, eq);
        else
            inc(ParticleComparable{p} < p1[i1] ? ++i0 : ++i1);
    }

    return eq;
}

} // namespace PHARE::core

#endif /* PHARE_CORE_DATA_TEST_PARTICLES_HPP */
