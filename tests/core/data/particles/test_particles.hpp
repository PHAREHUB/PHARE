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

    if constexpr (std::is_same_v<Particle_t, CountedParticle<dim>>)
        return std::make_from_tuple<Particle_t>(std::tuple_cat(params, std::make_tuple(id)));
    else
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

template<typename Particles, typename GridLayout, typename Box>
void add_ghost_particles(Particles& particles, GridLayout const& layout, Box const& box,
                         std::size_t const ppc)
{
    auto constexpr nghosts = GridLayout::nbrParticleGhosts();

    // PHARE_LOG_LINE_SS(box);
    if constexpr (is_tiled(Particles::layout_mode))
    {
        auto const amr_ghost_box = grow(layout.AMRBox(), nghosts);
        auto const lcl_box       = as_unsigned(shift(box, amr_ghost_box.lower * -1));
        std::size_t id           = particles.size();

        auto tiles_per_cell = make_particle_nd_span_set_from(particles(), layout);
        check_tile_span_set(tiles_per_cell);

        for (auto const& bix : box)
        {
            auto& arr = (*tiles_per_cell.cells((bix - amr_ghost_box.lower).as_unsigned())[0])();
            for (std::size_t i = 0; i < ppc; ++i)
                arr.push_back(particle(*bix, id++));
        }

        particles.template sync<2, ParticleType::Ghost>();
    }
    else
        add_particles_in<ParticleType::Ghost>(particles, box, ppc);
}

template<typename Particles, typename Box>
void add_ghost_particles(Particles& particles, Box const& box, std::size_t const ppc,
                         std::size_t const ghosts)
{
    auto constexpr type = ParticleType::Ghost;

    using enum LayoutMode;
    if constexpr (any_in(Particles::layout_mode, AoSPC, SoAPC))
        particles.template reserve_ppc<ParticleType::Ghost>(ppc);

    std::size_t id = particles.size();
    for (auto const& bix : grow(box, ghosts))
        if (not isIn(bix, box)) // order guaranteed this way
            for (std::size_t i = 0; i < ppc; ++i)
                particles.push_back(particle(*bix, id++));

    if constexpr (any_in(Particles::layout_mode, AoSPC, SoAPC, AoSTS, AoSPCTS))
        particles.template sync<2, type>();
}


template<auto type = ParticleType::Domain, typename Particles, typename Box>
void add_particles_in(Particles& particles, Box const& box, std::size_t const ppc)
{
    static_assert(std::is_same_v<decltype(type), ParticleType>);

    // using enum LayoutMode;
    // if constexpr (any_in(Particles::layout_mode, AoSPC, SoAPC, AoSTS))
    //     particles.template reserve_ppc<type>(ppc);

    std::size_t id = particles.size();
    for (auto const& bix : box)
        for (std::size_t i = 0; i < ppc; ++i)
            particles.emplace_back(particle(*bix, id++));

    if constexpr (any_in(Particles::layout_mode, AoSPC, SoAPC, AoSTS, AoSPCTS))
        particles.template sync<2, type>();
}


template<auto type, typename Src, typename Dst>
void add_particles_from(Src const& src, Dst& dst)
{
    static_assert(std::is_same_v<decltype(type), ParticleType>);
    // SLOW!
    for (auto const& p : src)
        dst.emplace_back(p);

    using enum LayoutMode;
    if constexpr (any_in(Dst::layout_mode, AoSPC, SoAPC, AoSTS))
        dst.template sync<2, type>();
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
