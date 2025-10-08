#ifndef PHARE_BENCH_CORE_BENCH_H
#define PHARE_BENCH_CORE_BENCH_H

#include "phare_core.hpp"

#include "core/utilities/box/box.hpp"

#include "tests/core/data/vecfield/test_vecfield_fixtures.hpp"
#include "tests/core/data/electromag/test_electromag_fixtures.hpp"

#include "benchmark/benchmark.h"

#include <fstream>

namespace PHARE::core::bench
{

template<std::size_t dim>
PHARE::core::Particle<dim> particle(int icell = 15)
{
    return {/*.weight = */ 0,
            /*.charge = */ 1,
            /*.iCell  = */ PHARE::core::ConstArray<int, dim>(icell),
            /*.delta  = */ PHARE::core::ConstArray<double, dim>(.5),
            /*.v      = */ {{.00001, .00001, .00001}}};
}

template<typename Particles>
void write_raw_to_file(Particles const& particles, std::string const& filename)
{
    using Particle_t = typename Particles::value_type;
    std::ofstream f{filename};
    f.write(reinterpret_cast<char const*>(particles.vector().data()),
            particles.size() * sizeof(Particle_t));
}

template<typename Particles>
void read_raw_from_file(Particles& particles, std::string const& filename)
{
    using Particle_t = typename Particles::value_type;

    std::ifstream f{filename};

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
}

template<std::size_t dim, typename ParticleArray_t = PHARE::core::ParticleArray<dim>>
auto make_particles(std::size_t n_particles)
{
    auto particles = ParticleArray_t{Box<int, dim>{}};
    particles.vector().resize(n_particles, particle<dim>());
    return particles;
}

template<typename Particles, typename Point>
void disperse(Particles& particles, Point lo, Point up, std::optional<int> seed = std::nullopt)
{
    auto gen = [&]() {
        if (!seed.has_value())
        {
            std::random_device rd;
            std::seed_seq seed_seq{rd(), rd(), rd(), rd(), rd(), rd(), rd()};
            return std::mt19937_64(seed_seq);
        }
        return std::mt19937_64(*seed);
    }();
    for (std::size_t i = 0; i < Particles::dimension; i++)
    {
        std::uniform_int_distribution<> distrib(lo[i], up[i]);
        for (auto& particle : particles)
            particle.iCell[i] = distrib(gen);
    }
}
template<typename Particles>
void disperse(Particles& particles, std::size_t lo, std::size_t up,
              std::optional<int> seed = std::nullopt)
{
    auto constexpr static dim = Particles::dimension;

    disperse(particles, core::ConstArray<int, dim>(lo), core::ConstArray<int, dim>(up), seed);
}

template<std::size_t dim, typename Box, typename ParticleArray_t = PHARE::core::ParticleArray<dim>>
auto make_particles(std::size_t ppc, Box disperse_in, std::optional<int> seed = std::nullopt)
{
    auto particles = ParticleArray_t{disperse_in};
    particles.vector().resize(ppc * disperse_in.size(), particle<dim>());
    disperse(particles, disperse_in.lower, disperse_in.upper, seed);
    return particles;
}



} // namespace PHARE::core::bench


namespace PHARE::core
{
template<typename Box_t, typename RValue = std::uint32_t>
class LocalisedCellFlattener
{
public:
    static constexpr std::size_t dim = Box_t::dimension;

    LocalisedCellFlattener(Box_t const& b)
        : box{b}
        , shape{box.shape().toArray()}
    {
    }

    RValue operator()(std::array<int, dim> icell) const
    {
        for (std::size_t i = 0; i < dim; ++i)
            icell[i] -= box.lower[i];

        if constexpr (dim == 2)
            return icell[1] + icell[0] * shape[1];
        if constexpr (dim == 3)
            return icell[2] + icell[1] * shape[2] + icell[0] * shape[1] * shape[2];
        return icell[0];
    }
    template<typename Particle>
    RValue operator()(Particle const& particle) const
    {
        return (*this)(particle.iCell);
    }


    Box_t const box;

private:
    std::array<int, dim> const shape;
};
} // namespace PHARE::core


namespace std
{

template<std::size_t dim>
auto& sort(PHARE::core::ParticleArray<dim>& particles)
{
    using box_t = typename PHARE::core::ParticleArray<dim>::box_t;
    PHARE::core::LocalisedCellFlattener<box_t> cell_flattener{grow(particles.box(), 1)};
    std::sort(particles.vector().begin(), particles.vector().end(),
              [&](auto const& a, auto const& b) {
                  return cell_flattener(a.iCell) < cell_flattener(b.iCell);
              });
    return particles;
}

template<std::size_t dim>
auto sort(PHARE::core::ParticleArray<dim>&& particles)
{
    return sort(particles);
}


} // namespace std


#endif /*PHARE_BENCH_CORE_BENCH_H*/
