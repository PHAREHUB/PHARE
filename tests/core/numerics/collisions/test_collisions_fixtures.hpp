#ifndef PHARE_TEST_CORE_NUMERICS_COLLISIONS_FIXTURES_HPP
#define PHARE_TEST_CORE_NUMERICS_COLLISIONS_FIXTURES_HPP

#include "core/data/particles/particle_array.hpp"
#include "core/numerics/collisions/collisions.hpp"
#include "initializer/data_provider.hpp"

#include <array>
#include <cmath>
#include <random>
#include <string>
#include <vector>


namespace PHARE::core
{
// the population interface Collisions uses: mass, charge and domain particles
struct TestPopulation
{
    using ParticleArray_t = ParticleArray<1>;

    TestPopulation(double const mass, double const charge, int const nbrCells)
        : mass_{mass}
        , particles_{Box<int, 1>{{0}, {nbrCells - 1}}}
    {
        info_["charge"] = charge;
    }

    NO_DISCARD double mass() const { return mass_; }
    NO_DISCARD auto const& particleInitializerInfo() const { return info_; }
    NO_DISCARD auto& domainParticles() { return particles_; }
    NO_DISCARD auto const& domainParticles() const { return particles_; }

    // `ppc` particles per cell of density `density`, drifting maxwellian of thermal speeds vth
    void load(std::size_t const ppc, double const density, std::array<double, 3> const& drift,
              std::array<double, 3> const& vth, std::mt19937_64& rng, bool const randomWeights)
    {
        std::uniform_real_distribution<double> uniform{0., 1.};
        std::normal_distribution<double> normal{0., 1.};
        auto const nbrCells = particles_.box().upper[0] + 1;

        for (int cell = 0; cell < nbrCells; ++cell)
        {
            std::vector<double> weights(ppc, 1.);
            if (randomWeights) // log-uniform over two decades
                for (auto& w : weights)
                    w = std::pow(10., 2. * uniform(rng));
            double total = 0;
            for (auto const w : weights)
                total += w;

            for (auto const w : weights)
            {
                Particle<1> particle;
                particle.weight = w * density / total;
                particle.iCell  = {cell};
                particle.delta  = {uniform(rng)};
                for (std::size_t i = 0; i < 3; ++i)
                    particle.v[i] = drift[i] + vth[i] * normal(rng);
                particles_.push_back(particle);
            }
        }
    }

    double mass_;
    initializer::PHAREDict info_;
    ParticleArray_t particles_;
};

struct TestIons
{
    using particle_array_type = ParticleArray<1>;

    NO_DISCARD auto begin() { return populations.begin(); }
    NO_DISCARD auto end() { return populations.end(); }
    NO_DISCARD auto begin() const { return populations.begin(); }
    NO_DISCARD auto end() const { return populations.end(); }

    std::vector<TestPopulation> populations;
};


// density weighted moments of one population, over all its particles
struct Moments
{
    double density = 0;
    std::array<double, 3> velocity{};
    std::array<double, 3> temperature{}; // m <(v - V)^2> per direction

    Moments(TestPopulation const& pop, std::size_t const nbrCells)
    {
        double w = 0;
        std::array<double, 3> sum{}, sq{};
        for (auto const& p : pop.domainParticles())
        {
            w += p.weight;
            for (std::size_t i = 0; i < 3; ++i)
            {
                sum[i] += p.weight * p.v[i];
                sq[i] += p.weight * p.v[i] * p.v[i];
            }
        }
        density = w / nbrCells;
        for (std::size_t i = 0; i < 3; ++i)
        {
            velocity[i]    = sum[i] / w;
            temperature[i] = pop.mass() * (sq[i] / w - velocity[i] * velocity[i]);
        }
    }

    NO_DISCARD double parallel() const { return temperature[0]; }
    NO_DISCARD double perpendicular() const { return 0.5 * (temperature[1] + temperature[2]); }
    NO_DISCARD double isotropic() const
    {
        return (temperature[0] + temperature[1] + temperature[2]) / 3.;
    }
};


// classical 4th order Runge-Kutta for y' = f(y), y an array
template<std::size_t N>
auto rk4(std::array<double, N> y, double const dt, auto&& f)
{
    auto axpy = [](auto const& a, double const s, auto const& b) {
        auto r = a;
        for (std::size_t i = 0; i < N; ++i)
            r[i] += s * b[i];
        return r;
    };
    auto const k1 = f(y);
    auto const k2 = f(axpy(y, 0.5 * dt, k1));
    auto const k3 = f(axpy(y, 0.5 * dt, k2));
    auto const k4 = f(axpy(y, dt, k3));
    for (std::size_t i = 0; i < N; ++i)
        y[i] += dt / 6. * (k1[i] + 2. * k2[i] + 2. * k3[i] + k4[i]);
    return y;
}

} // namespace PHARE::core

#endif
