// Accuracy of reduced precision particle storage (see core/utilities/multi_precision.hpp)
//  against plain double storage, for the boris pusher in uniform fields.
//  All arithmetic is in double, only the stored particle attributes lose precision.
//
//  The parametric study table is printed to stdout, run with
//    ./test-particle-precision --gtest_filter='*Study*'

#include "gtest/gtest.h"

#include <array>
#include <cmath>
#include <tuple>
#include <random>
#include <string>
#include <vector>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <algorithm>

#include "core/data/particles/particle.hpp"
#include "core/data/particles/particle_array.hpp"
#include "core/data/particles/particle_storage.hpp"
#include "core/numerics/boundary_condition/boundary_condition.hpp"
#include "core/numerics/pusher/boris.hpp"
#include "core/utilities/box/box.hpp"
#include "core/utilities/types.hpp"
#include "core/utilities/range/range.hpp"

using namespace PHARE::core;



// storage configurations under test

using F64 = double;
template<std::size_t B>
using Fixed = FixedPointUnit<B>;
template<std::size_t B>
using Trunc = TruncatedDouble<B>;

template<std::size_t dim>
using ReferenceParticle = Particle<dim, F64, F64>;

// layout "B": 3d particle in 64 bytes
template<std::size_t dim>
using CompactParticle = Particle<dim, Fixed<6>, Trunc<6>>;



template<typename T>
std::string storage_name()
{
    if constexpr (std::is_same_v<T, double>)
        return "f64";
    else if constexpr (std::is_same_v<T, Fixed<T::bytes>>)
        return "fixed" + std::to_string(8 * T::bytes);
    else
        return "trunc" + std::to_string(8 * T::bytes);
}



// ---------------------------------------------------------------------------------------------
// storage types

TEST(ParticleStorage, compactParticleSizes)
{
    EXPECT_EQ(sizeof(CompactParticle<1>), 48u);
    EXPECT_EQ(sizeof(CompactParticle<2>), 56u);
    EXPECT_EQ(sizeof(CompactParticle<3>), 64u);

    EXPECT_EQ(sizeof(ReferenceParticle<3>), 80u);

    static_assert(sizeof(Fixed<6>) == 6 and alignof(Fixed<6>) == 2);
    static_assert(sizeof(Trunc<6>) == 6 and alignof(Trunc<6>) == 2);
    static_assert(std::is_trivially_copyable_v<CompactParticle<3>>); // SAMRAI streams memcpy
}


template<typename Storage, typename Dist>
void check_roundtrip(Dist&& dist, double const max_error, bool const relative)
{
    std::mt19937_64 gen{1};
    double error_sum = 0, max_seen = 0;
    std::size_t constexpr N = 1'000'000;

    for (std::size_t i = 0; i < N; ++i)
    {
        double const x  = dist(gen);
        Storage const s = x;
        double const d  = s;

        // decoded values are exactly representable, re-encoding is lossless
        ASSERT_EQ(Storage{d}.data, s.data) << x;

        double const error = relative ? (d - x) / std::abs(x) : d - x;
        error_sum += error;
        max_seen = std::max(max_seen, std::abs(error));
    }

    EXPECT_LE(max_seen, max_error);
    // round to nearest, no systematic bias
    EXPECT_LT(std::abs(error_sum / N), max_error * 1e-2);
}

TEST(ParticleStorage, fixedPointRoundtripIsIdempotentAndUnbiased)
{
    std::uniform_real_distribution<double> unit{0, 1};
    // half a step, except within half a step of 1 where the clamp can cost a full step
    check_roundtrip<Fixed<6>>(unit, std::ldexp(1., -48), false);
    check_roundtrip<Fixed<4>>(unit, std::ldexp(1., -32), false);
    check_roundtrip<Fixed<2>>(unit, std::ldexp(1., -16), false);
}

TEST(ParticleStorage, truncatedDoubleRoundtripIsIdempotentAndUnbiased)
{
    std::uniform_real_distribution<double> vel{-10, 10};
    // 1 sign + 11 exponent + (8B - 12) mantissa bits, half ulp relative error
    check_roundtrip<Trunc<7>>(vel, std::ldexp(1., -44), true);
    check_roundtrip<Trunc<6>>(vel, std::ldexp(1., -36), true);
    check_roundtrip<Trunc<4>>(vel, std::ldexp(1., -20), true);
}

TEST(ParticleStorage, fixedPointNeverWrapsToZero)
{
    double const top = static_cast<double>(Fixed<6>{1. - 1e-17});
    EXPECT_LT(top, 1.);
    EXPECT_GT(top, 1. - std::ldexp(1., -47));
    EXPECT_EQ(static_cast<double>(Fixed<6>{1.}), top);
    EXPECT_EQ(static_cast<double>(Fixed<6>{-1e-17}), 0.);
    EXPECT_EQ(static_cast<double>(Fixed<6>{std::nan("")}), 0.);
}

TEST(ParticleStorage, truncatedDoubleSpecialValues)
{
    EXPECT_EQ(static_cast<double>(Trunc<6>{0.}), 0.);
    EXPECT_TRUE(std::signbit(static_cast<double>(Trunc<6>{-0.})));
    EXPECT_EQ(static_cast<double>(Trunc<6>{1.}), 1.);
    EXPECT_EQ(static_cast<double>(Trunc<6>{-2.5}), -2.5);
    // rounding up carries into the exponent
    EXPECT_EQ(static_cast<double>(Trunc<2>{std::nextafter(2., 0.)}), 2.);
}

TEST(ParticleStorage, particleConvertsFromAndComparesWithDoubleView)
{
    std::array<int, 3> cell{1, 2, 3};
    MultiPrecisionArray<double, 3> delta{{.25, .5, .75}}; // exactly representable in all storages
    MultiPrecisionArray<double, 3> v{{1., -2., .5}};
    double weight = 1, charge = 1;

    CompactParticle<3> const particle{weight, charge, cell, delta, v};
    ParticleView<3> const view{weight, charge, cell, delta, v};

    EXPECT_TRUE(particle == view);
    EXPECT_TRUE(view == particle);
    EXPECT_EQ(particle, (CompactParticle<3>{weight, charge, cell, delta, v}));
}



// ---------------------------------------------------------------------------------------------
// boris pusher in uniform fields

struct UniformFields
{
};

struct UniformInterpolator
{
    std::array<double, 3> E, B;

    template<typename Particle_t, typename Electromag, typename GridLayout>
    auto operator()(Particle_t&, Electromag const&, GridLayout&) const
    {
        return std::make_tuple(E, B);
    }
};

struct AllInSelector
{
    template<typename Range>
    Range operator()(Range& particles) const
    {
        return particles;
    }
};

template<std::size_t dim>
struct DummyLayout
{
    static constexpr std::size_t dimension = dim;
    auto AMRBox() const { return Box<int, dimension>{}; }
    auto levelNumber() const { return 0; }
};


struct Setup
{
    std::string name;
    std::array<double, 3> E, B;
    std::size_t n_particles = 64;
    std::size_t n_steps     = 10'000;
    double dt = 0.01, dx = 0.1, mass = 1;
};

// E x B drift + gyration around a tilted B
Setup const drift{"ExB", {0.01, -0.05, 0.05}, {1., 1., 1.}};
// pure gyration, |v| is exactly conserved by boris
Setup const gyration{"gyration", {0., 0., 0.}, {0., 0., 1.}};


template<std::size_t dim>
struct State
{
    std::vector<std::array<double, dim>> x; // iCell + delta, in cells
    std::vector<std::array<double, 3>> v;
};


template<typename Particle_t>
auto push(Setup const& setup)
{
    auto constexpr dim = Particle_t::dimension;
    using Array_t      = ParticleArray<dim, Particle_t>;
    using Range_t      = IndexRange<Array_t>;
    using Pusher_t     = BorisPusher<dim, Range_t, UniformFields, UniformInterpolator,
                                     BoundaryCondition<dim, 1>, DummyLayout<dim>>;

    DummyLayout<dim> layout;
    Array_t particles{layout.AMRBox()};

    // same initial conditions for every storage type
    std::mt19937_64 gen{42};
    std::uniform_real_distribution<double> delta_dist{0, 1}, v_dist{-1, 1};
    for (std::size_t i = 0; i < setup.n_particles; ++i)
    {
        std::array<double, dim> delta;
        std::array<double, 3> v;
        for (auto& d : delta)
            d = delta_dist(gen);
        for (auto& vi : v)
            vi = v_dist(gen);
        particles.push_back(Particle_t{1., 1., ConstArray<int, dim>(0), delta, v});
    }

    Pusher_t pusher;
    pusher.setMeshAndTimeStep(ConstArray<double, dim>(setup.dx), setup.dt);

    UniformFields const em;
    UniformInterpolator interpolator{setup.E, setup.B};
    AllInSelector const selector;
    auto range = makeIndexRange(particles);

    for (std::size_t step = 0; step < setup.n_steps; ++step)
        pusher.move(range, range, em, setup.mass, interpolator, layout, selector, selector);

    State<dim> state;
    for (auto const& particle : particles)
    {
        auto& x = state.x.emplace_back();
        auto& v = state.v.emplace_back();
        for (std::size_t i = 0; i < dim; ++i)
            x[i] = particle.iCell[i] + static_cast<double>(particle.delta[i]);
        for (std::size_t i = 0; i < 3; ++i)
            v[i] = particle.v[i];
    }
    return state;
}


struct StorageErrors
{
    double x     = 0; // max |x - x_ref|, in cells
    double v     = 0; // max |v - v_ref| / |v_ref|
    double speed = 0; // max ||v| - |v_ref|| / |v_ref|
};

template<std::size_t dim>
StorageErrors compare(State<dim> const& state, State<dim> const& ref)
{
    auto norm = [](auto const& a) { return std::sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]); };

    StorageErrors errors;
    for (std::size_t p = 0; p < ref.x.size(); ++p)
    {
        for (std::size_t i = 0; i < dim; ++i)
            errors.x = std::max(errors.x, std::abs(state.x[p][i] - ref.x[p][i]));

        auto const ref_speed = norm(ref.v[p]);
        std::array<double, 3> dv;
        for (std::size_t i = 0; i < 3; ++i)
            dv[i] = state.v[p][i] - ref.v[p][i];
        errors.v     = std::max(errors.v, norm(dv) / ref_speed);
        errors.speed = std::max(errors.speed, std::abs(norm(state.v[p]) - ref_speed) / ref_speed);
    }
    return errors;
}


template<typename Particle_t>
StorageErrors errors_for(Setup const& setup)
{
    auto constexpr dim = Particle_t::dimension;
    return compare(push<Particle_t>(setup), push<ReferenceParticle<dim>>(setup));
}



TEST(ParticleStoragePusher, referenceStorageIsExact)
{
    auto const errors = errors_for<ReferenceParticle<3>>(drift);
    EXPECT_EQ(errors.x, 0.);
    EXPECT_EQ(errors.v, 0.);
}

TEST(ParticleStoragePusher, compactParticleStaysCloseToDouble)
{
    for (auto const& setup : {drift, gyration})
    {
        auto const errors = errors_for<CompactParticle<3>>(setup);
        // velocity errors integrate into position, ~1e-9 * |v| * t / dx here
        EXPECT_LT(errors.x, 1e-6) << setup.name;
        EXPECT_LT(errors.v, 1e-8) << setup.name;
        EXPECT_LT(errors.speed, 1e-8) << setup.name;
    }
}



template<typename... Particles>
void study(Setup const& setup)
{
    std::cout << "\n"
              << setup.name << ": " << setup.n_particles << " particles, " << setup.n_steps
              << " steps, dt " << setup.dt << ", errors vs f64 storage\n";
    std::cout << std::setw(9) << "delta" << std::setw(9) << "v" << std::setw(7) << "bytes"
              << std::setw(12) << "max dx" << std::setw(12) << "max dv/v" << std::setw(12)
              << "max d|v|/v" << "\n";

    auto row = [&]<typename Particle_t>() {
        auto const errors = errors_for<Particle_t>(setup);
        std::cout << std::setw(9) << storage_name<typename Particle_t::delta_type>() << std::setw(9)
                  << storage_name<typename Particle_t::v_type>() << std::setw(7)
                  << sizeof(Particle_t) << std::scientific << std::setprecision(2) << std::setw(12)
                  << errors.x << std::setw(12) << errors.v << std::setw(12) << errors.speed
                  << std::defaultfloat << "\n";
    };
    (row.template operator()<Particles>(), ...);
}

TEST(ParticleStoragePusher, Study)
{
    for (auto const& setup : {drift, gyration})
        study<Particle<3, Fixed<2>, F64>, Particle<3, Fixed<3>, F64>, //
              Particle<3, Fixed<4>, F64>, Particle<3, Fixed<5>, F64>, //
              Particle<3, Fixed<6>, F64>,                             //
              Particle<3, F64, Trunc<3>>, Particle<3, F64, Trunc<4>>, //
              Particle<3, F64, Trunc<5>>, Particle<3, F64, Trunc<6>>, //
              Particle<3, F64, Trunc<7>>,                             //
              Particle<3, Fixed<4>, Trunc<6>>, CompactParticle<3>,    //
              Particle<3, Fixed<6>, Trunc<7>>>(setup);
}



int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
