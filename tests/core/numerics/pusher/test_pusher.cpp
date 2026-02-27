#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <array>
#include <core/data/electromag/electromag.hpp>
#include <core/data/grid/gridlayout.hpp>
#include <core/numerics/interpolator/interpolator.hpp>
#include <core/numerics/ion_updater/ion_updater.hpp>
#include <cstddef>
#include <fstream>
#include <iterator>
#include <memory>
#include <random>
#include <string>
#include <vector>

#include "core/data/particles/particle_array.hpp"
#include "core/numerics/boundary_condition/boundary_condition.hpp"
#include "core/numerics/pusher/boris.hpp"
#include "core/utilities/range/range.hpp"
#include "core/utilities/box/box.hpp"

#include "tests/core/data/gridlayout/test_gridlayout.hpp"
#include "tests/core/data/ion_population/test_ion_population_fixtures.hpp"

using namespace PHARE::core;



struct Trajectory
{
    std::vector<float> x, y, z;
    std::vector<float> vx, vy, vz;

    Trajectory(std::size_t size)
        : x(size)
        , y(size)
        , z(size)
        , vx(size)
        , vy(size)
        , vz(size)
    {
    }
};


Trajectory readExpectedTrajectory()
{
    std::ifstream file{"pusher_test_in.txt"};

    std::vector<double> data;

    std::copy(std::istream_iterator<double>(file), std::istream_iterator<double>(),
              std::back_inserter(data));

    auto nlines = data.size() / 6;
    Trajectory traj(nlines);

    for (auto iline = 0u; iline < nlines; ++iline)
    {
        traj.x[iline]  = static_cast<float>(data[iline * 6]);
        traj.y[iline]  = static_cast<float>(data[iline * 6 + 1]);
        traj.z[iline]  = static_cast<float>(data[iline * 6 + 2]);
        traj.vx[iline] = static_cast<float>(data[iline * 6 + 3]);
        traj.vy[iline] = static_cast<float>(data[iline * 6 + 4]);
        traj.vz[iline] = static_cast<float>(data[iline * 6 + 5]);
    }

    return traj;
}


// this is a mock of a true Interpolator
// we hard-code the fields that the particles will see to the value used
// in the python script that generates a trajectory
// this way, we don't need a proper Electromag, VecField, Field etc. objects
// to test the pusher.
class TestInterpolator
{
    using E_B_tuple = std::tuple<std::array<double, 3>, std::array<double, 3>>;

public:
    template<typename Particle_t, typename Electromag, typename GridLayout>
    auto operator()(Particle_t& particle, Electromag const&, GridLayout&)
    {
        E_B_tuple eb_interop;
        auto& [pE, pB]        = eb_interop;
        auto& [pEx, pEy, pEz] = pE;
        auto& [pBx, pBy, pBz] = pB;

        pEx = 0.01;
        pEy = -0.05;
        pEz = 0.05;
        pBx = 1.;
        pBy = 1.;
        pBz = 1.;

        return eb_interop;
    }

    void particleToMesh(auto&&... args) {}
};


// mock of electromag just so that the Pusher gives something to
// the Interpolator
class TestElectromag
{
};

template<typename GridLayout, typename ParticleArray_>
struct TestIons // Some tests have large domains but no need for fields
{
    using particle_array_type = ParticleArray_;

    struct TestIonPop
    {
        auto& domainParticles() { return domain; }
        auto& levelGhostParticles() { return levelGhost; }
        auto& patchGhostParticles() { return patchGhost; }
        double mass() const { return 1; }

        GridLayout const& layout;
        ParticleArray_ domain{layout.AMRBox()};
        ParticleArray_ patchGhost{layout.AMRBox()};
        ParticleArray_ levelGhost{layout.AMRBox()};
    };

    GridLayout layout;
    TestIonPop pop{layout};
};


// with this mock, all particles are found inside
class DummySelector
{
public:
    template<typename Range>
    Range operator()(Range& particles) const
    {
        return particles;
    }
};



template<std::size_t dim>
class APusher : public ::testing::Test
{
    static constexpr auto dimension    = dim;
    static constexpr auto interp_order = 1;
    constexpr static PHARE::SimOpts opts{dimension, interp_order};

    using PHARE_TYPES  = PHARE::core::PHARE_Types<opts>;
    using Interpolator = TestInterpolator;
    using Electromag   = TestElectromag;
    using GridLayout_t = TestGridLayout<typename PHARE_TYPES::GridLayout_t>;
    using Ions_t       = TestIons<GridLayout_t, ParticleArray<dim>>;
    using IonUpdater   = typename PHARE::core::IonUpdater<Ions_t, Electromag, GridLayout_t>;
    using Boxing_t     = PHARE::core::UpdaterSelectionBoxing<IonUpdater, GridLayout_t>;

public:
    using Pusher_ = BorisPusher<dim, ParticleArray<dim>, Electromag, Interpolator,
                                BoundaryCondition<dim, 1>, GridLayout_t>;

    APusher()
        : expectedTrajectory{readExpectedTrajectory()}
        , layout{30}
        , pusher{std::make_unique<Pusher_>()}
        , mass{1}
        , dt{0.0001}
        , tstart{0}
        , tend{10}
        , nt{static_cast<std::size_t>((tend - tstart) / dt + 1)}

    {
        particles.emplace_back(
            Particle{1., 1., ConstArray<int, dim>(5), ConstArray<double, dim>(0.), {0., 10., 0.}});
        dxyz.fill(0.05);
        for (std::size_t i = 0; i < dim; i++)
            actual[i].resize(nt, 0.05);
        pusher->setMeshAndTimeStep(dxyz, dt);

        std::transform(std::begin(dxyz), std::end(dxyz), std::begin(halfDtOverDl),
                       [dt = this->dt](double& x) { return 0.5 * dt / x; });
    }

protected:
    using Particle = typename ParticleArray<dim>::Particle_t;
    Trajectory expectedTrajectory;
    GridLayout_t layout;

    std::unique_ptr<Pusher_> pusher;
    double mass;
    double dt;
    double tstart, tend;
    std::size_t nt;
    Electromag em;
    Interpolator interpolator;
    // BoundaryCondition bc;

    std::array<std::vector<float>, dim> actual;
    std::array<double, dim> dxyz;
    Ions_t ions{layout};

    ParticleArray<dim>& particles = ions.pop.domainParticles();

    std::array<double, dim> halfDtOverDl;
    double const dto2m = 0.5 * dt / mass;

    void move()
    {
        for (auto& particle : particles)
        {
            particle.iCell = boris::advance(particle, halfDtOverDl);
            boris::accelerate(particle, interpolator(particle, em, layout), dto2m);
            particle.iCell = boris::advance(particle, halfDtOverDl);
        }
    }
};


using APusher1D = APusher<1>;
using APusher2D = APusher<2>;
using APusher3D = APusher<3>;

TEST_F(APusher3D, trajectoryIsOk)
{
    for (decltype(nt) i = 0; i < nt; ++i)
    {
        actual[0][i] = (particles[0].iCell[0] + particles[0].delta[0]) * dxyz[0];
        actual[1][i] = (particles[0].iCell[1] + particles[0].delta[1]) * dxyz[1];
        actual[2][i] = (particles[0].iCell[2] + particles[0].delta[2]) * dxyz[2];

        move();
    }

    EXPECT_THAT(actual[0], ::testing::Pointwise(::testing::DoubleNear(1e-5), expectedTrajectory.x));
    EXPECT_THAT(actual[1], ::testing::Pointwise(::testing::DoubleNear(1e-5), expectedTrajectory.y));
    EXPECT_THAT(actual[2], ::testing::Pointwise(::testing::DoubleNear(1e-5), expectedTrajectory.z));
}




TEST_F(APusher2D, trajectoryIsOk)
{
    for (decltype(nt) i = 0; i < nt; ++i)
    {
        actual[0][i] = (particles[0].iCell[0] + particles[0].delta[0]) * dxyz[0];
        actual[1][i] = (particles[0].iCell[1] + particles[0].delta[1]) * dxyz[1];

        move();
    }

    EXPECT_THAT(actual[0], ::testing::Pointwise(::testing::DoubleNear(1e-5), expectedTrajectory.x));
    EXPECT_THAT(actual[1], ::testing::Pointwise(::testing::DoubleNear(1e-5), expectedTrajectory.y));
}



TEST_F(APusher1D, trajectoryIsOk)
{
    for (decltype(nt) i = 0; i < nt; ++i)
    {
        actual[0][i] = (particles[0].iCell[0] + particles[0].delta[0]) * dxyz[0];

        move();
    }

    EXPECT_THAT(actual[0], ::testing::Pointwise(::testing::DoubleNear(1e-5), expectedTrajectory.x));
}



// the idea of this test is to create a 1D domain [0,1[, push the particles
// until the newEnd returned by the pusher is != the original end, which means
// some particles are out. Then we test the properties of the particles that leave
// and those that stay.
class APusherWithLeavingParticles : public ::testing::Test
{
    using Interpolator = TestInterpolator;
    using Electromag   = TestElectromag;

    static constexpr auto dimension    = 1;
    static constexpr auto interp_order = 1;
    constexpr static PHARE::SimOpts opts{dimension, interp_order};

    using PHARE_TYPES     = PHARE::core::PHARE_Types<opts>;
    using ParticleArray_t = ParticleArray<1>;
    using GridLayout_t    = TestGridLayout<typename PHARE_TYPES::GridLayout_t>;
    using Ions_t          = PHARE_TYPES::Ions_t;
    using IonUpdater      = typename PHARE::core::IonUpdater<Ions_t, Electromag, GridLayout_t>;
    using Boxing_t        = PHARE::core::UpdaterSelectionBoxing<IonUpdater, GridLayout_t>;

public:
    APusherWithLeavingParticles()
        : pusher{std::make_unique<BorisPusher<1, ParticleArray<1>, TestElectromag, TestInterpolator,
                                              BoundaryCondition<1, 1>, GridLayout_t>>()}
        , mass{1}
        , dt{0.001}
        , tstart{0}
        , tend{10}
        , nt{static_cast<std::size_t>((tend - tstart) / dt + 1)}
        , layout{10}
    // , particlesOut1{grow(domain, 1), 1000}
    // , particlesOut2{grow(domain, 1), 1000}
    {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> dis(0, 9);
        std::uniform_real_distribution<double> delta(0, 1);

        for (std::size_t iPart = 0; iPart < 1000; ++iPart)
        {
            particlesIn.emplace_back(
                Particle<1>{1., 1., {{dis(gen)}}, {{delta(gen)}}, {{5., 0., 0.}}});
        }
        pusher->setMeshAndTimeStep({{dx}}, dt);
    }


protected:
    std::unique_ptr<BorisPusher<1, ParticleArray<1>, Electromag, Interpolator,
                                BoundaryCondition<1, 1>, GridLayout_t>>
        pusher;
    double mass;
    double dt;
    double tstart;
    double tend;
    std::size_t nt;
    Electromag em;
    Interpolator interpolator;
    double dx = 0.1;
    GridLayout_t layout;
    Box<int, 1> domain = layout.AMRBox();
    // BoundaryCondition<1, 1> bc;
    // ParticleArray<1> particlesOut1;
    // ParticleArray<1> particlesOut2;

    // TestIons<GridLayout_t, ParticleArray<dimension>> ions{layout};
    UsableIons<ParticleArray_t> ions{layout};

    ParticleArray<dimension>& particlesIn = ions[0].domainParticles();

    Boxing_t const boxing{layout, {layout.AMRBox()}};
};




TEST_F(APusherWithLeavingParticles, splitLeavingFromNonLeavingParticles)
{
    for (decltype(nt) i = 0; i < nt; ++i)
        pusher->move(ions[0], em, interpolator, boxing);

    auto& patchGhost = ions[0].patchGhostParticles();
    EXPECT_TRUE(
        std::none_of(std::begin(patchGhost), std::end(patchGhost),
                     [this](auto const& particle) { return PHARE::core::isIn(particle, domain); }));
    EXPECT_TRUE(
        std::all_of(std::begin(particlesIn), std::end(particlesIn),
                    [this](auto const& particle) { return PHARE::core::isIn(particle, domain); }));
}


// removed boundary condition partitioner, fix that when BCs are implemented
#if 0
TEST_F(APusherWithLeavingParticles, pusherWithOrWithoutBCReturnsSameNbrOfStayingParticles)
{
    auto rangeIn   = makeRange(std::begin(particlesIn), std::end(particlesIn));
    auto rangeOut1 = makeRange(std::begin(particlesOut1), std::end(particlesOut1));
    auto rangeOut2 = makeRange(std::begin(particlesOut2), std::end(particlesOut2));
    std::copy(rangeIn.begin(), rangeIn.end(), rangeOut1.begin());
    std::copy(rangeIn.begin(), rangeIn.end(), rangeOut2.begin());

    auto newEndWithBC    = std::end(particlesOut1);
    auto newEndWithoutBC = std::end(particlesOut2);

    bc.setBoundaryBoxes(std::vector<Box<int, 1>>{});

    auto selector
        = [this](Particle<1> const& part) { return PHARE::core::isIn(cellAsPoint(part), cells); };

    for (decltype(nt) i = 0; i < nt; ++i)
    {
        auto layout = DummyLayout<1>{};
        newEndWithBC
            = pusher->move(rangeIn, rangeOut1, em, mass, interpolator, selector, bc, layout);
        newEndWithoutBC
            = pusher->move(rangeIn, rangeOut2, em, mass, interpolator, selector, layout);

        if (newEndWithBC != std::end(particlesOut1) || newEndWithoutBC != std::end(particlesOut2))
        {
            std::cout << "stopping integration at i = " << i << "\n";
            std::cout << std::distance(std::begin(particlesOut1), newEndWithBC) << " in domain\n";
            std::cout << std::distance(newEndWithBC, std::end(particlesOut1)) << " leaving\n";
            break;
        }
    }
    auto d1 = std::distance(std::begin(particlesOut1), newEndWithBC);
    auto d2 = std::distance(std::begin(particlesOut2), newEndWithoutBC);
    EXPECT_EQ(d1, d2);
}
#endif




int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
