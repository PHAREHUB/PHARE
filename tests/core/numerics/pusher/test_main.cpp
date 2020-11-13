#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <cstddef>
#include <fstream>
#include <iterator>
#include <memory>
#include <random>
#include <string>
#include <vector>

#include "core/data/particles/particle_array.h"
#include "core/numerics/boundary_condition/boundary_condition.h"
#include "core/numerics/pusher/boris.h"
#include "core/numerics/pusher/pusher_factory.h"
#include "core/utilities/range/range.h"
#include "core/utilities/box/box.h"

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
class Interpolator
{
public:
    template<typename PartIterator, typename Electromag, typename GridLayout>
    void operator()(PartIterator begin, PartIterator end, Electromag const&, GridLayout&)
    {
        for (auto currPart = begin; currPart != end; ++currPart)
        {
            currPart->Ex = 0.01;
            currPart->Ey = -0.05;
            currPart->Ez = 0.05;
            currPart->Bx = 1.;
            currPart->By = 1.;
            currPart->Bz = 1.;
        }
    }
};


// mock of electromag just so that the Pusher gives something to
// the Interpolator
class Electromag
{
};


// with this mock, all particles are found inside
class DummySelector
{
public:
    template<typename Particle>
    bool operator()(Particle) const
    {
        return true;
    }
};


template<std::size_t dimension_>
struct DummyLayout
{
    static constexpr std::size_t dimension = dimension_;
};

template<std::size_t dim>
class APusher : public ::testing::Test
{
public:
    using Pusher_ = BorisPusher<dim, typename ParticleArray<dim>::iterator, Electromag,
                                Interpolator, BoundaryCondition<dim, 1>, DummyLayout<dim>>;

    APusher()
        : expectedTrajectory{readExpectedTrajectory()}
        , particlesIn(1)
        , particlesOut(1)
        , pusher{std::make_unique<Pusher_>()}
        , mass{1}
        , dt{0.0001}
        , tstart{0}
        , tend{10}
        , nt{static_cast<std::size_t>((tend - tstart) / dt + 1)}
    {
        particlesIn[0].charge = 1;
        particlesIn[0].v      = {{0, 10., 0}};
        particlesIn[0].iCell.fill(5);
        particlesIn[0].delta.fill(0.0);
        dxyz.fill(0.05);
        for (std::size_t i = 0; i < dim; i++)
            actual[i].resize(nt, 0.05);
        pusher->setMeshAndTimeStep(dxyz, dt);
    }

protected:
    Trajectory expectedTrajectory;
    ParticleArray<dim> particlesIn;
    ParticleArray<dim> particlesOut;
    std::unique_ptr<Pusher_> pusher;
    double mass;
    double dt;
    double tstart, tend;
    std::size_t nt;
    Electromag em;
    Interpolator interpolator;
    DummySelector selector;
    DummyLayout<dim> layout;
    // BoundaryCondition bc;

    std::array<std::vector<float>, dim> actual;
    std::array<double, dim> dxyz;
};


using APusher1D = APusher<1>;
using APusher2D = APusher<2>;
using APusher3D = APusher<3>;

TEST_F(APusher3D, trajectoryIsOk)
{
    auto rangeIn  = makeRange(std::begin(particlesIn), std::end(particlesIn));
    auto rangeOut = makeRange(std::begin(particlesOut), std::end(particlesOut));
    std::copy(rangeIn.begin(), rangeIn.end(), rangeOut.begin());

    for (decltype(nt) i = 0; i < nt; ++i)
    {
        actual[0][i]
            = (particlesOut[0].iCell[0] + particlesOut[0].delta[0]) * static_cast<float>(dxyz[0]);
        actual[1][i]
            = (particlesOut[0].iCell[1] + particlesOut[0].delta[1]) * static_cast<float>(dxyz[1]);
        actual[2][i]
            = (particlesOut[0].iCell[2] + particlesOut[0].delta[2]) * static_cast<float>(dxyz[2]);

        pusher->move(rangeIn, rangeOut, em, mass, interpolator, selector, layout);

        std::copy(rangeOut.begin(), rangeOut.end(), rangeIn.begin());
    }

    EXPECT_THAT(actual[0], ::testing::Pointwise(::testing::DoubleNear(1e-5), expectedTrajectory.x));
    EXPECT_THAT(actual[1], ::testing::Pointwise(::testing::DoubleNear(1e-5), expectedTrajectory.y));
    EXPECT_THAT(actual[2], ::testing::Pointwise(::testing::DoubleNear(1e-5), expectedTrajectory.z));
}




TEST_F(APusher2D, trajectoryIsOk)
{
    auto rangeIn  = makeRange(std::begin(particlesIn), std::end(particlesIn));
    auto rangeOut = makeRange(std::begin(particlesOut), std::end(particlesOut));
    std::copy(rangeIn.begin(), rangeIn.end(), rangeOut.begin());

    for (decltype(nt) i = 0; i < nt; ++i)
    {
        actual[0][i]
            = (particlesOut[0].iCell[0] + particlesOut[0].delta[0]) * static_cast<float>(dxyz[0]);
        actual[1][i]
            = (particlesOut[0].iCell[1] + particlesOut[0].delta[1]) * static_cast<float>(dxyz[1]);

        pusher->move(rangeIn, rangeOut, em, mass, interpolator, selector, layout);

        std::copy(rangeOut.begin(), rangeOut.end(), rangeIn.begin());
    }

    EXPECT_THAT(actual[0], ::testing::Pointwise(::testing::DoubleNear(1e-5), expectedTrajectory.x));
    EXPECT_THAT(actual[1], ::testing::Pointwise(::testing::DoubleNear(1e-5), expectedTrajectory.y));
}



TEST_F(APusher1D, trajectoryIsOk)
{
    auto rangeIn  = makeRange(std::begin(particlesIn), std::end(particlesIn));
    auto rangeOut = makeRange(std::begin(particlesOut), std::end(particlesOut));
    std::copy(rangeIn.begin(), rangeIn.end(), rangeOut.begin());

    for (decltype(nt) i = 0; i < nt; ++i)
    {
        actual[0][i]
            = (particlesOut[0].iCell[0] + particlesOut[0].delta[0]) * static_cast<float>(dxyz[0]);

        pusher->move(rangeIn, rangeOut, em, mass, interpolator, selector, layout);

        std::copy(rangeOut.begin(), rangeOut.end(), rangeIn.begin());
    }

    EXPECT_THAT(actual[0], ::testing::Pointwise(::testing::DoubleNear(1e-5), expectedTrajectory.x));
}



// the idea of this test is to create a 1D domain [0,1[, push the particles
// until the newEnd returned by the pusher is != the original end, which means
// some particles are out. Then we test the properties of the particles that leave
// and those that stay.
class APusherWithLeavingParticles : public ::testing::Test
{
public:
    APusherWithLeavingParticles()
        : particlesIn(1000)
        , particlesOut1(1000)
        , particlesOut2(1000)
        , pusher{std::make_unique<
              BorisPusher<1, ParticleArray<1>::iterator, Electromag, Interpolator,
                          BoundaryCondition<1, 1>, DummyLayout<1>>>()}
        , mass{1}
        , dt{0.001}
        , tstart{0}
        , tend{10}
        , nt{static_cast<std::size_t>((tend - tstart) / dt + 1)}
        , domain{Point<double, 1>{0.}, Point<double, 1>{1.}}
        , cells{Point{0}, Point{9}}
    {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> dis(0, 9);
        std::uniform_real_distribution<float> delta(0, 1);

        for (auto& part : particlesIn)
        {
            part.charge = 1;
            part.v      = {{5., 0., 0.}};
            part.iCell  = {{dis(gen)}};
            part.delta  = {{delta(gen)}};
        }
        pusher->setMeshAndTimeStep({{dx}}, dt);
    }


protected:
    ParticleArray<1> particlesIn;
    ParticleArray<1> particlesOut1;
    ParticleArray<1> particlesOut2;
    std::unique_ptr<BorisPusher<1, ParticleArray<1>::iterator, Electromag, Interpolator,
                                BoundaryCondition<1, 1>, DummyLayout<1>>>
        pusher;
    double mass;
    double dt;
    double tstart;
    double tend;
    std::size_t nt;
    Electromag em;
    Interpolator interpolator;
    double dx = 0.1;
    Box<double, 1> domain;
    Box<int, 1> cells;
    BoundaryCondition<1, 1> bc;
};




TEST_F(APusherWithLeavingParticles, splitLeavingFromNonLeavingParticles)
{
    auto rangeIn = makeRange(std::begin(particlesIn), std::end(particlesIn));
    auto newEnd  = std::end(particlesIn);
    auto selector
        = [this](Particle<1> const& part) { return PHARE::core::isIn(cellAsPoint(part), cells); };


    for (decltype(nt) i = 0; i < nt; ++i)
    {
        auto layout = DummyLayout<1>{};
        newEnd      = pusher->move(rangeIn, rangeIn, em, mass, interpolator, selector, layout);

        if (newEnd != std::end(particlesIn))
        {
            std::cout << "stopping integration at i = " << i << "\n";
            std::cout << std::distance(std::begin(particlesIn), newEnd) << " in domain\n";
            std::cout << std::distance(newEnd, std::end(particlesIn)) << " leaving\n";
            break;
        }
    }
    EXPECT_TRUE(std::none_of(newEnd, std::end(particlesIn), selector));
    EXPECT_TRUE(std::all_of(std::begin(particlesIn), newEnd, selector));
}



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




TEST_F(APusherWithLeavingParticles, pusherWithOrWithoutBCReturnsReturnEqualStayingParticles)
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
        = [this](auto const& part) { return PHARE::core::isIn(cellAsPoint(part), cells); };

    for (decltype(nt) i = 0; i < nt; ++i)
    {
        auto layout = DummyLayout<1>{};
        newEndWithBC
            = pusher->move(rangeIn, rangeOut1, em, mass, interpolator, selector, bc, layout);
        newEndWithoutBC
            = pusher->move(rangeIn, rangeOut2, em, mass, interpolator, selector, layout);
        [[maybe_unused]] auto s2 = rangeOut2.size();
        [[maybe_unused]] auto s1 = rangeOut1.size();
        [[maybe_unused]] auto s  = rangeIn.size();

        if (newEndWithBC != std::end(particlesOut1) || newEndWithoutBC != std::end(particlesOut2))
        {
            std::cout << "stopping integration at i = " << i << "\n";
            std::cout << std::distance(std::begin(particlesOut1), newEndWithBC) << " in domain\n";
            std::cout << std::distance(newEndWithBC, std::end(particlesOut1)) << " leaving\n";
            break;
        }
    }
    auto part1 = std::begin(particlesOut1);
    auto part2 = std::begin(particlesOut2);

    for (; part1 < newEndWithBC && part2 < newEndWithoutBC; ++part1, ++part2)
    {
        EXPECT_FLOAT_EQ(part1->delta[0], part2->delta[0]);
        EXPECT_FLOAT_EQ(part1->delta[1], part2->delta[1]);
        EXPECT_FLOAT_EQ(part1->delta[2], part2->delta[2]);

        EXPECT_DOUBLE_EQ(part1->v[0], part2->v[0]);
        EXPECT_DOUBLE_EQ(part1->v[1], part2->v[1]);
        EXPECT_DOUBLE_EQ(part1->v[2], part2->v[2]);
    }
}



TEST(APusherFactory, canReturnABorisPusher)
{
    auto pusher
        = PusherFactory::makePusher<1, ParticleArray<1>::iterator, Electromag, Interpolator,
                                    BoundaryCondition<1, 1>, DummyLayout<1>>("modified_boris");

    EXPECT_NE(nullptr, pusher);
}



int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
