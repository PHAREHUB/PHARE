#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <cstddef>
#include <fstream>
#include <iterator>
#include <random>
#include <string>
#include <vector>

#include "data/particles/particle_array.h"
#include "numerics/pusher/boris.h"
#include "utilities/range/range.h"

using namespace PHARE;


struct Trajectory
{
    std::vector<float> x;
    std::vector<float> y;
    std::vector<float> z;
    std::vector<float> vx;
    std::vector<float> vy;
    std::vector<float> vz;

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


Trajectory readFile()
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

class Interpolator
{
public:
    template<typename PartIterator, typename Electromag>
    void operator()(PartIterator begin, PartIterator end, Electromag const& em)
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

class VecField
{
};

class Electromag
{
public:
    VecField E;
    VecField B;
};

class DummySelector
{
public:
    template<typename Particle>
    bool operator()(Particle) const
    {
        return true;
    }
};

class BoundaryCondition
{
public:
    template<typename ParticleIterator>
    ParticleIterator applyOutgoingParticleBC(ParticleIterator begin, ParticleIterator end)
    {
    }
};


class APusher3D : public ::testing::Test
{
public:
    APusher3D()
        : traj{readFile()}
        , particlesIn(1)
        , particlesOut(1)
        , pusher{std::make_unique<BorisPusher<3>>()}
        , mass{1}
        , dt{0.05}
        , tstart{0}
        , tend{10}
        , nt{static_cast<std::size_t>((tend - tstart) / dt + 1)}
    {
        particlesIn[0].charge = 1;
        particlesIn[0].iCell  = {{5, 5, 5}}; // arbitrary we don't care
        particlesIn[0].v      = {{0, 1, 0}};
        particlesIn[0].delta  = {{0.0, 0.0, 0.0}};
    }


protected:
    Trajectory traj;
    ParticleArray<3> particlesIn;
    ParticleArray<3> particlesOut;
    std::unique_ptr<BorisPusher<3>> pusher;
    double mass;
    double dt;
    double tstart;
    double tend;
    std::size_t nt;
    Electromag em;
    Interpolator interpolator;
    DummySelector selector;
    BoundaryCondition bc;
};



TEST_F(APusher3D, trajectoryIsOk)
{
    /*Trajectory traj = readFile();
    std::cout << traj.x[200] << " " << traj.y[200] << " " << traj.z[200] << " " << traj.vx[200]
              << " " << traj.vy[200] << " " << traj.vz[200] << "\n";*/

    auto rangeIn  = makeRange(std::begin(particlesIn), std::end(particlesIn));
    auto rangeOut = makeRange(std::begin(particlesOut), std::end(particlesOut));

    std::vector<float> x(nt);
    std::vector<float> y(nt);
    std::vector<float> z(nt);

    double dx = 0.1;
    double dy = 0.1;
    double dz = 0.1;

    pusher->setMeshAndTimeStep({{dx, dy, dz}}, dt);

    for (decltype(nt) i = 0; i < nt; ++i)
    {
        pusher->move(rangeIn, rangeOut, em, mass, interpolator, selector, bc);

        x[i] = (particlesOut[0].iCell[0] + particlesOut[0].delta[0]) * static_cast<float>(dx);
        y[i] = (particlesOut[0].iCell[1] + particlesOut[0].delta[1]) * static_cast<float>(dy);
        z[i] = (particlesOut[0].iCell[2] + particlesOut[0].delta[2]) * static_cast<float>(dz);
    }

    EXPECT_THAT(traj.x, ::testing::Pointwise(::testing::DoubleNear(1e-5), x));
    EXPECT_THAT(traj.y, ::testing::Pointwise(::testing::DoubleNear(1e-5), y));
    EXPECT_THAT(traj.z, ::testing::Pointwise(::testing::DoubleNear(1e-5), z));
}




int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
