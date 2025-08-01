#ifndef PHARE_TEST_MAXWELLIAN_PARTICLE_INITIALIZER_TEST_INIT_FUNCS_HPP
#define PHARE_TEST_MAXWELLIAN_PARTICLE_INITIALIZER_TEST_INIT_FUNCS_HPP


#include <memory>
#include <vector>
#include <optional>

#include "test_funcs.hpp"
#include "core/utilities/types.hpp"
#include "core/data/ions/particle_initializers/maxwellian_particle_initializer.hpp"

// These functions were extracted from the MaxwellianParticleInitializer when
//  the ScalarFunctions were replaced


namespace PHARE::core::test_fn::func_1d
{
template<typename ParticleArray, typename GridLayout>
class TestMaxwellianParticleInitializer
{
    constexpr static std::size_t dimension           = 1;
    constexpr static std::size_t nbrParticlePerCell_ = 10;

    using MaxInit = MaxwellianParticleInitializer<ParticleArray, GridLayout>;

public:
    TestMaxwellianParticleInitializer(ParticleArray& particles, GridLayout const& layout);

private:
    Basis basis_           = Basis::Cartesian;
    double particleCharge_ = 1.1;
    std::optional<std::size_t> rngSeed_{1337};
};

template<typename ParticleArray, typename GridLayout>
TestMaxwellianParticleInitializer<ParticleArray, GridLayout> //
    ::TestMaxwellianParticleInitializer(ParticleArray& particles, GridLayout const& layout)
{
    auto const meshSize = layout.meshSize();
    auto const& dx      = meshSize[0];

    /* get indices start and stop. we take primal/primal/primal because
       that is what GridLayout::cellCenteredCoordinate() requires */
    auto const [ix0, ix1] = layout.physicalStartToEnd(QtyCentering::primal, Direction::X);

    double cellVolume = dx;

    // random seed and generator needed to load maxwellian velocity
    // and random position with the cell
    auto generator = MaxInit::getRNG(rngSeed_);

    // beware: we're looping over the cell but use primal indices because of
    // GridLayout::cellCenteredCoordinates
    // therefore i(x,y,z)1 must be excluded.

    // gab references for convenience
    // auto& density         = *density_;
    // auto& bulkVelocity    = *bulkVelocity_;
    // auto& thermalVelocity = *thermalVelocity_;

    for (std::uint32_t ix = ix0; ix < ix1; ++ix)
    {
        double n; // cell centered density
        std::array<double, 3> particleVelocity;
        std::array<std::array<double, 3>, 3> basis;

        // get the coordinate of the current cell
        auto coord = layout.cellCenteredCoordinates(ix);
        auto x     = coord[0];

        // now get density, velocity and thermal speed values
        n       = density_(x);
        auto Vx = vx(x);
        auto Vy = vy(x);
        auto Vz = vz(x);

        auto Vthx = vthx(x);
        auto Vthy = vthy(x);
        auto Vthz = vthz(x);

        // weight for all particles in this cell
        auto cellWeight = n * cellVolume / nbrParticlePerCell_;

        ParticleDeltaDistribution randPosX;

        if (basis_ == Basis::Magnetic)
        {
            auto Bx = bx(x);
            auto By = by(x);
            auto Bz = bz(x);

            localMagneticBasis({Bx, By, Bz}, basis);
        }

        for (std::uint32_t ipart = 0; ipart < nbrParticlePerCell_; ++ipart)
        {
            maxwellianVelocity({Vx, Vy, Vz}, {Vthx, Vthy, Vthz}, generator, particleVelocity);

            if (basis_ == Basis::Magnetic)
            {
                particleVelocity = basisTransform(basis, particleVelocity);
            }

            std::array<float, dimension> delta = {{randPosX(generator)}};

            Particle<dimension> tmpParticle;

            // particle iCell is in AMR index
            auto AMRCellIndex = layout.localToAMR(Point{ix});

            tmpParticle.weight = cellWeight;
            tmpParticle.charge = particleCharge_;
            tmpParticle.iCell  = AMRCellIndex.template toArray<int>();
            tmpParticle.delta  = delta;
            tmpParticle.v      = particleVelocity;

            particles.push_back(std::move(tmpParticle));
        }
    }
}
} // namespace PHARE::core::test_fn::func_1d

namespace PHARE::core::test_fn::func_2d
{
template<typename ParticleArray, typename GridLayout>
class TestMaxwellianParticleInitializer
{
    constexpr static std::size_t dimension           = 2;
    constexpr static std::size_t nbrParticlePerCell_ = 10;

    using MaxInit = MaxwellianParticleInitializer<ParticleArray, GridLayout>;

public:
    TestMaxwellianParticleInitializer(ParticleArray& particles, GridLayout const& layout);

private:
    Basis basis_           = Basis::Cartesian;
    double particleCharge_ = 1.1;
    std::optional<std::size_t> rngSeed_{1337};
};

template<typename ParticleArray, typename GridLayout>
TestMaxwellianParticleInitializer<ParticleArray, GridLayout>::TestMaxwellianParticleInitializer(
    ParticleArray& particles, GridLayout const& layout)
{
    auto const meshSize = layout.meshSize();
    auto const& dx      = meshSize[0];
    auto const& dy      = meshSize[1];

    /* get indices start and stop. we take primal/primal/primal because
       that is what GridLayout::cellCenteredCoordinate() requires */
    auto const [ix0, ix1] = layout.physicalStartToEnd(QtyCentering::primal, Direction::X);
    auto const [iy0, iy1] = layout.physicalStartToEnd(QtyCentering::primal, Direction::Y);

    auto cellVolume = dx * dy;

    // random seed and generator needed to load maxwellian velocity
    // and random position with the cell
    auto generator = MaxInit::getRNG(rngSeed_);

    // beware: we're looping over the cell but use primal indices because of
    // GridLayout::cellCenteredCoordinates
    // therefore i(x,y,z)1 must be excluded.
    // gab references for convenience
    // auto& density         = *density_;
    // auto& bulkVelocity    = *bulkVelocity_;
    // auto& thermalVelocity = *thermalVelocity_;


    for (std::uint32_t ix = ix0; ix < ix1; ++ix)
    {
        for (std::uint32_t iy = iy0; iy < iy1; ++iy)
        {
            double n; // cell centered density
            std::array<double, 3> particleVelocity;
            std::array<std::array<double, 3>, 3> basis;

            // get the coordinate of the current cell
            auto coord = layout.cellCenteredCoordinates(ix, iy);
            auto x     = coord[0];
            auto y     = coord[1];

            // now get density, velocity and thermal speed values
            n         = density_(x, y);
            auto Vx   = vx(x, y);
            auto Vy   = vy(x, y);
            auto Vz   = vz(x, y);
            auto Vthx = vthx(x, y);
            auto Vthy = vthy(x, y);
            auto Vthz = vthz(x, y);

            // weight for all particles in this cell
            auto cellWeight = n * cellVolume / nbrParticlePerCell_;

            ParticleDeltaDistribution randPosX;
            ParticleDeltaDistribution randPosY;

            if (basis_ == Basis::Magnetic)
            {
                auto Bx = bx(x, y);
                auto By = by(x, y);
                auto Bz = bz(x, y);

                localMagneticBasis({Bx, By, Bz}, basis);
            }


            for (std::uint32_t ipart = 0; ipart < nbrParticlePerCell_; ++ipart)
            {
                maxwellianVelocity({Vx, Vy, Vz}, {Vthx, Vthy, Vthz}, generator, particleVelocity);

                if (basis_ == Basis::Magnetic)
                {
                    particleVelocity = basisTransform(basis, particleVelocity);
                }

                std::array<float, dimension> delta = {{randPosX(generator), randPosY(generator)}};

                Particle<dimension> tmpParticle;

                // particle iCell is in AMR index
                auto AMRCellIndex = layout.localToAMR(Point{ix, iy});


                tmpParticle.weight = cellWeight;
                tmpParticle.charge = particleCharge_;
                tmpParticle.iCell  = AMRCellIndex.template toArray<int>();
                tmpParticle.delta  = delta;
                tmpParticle.v      = particleVelocity;

                particles.push_back(std::move(tmpParticle));
            }
        }
    }
}

} // namespace PHARE::core::test_fn::func_2d

namespace PHARE::core::test_fn::func_3d
{
template<typename ParticleArray, typename GridLayout>
class TestMaxwellianParticleInitializer
{
    constexpr static std::size_t dimension           = 3;
    constexpr static std::size_t nbrParticlePerCell_ = 10;

    using MaxInit = MaxwellianParticleInitializer<ParticleArray, GridLayout>;

public:
    TestMaxwellianParticleInitializer(ParticleArray& particles, GridLayout const& layout);

private:
    Basis basis_           = Basis::Cartesian;
    double particleCharge_ = 1.1;
    std::optional<std::size_t> rngSeed_{1337};
};

template<typename ParticleArray, typename GridLayout>
TestMaxwellianParticleInitializer<ParticleArray, GridLayout>::TestMaxwellianParticleInitializer(
    ParticleArray& particles, GridLayout const& layout)
{
    auto const meshSize = layout.meshSize();
    auto const& dx      = meshSize[0];
    auto const& dy      = meshSize[1];
    auto const& dz      = meshSize[2];

    /* get indices start and stop. we take primal/primal/primal because
       that is what GridLayout::cellCenteredCoordinate() requires */
    auto const [ix0, ix1] = layout.physicalStartToEnd(QtyCentering::primal, Direction::X);
    auto const [iy0, iy1] = layout.physicalStartToEnd(QtyCentering::primal, Direction::Y);
    auto const [iz0, iz1] = layout.physicalStartToEnd(QtyCentering::primal, Direction::Z);

    double cellVolume = dx * dy * dz;

    // random seed and generator needed to load maxwellian velocity
    // and random position with the cell
    auto generator = MaxInit::getRNG(rngSeed_);

    // beware: we're looping over the cell but use primal indices because of
    // GridLayout::cellCenteredCoordinates
    // therefore i(x,y,z)1 must be excluded.

    // gab references for convenience
    // auto& density         = *density_;
    // auto& bulkVelocity    = *bulkVelocity_;
    // auto& thermalVelocity = *thermalVelocity_;


    for (std::uint32_t ix = ix0; ix < ix1; ++ix)
    {
        for (std::uint32_t iy = iy0; iy < iy1; ++iy)
        {
            for (std::uint32_t iz = iz0; iz < iz1; ++iz)
            {
                double n; // cell centered density
                std::array<double, 3> particleVelocity;
                std::array<std::array<double, 3>, 3> basis;

                // get the coordinate of the current cell
                auto coord = layout.cellCenteredCoordinates(ix, iy, iz);
                auto x     = coord[0];
                auto y     = coord[1];
                auto z     = coord[2];

                // now get density, velocity and thermal speed values
                n         = density_(x, y, z);
                auto Vx   = vx(x, y, z);
                auto Vy   = vy(x, y, z);
                auto Vz   = vz(x, y, z);
                auto Vthx = vthx(x, y, z);
                auto Vthy = vthy(x, y, z);
                auto Vthz = vthz(x, y, z);

                // weight for all particles in this cell
                auto cellWeight = n * cellVolume / nbrParticlePerCell_;

                ParticleDeltaDistribution randPosX;
                ParticleDeltaDistribution randPosY;
                ParticleDeltaDistribution randPosZ;

                if (basis_ == Basis::Magnetic)
                {
                    auto Bx = bx(x, y, z);
                    auto By = by(x, y, z);
                    auto Bz = bz(x, y, z);

                    localMagneticBasis({Bx, By, Bz}, basis);
                }

                for (std::uint32_t ipart = 0; ipart < nbrParticlePerCell_; ++ipart)
                {
                    maxwellianVelocity({Vx, Vy, Vz}, {Vthx, Vthy, Vthz}, generator,
                                       particleVelocity);

                    if (basis_ == Basis::Magnetic)
                    {
                        particleVelocity = basisTransform(basis, particleVelocity);
                    }

                    std::array<float, dimension> delta
                        = {{randPosX(generator), randPosY(generator), randPosZ(generator)}};


                    Particle<dimension> tmpParticle;

                    // particle iCell is in AMR index
                    auto AMRCellIndex = layout.localToAMR(Point{ix, iy, iz});

                    tmpParticle.weight = cellWeight;
                    tmpParticle.charge = particleCharge_;
                    tmpParticle.iCell  = AMRCellIndex.template toArray<int>();
                    tmpParticle.delta  = delta;
                    tmpParticle.v      = particleVelocity;

                    particles.push_back(std::move(tmpParticle));
                } // end particle looop
            } // end z
        } // end y
    } // end x
}



} // namespace PHARE::core::test_fn::func_3d


#endif // PHARE_TEST_MAXWELLIAN_PARTICLE_INITIALIZER_TEST_INIT_FUNCS_HPP
