#ifndef PHARE_FLUID_PARTICLE_INITIALIZER_H
#define PHARE_FLUID_PARTICLE_INITIALIZER_H

#include <functional>
#include <memory>
#include <random>

#include "core/data/grid/gridlayoutdefs.h"
#include "core/hybrid/hybrid_quantities.h"
#include "core/utilities/types.h"
#include "core/data/ions/particle_initializers/particle_initializer.h"
#include "core/data/particles/particle.h"
#include "initializer/data_provider.h"
#include "core/utilities/point/point.h"

namespace PHARE
{
namespace core
{
    void maxwellianVelocity(std::array<double, 3> V, std::array<double, 3> Vth,
                            std::mt19937_64& generator, std::array<double, 3>& partVelocity);


    std::array<double, 3> basisTransform(std::array<std::array<double, 3>, 3> const& basis,
                                         std::array<double, 3> const& vec);

    void localMagneticBasis(std::array<double, 3> B, std::array<std::array<double, 3>, 3>& basis);



    /** @brief a MaxwellianParticleInitializer is a ParticleInitializer that loads particles from a
     * local Maxwellian distribution given density, bulk velocity and thermal velocity profiles.
     */
    template<typename ParticleArray, typename GridLayout>
    class MaxwellianParticleInitializer : public ParticleInitializer<ParticleArray, GridLayout>
    {
    public:
        static constexpr auto dimension = GridLayout::dimension;

        MaxwellianParticleInitializer(
            PHARE::initializer::ScalarFunction<dimension> density,
            std::array<PHARE::initializer::ScalarFunction<dimension>, 3> bulkVelocity,
            std::array<PHARE::initializer::ScalarFunction<dimension>, 3> thermalVelocity,
            double particleCharge, std::uint32_t nbrParticlesPerCell,
            std::optional<std::size_t> seed = {}, Basis basis = Basis::Cartesian,
            std::array<PHARE::initializer::ScalarFunction<dimension>, 3> magneticField
            = {nullptr, nullptr, nullptr})
            : density_{density}
            , bulkVelocity_{bulkVelocity}
            , thermalVelocity_{thermalVelocity}
            , magneticField_{magneticField}
            , particleCharge_{particleCharge}
            , nbrParticlePerCell_{nbrParticlesPerCell}
            , basis_{basis}
            , rngSeed_{seed}
        {
        }


        /**
         * @brief load particles in a ParticleArray in a domain defined by the given layout
         */
        virtual void loadParticles(ParticleArray& particles,
                                   GridLayout const& layout) const override
        {
            if constexpr (dimension == 1)
            {
                loadParticles1D_(particles, layout);
            }
            else if constexpr (dimension == 2)
            {
                loadParticles2D_(particles, layout);
            }
            else if constexpr (dimension == 3)
            {
                loadParticles3D_(particles, layout);
            }
        }


        virtual ~MaxwellianParticleInitializer() = default;



    private:
        // can be relocated if needed
        static std::mt19937_64 getRNG(std::optional<std::size_t> const& seed)
        {
            if (!seed.has_value())
            {
                std::random_device randSeed;
                std::seed_seq seed_seq{randSeed(), randSeed(), randSeed(), randSeed(),
                                       randSeed(), randSeed(), randSeed(), randSeed()};
                return std::mt19937_64(seed_seq);
            }
            return std::mt19937_64(*seed);
        }



        void loadParticles1D_(ParticleArray& particles, GridLayout const& layout) const
        {
            auto const meshSize = layout.meshSize();
            auto const& dx      = meshSize[0];

            /* get indices start and stop. we take primal/primal/primal because
               that is what GridLayout::cellCenteredCoordinate() requires */
            std::uint32_t ix0 = layout.physicalStartIndex(QtyCentering::primal, Direction::X);
            std::uint32_t ix1 = layout.physicalEndIndex(QtyCentering::primal, Direction::X);

            double cellVolume = dx;

            // random seed and generator needed to load maxwellian velocity
            // and random position with the cell
            auto generator = getRNG(rngSeed_);

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
                auto Vx = bulkVelocity_[0](x);
                auto Vy = bulkVelocity_[1](x);
                auto Vz = bulkVelocity_[2](x);

                auto Vthx = thermalVelocity_[0](x);
                auto Vthy = thermalVelocity_[1](x);
                auto Vthz = thermalVelocity_[2](x);

                // weight for all particles in this cell
                auto cellWeight = n * cellVolume / nbrParticlePerCell_;

                ParticleDeltaDistribution randPosX;

                if (basis_ == Basis::Magnetic)
                {
                    auto Bx = magneticField_[0](x);
                    auto By = magneticField_[1](x);
                    auto Bz = magneticField_[2](x);

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




        void loadParticles2D_(ParticleArray& particles, GridLayout const& layout) const
        {
            auto const meshSize = layout.meshSize();
            auto const& dx      = meshSize[0];
            auto const& dy      = meshSize[1];

            /* get indices start and stop. we take primal/primal/primal because
               that is what GridLayout::cellCenteredCoordinate() requires */
            auto ix0 = layout.physicalStartIndex(QtyCentering::primal, Direction::X);
            auto ix1 = layout.physicalEndIndex(QtyCentering::primal, Direction::X);
            auto iy0 = layout.physicalStartIndex(QtyCentering::primal, Direction::Y);
            auto iy1 = layout.physicalEndIndex(QtyCentering::primal, Direction::Y);

            auto cellVolume = dx * dy;

            // random seed and generator needed to load maxwellian velocity
            // and random position with the cell
            auto generator = getRNG(rngSeed_);

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
                    auto Vx   = bulkVelocity_[0](x, y);
                    auto Vy   = bulkVelocity_[1](x, y);
                    auto Vz   = bulkVelocity_[2](x, y);
                    auto Vthx = thermalVelocity_[0](x, y);
                    auto Vthy = thermalVelocity_[1](x, y);
                    auto Vthz = thermalVelocity_[2](x, y);

                    // weight for all particles in this cell
                    auto cellWeight = n * cellVolume / nbrParticlePerCell_;

                    ParticleDeltaDistribution randPosX;
                    ParticleDeltaDistribution randPosY;

                    if (basis_ == Basis::Magnetic)
                    {
                        auto Bx = magneticField_[0](x, y);
                        auto By = magneticField_[1](x, y);
                        auto Bz = magneticField_[2](x, y);

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
                            = {{randPosX(generator), randPosY(generator)}};

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




        void loadParticles3D_(ParticleArray& particles, GridLayout const& layout) const
        {
            auto const meshSize = layout.meshSize();
            auto const& dx      = meshSize[0];
            auto const& dy      = meshSize[1];
            auto const& dz      = meshSize[2];

            /* get indices start and stop. we take primal/primal/primal because
               that is what GridLayout::cellCenteredCoordinate() requires */
            auto ix0 = layout.physicalStartIndex(QtyCentering::primal, Direction::X);
            auto ix1 = layout.physicalEndIndex(QtyCentering::primal, Direction::X);
            auto iy0 = layout.physicalStartIndex(QtyCentering::primal, Direction::Y);
            auto iy1 = layout.physicalEndIndex(QtyCentering::primal, Direction::Y);
            auto iz0 = layout.physicalStartIndex(QtyCentering::primal, Direction::Z);
            auto iz1 = layout.physicalEndIndex(QtyCentering::primal, Direction::Z);

            double cellVolume = dx * dy * dz;

            // random seed and generator needed to load maxwellian velocity
            // and random position with the cell
            auto generator = getRNG(rngSeed_);

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
                        auto Vx   = bulkVelocity_[0](x, y, z);
                        auto Vy   = bulkVelocity_[1](x, y, z);
                        auto Vz   = bulkVelocity_[2](x, y, z);
                        auto Vthx = thermalVelocity_[0](x, y, z);
                        auto Vthy = thermalVelocity_[1](x, y, z);
                        auto Vthz = thermalVelocity_[2](x, y, z);

                        // weight for all particles in this cell
                        auto cellWeight = n * cellVolume / nbrParticlePerCell_;

                        ParticleDeltaDistribution randPosX;
                        ParticleDeltaDistribution randPosY;
                        ParticleDeltaDistribution randPosZ;

                        if (basis_ == Basis::Magnetic)
                        {
                            auto Bx = magneticField_[0](x, y, z);
                            auto By = magneticField_[0](x, y, z);
                            auto Bz = magneticField_[0](x, y, z);

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

                            std::array<float, 3> delta
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
                    }     // end z
                }         // end y
            }             // end x
        }



        PHARE::initializer::ScalarFunction<dimension> density_;
        std::array<PHARE::initializer::ScalarFunction<dimension>, 3> bulkVelocity_;
        std::array<PHARE::initializer::ScalarFunction<dimension>, 3> thermalVelocity_;
        std::array<PHARE::initializer::ScalarFunction<dimension>, 3> magneticField_;

        double particleCharge_;
        std::uint32_t nbrParticlePerCell_;
        Basis basis_;
        std::optional<std::size_t> rngSeed_;
    };
} // namespace core
} // namespace PHARE


#endif
