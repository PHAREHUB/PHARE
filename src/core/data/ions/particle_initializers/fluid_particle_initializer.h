#ifndef PHARE_FLUID_PARTICLE_INITIALIZER_H
#define PHARE_FLUID_PARTICLE_INITIALIZER_H

#include <functional>
#include <memory>

#include "data/ions/particle_initializers/particle_initializer.h"

namespace PHARE
{
template<typename ParticleArray, typename GridLayout>
class FluidParticleInitializer : public ParticleInitializer<ParticleArray, GridLayout>
{
public:
    virtual void loadParticles(ParticleArray& particles, GridLayout const& laout) const override {}

private:
    /* template<typename... Coords>
     std::function<double(Coords...)> density_;*/
    /*
    std::unique_ptr<ScalarFunction> density_;
    std::unique_ptr<VectorFunction> bulkVelocity_;
    std::unique_ptr<VectorFunction> thermalSpeed_;
    std::unique_ptr<VectorFunction> magneticField_;
    double particleCharge_;
    uint32 nbrParticlePerCell_;
    Basis base_;


    void loadParticles1D_(std::vector<Particle>& particles) const;
    void loadParticles2D_(std::vector<Particle>& particles) const;
    void loadParticles3D_(std::vector<Particle>& particles) const;

public:
    FluidParticleInitializer(GridLayout const& layout,
                             std::unique_ptr<ScalarFunction> densityProfile,
                             std::unique_ptr<VectorFunction> bulkVelocityProfile,
                             std::unique_ptr<VectorFunction> thermalSpeedProfile,
                             uint32 nbrPartPerCell, double particleCharge,
                             Basis base                                    = Basis::Cartesian,
                             std::unique_ptr<VectorFunction> magneticField = nullptr);

    virtual void loadParticles(std::vector<Particle>& particles) const override;
    */
};

} // namespace PHARE


#endif
