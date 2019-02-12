#ifndef PHARE_PARTICLE_INITIALIZER_H
#define PHARE_PARTICLE_INITIALIZER_H


namespace PHARE
{
namespace core
{
    template<typename ParticleArray, typename GridLayout>
    class ParticleInitializer
    {
    public:
        virtual void loadParticles(ParticleArray& particles, GridLayout const& layout) const = 0;
        virtual ~ParticleInitializer() = default;
    };

} // namespace core
} // namespace PHARE

#endif
