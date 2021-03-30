#ifndef PHARE_PARTICLE_INITIALIZER_H
#define PHARE_PARTICLE_INITIALIZER_H


#include <optional>


namespace PHARE::core
{
struct ParticleInitiazationInfo
{
    ParticleInitiazationInfo(std::string const& _seed_mode = "none")
        : seed_mode{_seed_mode}
    {
    }

    std::string seed_mode;
    std::optional<std::size_t> seed;
};

template<typename ParticleArray, typename GridLayout>
class ParticleInitializer
{
public:
    virtual void loadParticles(ParticleArray& particles, GridLayout const& layout) const = 0;
    virtual ~ParticleInitializer()                                                       = default;
};

} // namespace PHARE::core

#endif
