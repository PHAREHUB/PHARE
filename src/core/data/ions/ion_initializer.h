#ifndef PHARE_IONS_INITIALIZER_H
#define PHARE_IONS_INITIALIZER_H


#include <memory>
#include <string>
#include <vector>

//#include "particleinitializer.h"
#include "utilities/types.h"

namespace PHARE
{
template<typename ParticleArray>
class ParticleInitializer
{
protected:
public:
    virtual void loadParticles(ParticleArray& particles) const = 0;

    virtual ~ParticleInitializer() = 0;
};




/**
 * @brief IonInitializer is need for Ions construction so that Ions
 *       may be initialized
 *
 *  IonInitializer objects are created by InitializerFactory objects
 *  which build and initialize them either from user input parameters when
 *  the simulation is being initialized or at regridding step.
 *
 * It contains all necessary information of Ions to build their populations. In
 * particular, it gives access for each population to
 *
 *  - its ParticleInitializer.
 *  - the mass of its particles
 *  - its name
 *
 */
template<typename ParticleArray>
struct IonsInitializer
{
    std::string name;
    std::vector<std::unique_ptr<ParticleInitializer<ParticleArray>>> particleInitializers;
    std::vector<double> masses;
    std::vector<std::string> names;
    uint32 nbrPopulations;
};
} // namespace PHARE


#endif // PHARE_IONS_INITIALIZER_H
