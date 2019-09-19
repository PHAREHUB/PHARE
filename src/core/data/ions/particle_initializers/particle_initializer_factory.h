#ifndef PHARE_PARTICLE_INITIALIZER_FACTORY_H
#define PHARE_PARTICLE_INITIALIZER_FACTORY_H


#include "data_provider.h"
#include "maxwellian_particle_initializer.h"
#include "particle_initializer.h"
#include "utilities/types.h"

namespace PHARE
{
namespace core
{
    template<typename ParticleArray, typename GridLayout>
    class ParticleInitializerFactory
    {
        using ParticleInitializerT      = ParticleInitializer<ParticleArray, GridLayout>;
        static constexpr auto dimension = GridLayout::dimension;

    public:
        static std::unique_ptr<ParticleInitializerT> create(PHARE::initializer::PHAREDict<1>& dict)
        {
            auto initializerName = dict["name"].to<std::string>();

            std::string toto;

            if (initializerName == "MaxwellianParticleInitializer")
            {
                auto& density = dict["density"].to<PHARE::initializer::ScalarFunction<dimension>>();
                auto& bulkVel
                    = dict["bulkVelocity"].to<PHARE::initializer::VectorFunction<dimension>>();

                auto& thermalVel
                    = dict["thermalVelocity"].to<PHARE::initializer::VectorFunction<dimension>>();

                auto charge = dict["charge"].to<double>();

                auto nbrPartPerCell = dict["nbrPartPerCell"].to<std::size_t>();

                auto basisName = dict["basis"].to<std::string>();


                if (basisName == "Cartesian")
                {
                    [[maybe_unused]] Basis basis = Basis::Cartesian;
                    return std::make_unique<
                        MaxwellianParticleInitializer<ParticleArray, GridLayout>>(
                        density, bulkVel, thermalVel, charge, nbrPartPerCell);
                }
                else if (basisName == "Magnetic")
                {
                    [[maybe_unused]] Basis basis = Basis::Magnetic;
                    [[maybe_unused]] auto& magnetic
                        = dict["magnetic"].to<PHARE::initializer::VectorFunction<dimension>>();
                    return std::make_unique<
                        MaxwellianParticleInitializer<ParticleArray, GridLayout>>(
                        density, bulkVel, thermalVel, charge, nbrPartPerCell);
                }
            }
            // TODO throw?
            return nullptr;
        }
    };

} // namespace core

} // namespace PHARE


#endif
