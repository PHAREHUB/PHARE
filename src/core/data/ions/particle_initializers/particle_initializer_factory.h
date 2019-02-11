#ifndef PHARE_PARTICLE_INITIALIZER_FACTORY_H
#define PHARE_PARTICLE_INITIALIZER_FACTORY_H


#include "data_provider.h"
#include "fluid_particle_initializer.h"
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

            if (initializerName == "FluidParticleInitializer")
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
                    Basis basis = Basis::Cartesian;
                    return std::make_unique<FluidParticleInitializer<ParticleArray, GridLayout>>(
                        density, bulkVel, thermalVel, charge, nbrPartPerCell);
                }
                else if (basisName == "Magnetic")
                {
                    Basis basis = Basis::Magnetic;
                    auto& magnetic
                        = dict["magnetic"].to<PHARE::initializer::VectorFunction<dimension>>();
                    return std::make_unique<FluidParticleInitializer<ParticleArray, GridLayout>>(
                        density, bulkVel, thermalVel, charge, nbrPartPerCell);
                }
                /*

                auto& magnetic
                    =
 dict["magnetic"].to<PHARE::initializer::VectorFunction<dimension>>();*/
            }

            return nullptr;
        }
    };

} // namespace core

} // namespace PHARE


#endif
