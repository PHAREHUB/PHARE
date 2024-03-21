#ifndef PHARE_PARTICLE_INITIALIZER_FACTORY_HPP
#define PHARE_PARTICLE_INITIALIZER_FACTORY_HPP


#include "core/def.hpp"
#include "core/utilities/types.hpp"
#include "initializer/data_provider.hpp"
#include "maxwellian_particle_initializer.hpp"
#include "particle_initializer.hpp"

#include <memory>

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
        NO_DISCARD static std::unique_ptr<ParticleInitializerT>
        create(initializer::PHAREDict const& dict)
        {
            using FunctionType = initializer::InitFunction<dimension>;

            auto initializerName = dict["name"].template to<std::string>();

            if (initializerName == "maxwellian")
            {
                auto& density = dict["density"].template to<FunctionType>();

                auto& bulkVelx = dict["bulk_velocity_x"].template to<FunctionType>();

                auto& bulkVely = dict["bulk_velocity_y"].template to<FunctionType>();

                auto& bulkVelz = dict["bulk_velocity_z"].template to<FunctionType>();

                auto& vthx = dict["thermal_velocity_x"].template to<FunctionType>();

                auto& vthy = dict["thermal_velocity_y"].template to<FunctionType>();

                auto& vthz = dict["thermal_velocity_z"].template to<FunctionType>();

                auto charge = dict["charge"].template to<double>();

                auto nbrPartPerCell = cppdict::get_value(dict, "nbr_part_per_cell", int{0});
                FunctionType nbrPartPerCellFn
                    = cppdict::get_value(dict, "nbr_part_per_cell_fn", FunctionType{nullptr});
                if (not nbrPartPerCellFn and nbrPartPerCell == 0)
                {
                    throw std::runtime_error("PPC cannot be 0");
                }

                auto basisName = dict["basis"].template to<std::string>();

                std::array<FunctionType, 3> v = {bulkVelx, bulkVely, bulkVelz};

                std::array<FunctionType, 3> vth = {vthx, vthy, vthz};

                std::optional<std::size_t> seed;
                if (dict.contains("init") && dict["init"].contains("seed"))
                    seed = dict["init"]["seed"].template to<std::optional<std::size_t>>();

                std::array<FunctionType, 3> magneticField = {nullptr, nullptr, nullptr};

                if (basisName == "cartesian")
                {
                    return std::make_unique<
                        MaxwellianParticleInitializer<ParticleArray, GridLayout>>(
                        density, v, vth, charge, nbrPartPerCell, seed, Basis::Cartesian,
                        magneticField, nbrPartPerCellFn);
                }
                else if (basisName == "magnetic")
                {
                    magneticField[0] = dict["magnetic_x"].template to<FunctionType>();
                    magneticField[1] = dict["magnetic_x"].template to<FunctionType>();
                    magneticField[2] = dict["magnetic_x"].template to<FunctionType>();

                    return std::make_unique<
                        MaxwellianParticleInitializer<ParticleArray, GridLayout>>(
                        density, v, vth, charge, nbrPartPerCell, seed, Basis::Magnetic,
                        magneticField, nbrPartPerCellFn);
                }
            }
            // TODO throw?
            return nullptr;
        }
    };

} // namespace core

} // namespace PHARE

#endif
