#ifndef PHARE_CORE_MODELS_EXTERNAL_FIELD_UPDATER_FACTORY_HPP
#define PHARE_CORE_MODELS_EXTERNAL_FIELD_UPDATER_FACTORY_HPP

#include "core/models/external_field_updater.hpp"
#include "core/models/external_field_updater_defs.hpp"
#include "core/models/external_field_updater_dipole.hpp"
#include "core/models/external_field_updater_user_defined.hpp"
#include "core/models/external_field_updater_zero.hpp"

#include "initializer/data_provider.hpp"
#include "initializer/dict_utils.hpp"

#include <memory>
#include <optional>
#include <stdexcept>

namespace PHARE::core
{
/**
 * @brief Factory for the external field updater.
 */
template<typename VecFieldT, typename GridLayoutT>
class ExternalFieldUpdaterFactory
{
private:
    using Interface   = IExternalFieldUpdater<VecFieldT, GridLayoutT>;
    using Dipole      = ExternalFieldUpdaterDipole<VecFieldT, GridLayoutT>;
    using Zero        = ExternalFieldUpdaterZero<VecFieldT, GridLayoutT>;
    using UserDefined = ExternalFieldUpdaterUserDefined<VecFieldT, GridLayoutT>;

public:
    using point_type               = Interface::point_type;
    using value_type               = Interface::value_type;
    using space_time_function_type = initializer::SpaceTimeFunction<GridLayoutT::dimension>;

    static constexpr std::size_t dimension = GridLayoutT::dimension;
    static constexpr std::size_t N         = VecFieldT::size();

    ExternalFieldUpdaterFactory() = delete;

    static std::unique_ptr<Interface> createZero() { return std::make_unique<Zero>(); }

    static std::unique_ptr<Interface> createDipole(initializer::PHAREDict const& dict)
    {
        auto position
            = point_type{initializer::parseDimXYZType<double, dimension>(dict, "position")};
        auto moment = typename Dipole::vector_type{
            initializer::parseDimXYZType<value_type, dimension>(dict, "moment")};
        if (!dict.contains("radius")) // no default: 0 (a point dipole) must be explicit
            throw std::runtime_error("dipole external field requires a 'radius'");
        auto const radius = dict["radius"].template to<double>();
        return std::make_unique<Dipole>(position, moment, radius);
    }

    static std::unique_ptr<Interface> createUserDefined(initializer::PHAREDict const& dict)
    {
        auto potential
            = initializer::parseDimXYZType<space_time_function_type, N>(dict, "potential");

        std::optional<std::array<space_time_function_type, N>> derivative;
        if (dict["is_time_dependent"].template to<bool>())
            derivative = initializer::parseDimXYZType<space_time_function_type, N>(
                dict, "potential_time_derivative");

        return std::make_unique<UserDefined>(std::move(potential), std::move(derivative));
    }

    static std::unique_ptr<Interface> create(initializer::PHAREDict const& dict)
    {
        auto const type = cppdict::get_value(dict, "type", ExternalFieldUpdaterType::Zero);

        switch (type)
        {
            case ExternalFieldUpdaterType::Zero: return createZero();
            case ExternalFieldUpdaterType::Dipole: return createDipole(dict);
            case ExternalFieldUpdaterType::UserDefined: return createUserDefined(dict);
        }
        throw std::runtime_error("external field updater: unknown type");
    }
};

} // namespace PHARE::core


#endif // PHARE_CORE_MODELS_EXTERNAL_FIELD_UPDATER_FACTORY_HPP
