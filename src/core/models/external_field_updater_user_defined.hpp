#ifndef PHARE_CORE_MODELS_EXTERNAL_FIELD_UPDATER_USER_DEFINED_HPP
#define PHARE_CORE_MODELS_EXTERNAL_FIELD_UPDATER_USER_DEFINED_HPP

#include "core/models/external_field_updater.hpp"
#include "core/data/user/user_field_updater.hpp"
#include "initializer/data_provider.hpp"

#include <array>
#include <cstddef>
#include <optional>
#include <stdexcept>
#include <utility>

namespace PHARE::core
{


/**
 * @brief External field updater with the vector potential and its derivative with respect to time
 * specified by python user-defined functions.
 *
 * @tparam VecFieldT vecfield implementation
 * @tparam GridLayoutT grid layout implementation
 *
 */
template<typename VecFieldT, typename GridLayoutT>
class ExternalFieldUpdaterUserDefined : public IExternalFieldUpdater<VecFieldT, GridLayoutT>
{
public:
    using Super         = IExternalFieldUpdater<VecFieldT, GridLayoutT>;
    using vecfield_type = Super::vecfield_type;

    static constexpr std::size_t dimension = Super::dimension;
    static constexpr std::size_t N         = vecfield_type::size();

    using space_time_function_type       = initializer::SpaceTimeFunction<dimension>;
    using space_time_function_array_type = std::array<space_time_function_type, N>;

    // constructor parameters have rvalue reference type to force calling the constructor with
    // arguments wrapped in std::move, thus avoiding any copy
    ExternalFieldUpdaterUserDefined(space_time_function_array_type&& potential,
                                    std::optional<space_time_function_array_type>&& derivative
                                    = std::nullopt)
        : potential_{std::move(potential)}
        , potential_time_derivative_{std::move(derivative)}
    {
    }

    NO_DISCARD bool isTimeDependent() const final { return potential_time_derivative_.has_value(); }

    void computePotential(vecfield_type& a0, double time, GridLayoutT const& layout) final
    {
        UserFieldUpdater::update(a0, layout, potential_, time);
    }

    void computePotentialTimeDerivative(vecfield_type& da0_dt, double time,
                                        GridLayoutT const& layout) final
    {
        if (potential_time_derivative_)
            UserFieldUpdater::update(da0_dt, layout, potential_time_derivative_.value(), time);
        else
            throw std::runtime_error(
                "computePotentialTimeDerivative called on a constant external field.");
    }

private:
    space_time_function_array_type potential_;
    std::optional<space_time_function_array_type> potential_time_derivative_;
};

} // namespace PHARE::core

#endif // PHARE_CORE_MODELS_EXTERNAL_FIELD_UPDATER_USER_DEFINED_HPP
