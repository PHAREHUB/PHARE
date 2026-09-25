#ifndef PHARE_CORE_UTILITIES_SPACE_TIME_FUNCTION_HPP
#define PHARE_CORE_UTILITIES_SPACE_TIME_FUNCTION_HPP

#include "core/utilities/span.hpp"

#include <cstddef>
#include <functional>
#include <memory>

namespace PHARE::core
{

/**
 * @brief A user function of space and time, f(x[, y[, z]], t).
 */
template<typename ReturnType, std::size_t dim>
struct SpaceTimeFunctionHelper
{
};

template<>
struct SpaceTimeFunctionHelper<double, 1>
{
    using return_type = std::shared_ptr<Span<double>>;
    using param_type  = core::Span<double const> const&;
    using type        = std::function<return_type(param_type, double)>;
};

template<>
struct SpaceTimeFunctionHelper<double, 2>
{
    using return_type = std::shared_ptr<Span<double>>;
    using param_type  = core::Span<double const> const&;
    using type        = std::function<return_type(param_type, param_type, double)>;
};

template<>
struct SpaceTimeFunctionHelper<double, 3>
{
    using return_type = std::shared_ptr<Span<double>>;
    using param_type  = core::Span<double const> const&;
    using type        = std::function<return_type(param_type, param_type, param_type, double)>;
};

template<std::size_t dim>
using SpaceTimeFunction = typename SpaceTimeFunctionHelper<double, dim>::type;

} // namespace PHARE::core

#endif // PHARE_CORE_UTILITIES_SPACE_TIME_FUNCTION_HPP
