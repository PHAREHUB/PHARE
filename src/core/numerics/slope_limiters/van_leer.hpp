#ifndef CORE_NUMERICS_SLOPE_LIMITER_VAN_LEER_HPP
#define CORE_NUMERICS_SLOPE_LIMITER_VAN_LEER_HPP

namespace PHARE::core
{
struct VanLeerLimiter
{
    static auto limit(auto const& Dil, auto const& Dir)
    {
        return Dil * Dir > 0.0 ? 2.0 * Dil * Dir / (Dil + Dir) : 0.0;
    }
};
} // namespace PHARE::core

#endif
