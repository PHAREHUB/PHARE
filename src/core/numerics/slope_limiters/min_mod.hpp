#ifndef CORE_NUMERICS_SLOPE_LIMITER_MIN_MOD_HPP
#define CORE_NUMERICS_SLOPE_LIMITER_MIN_MOD_HPP

namespace PHARE::core
{
struct MinModLimiter
{
    static auto limit(auto const& Dil, auto const& Dir)
    {
        return Dil * Dir < 0.0 ? 0.0 : fabs(Dir) < fabs(Dil) ? Dir : Dil;
    }
};
} // namespace PHARE::core

#endif
