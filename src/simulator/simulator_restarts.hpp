#ifndef PHARE_SIMULATOR_SIMULATOR_RESTARTS_HPP
#define PHARE_SIMULATOR_SIMULATOR_RESTARTS_HPP

#include "simulator/simulator_decl.hpp"

#include "restarts/restarts.hpp"

namespace PHARE
{

template<auto opts>
void Simulator<opts>::hybrid_restarts_init(initializer::PHAREDict const& dict)
{
    if constexpr (has_hybrid_v<opts>)
        rMan = restarts::RestartsManagerResolver::make_unique(*hierarchy_, *hyb_.resman_, dict);
}


template<auto opts>
void Simulator<opts>::mhd_restarts_init(initializer::PHAREDict const& dict)
{
    if constexpr (has_mhd_v<opts>)
        rMan = restarts::RestartsManagerResolver::make_unique(*hierarchy_, *mhd_.resman_, dict);
}

} // namespace PHARE

#endif /*PHARE_SIMULATOR_SIMULATOR_RESTARTS_HPP*/
