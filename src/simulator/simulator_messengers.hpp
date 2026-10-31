#ifndef PHARE_SIMULATOR_SIMULATOR_MESSENGERS_HPP
#define PHARE_SIMULATOR_SIMULATOR_MESSENGERS_HPP

#include "simulator/simulator_decl.hpp"

namespace PHARE
{

template<auto opts>
void Simulator<opts>::messengers_init()
{
    multiphysInteg_->registerAndSetupMessengers(messengerFactory_);
}

} // namespace PHARE

#endif /*PHARE_SIMULATOR_SIMULATOR_MESSENGERS_HPP*/
