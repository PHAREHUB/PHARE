#ifndef PHARE_ION_UPDATER_HPP
#define PHARE_ION_UPDATER_HPP


#include "core/numerics/ion_updater/ion_updater/ion_updater_impl0.hpp"


namespace PHARE::core
{

template<typename... Args>
using IonUpdater = IonUpdater0<Args...>;


} // namespace PHARE::core


#endif // ION_UPDATER_HPP
