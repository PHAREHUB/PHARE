#ifndef PHARE_ION_UPDATER_HPP
#define PHARE_ION_UPDATER_HPP


#include "core/numerics/ion_updater/ion_updater/ion_updater_impl0.hpp"
#include "core/numerics/ion_updater/ion_updater/ion_updater_impl1.hpp"

#ifndef PHARE_UPDATER_IMPL
#define PHARE_UPDATER_IMPL 1
#endif


namespace PHARE::core
{

template<typename... Args>
using IonUpdater
    = std::conditional_t<PHARE_UPDATER_IMPL == 0, IonUpdater0<Args...>, IonUpdater1<Args...>>;


} // namespace PHARE::core


#endif // ION_UPDATER_HPP
