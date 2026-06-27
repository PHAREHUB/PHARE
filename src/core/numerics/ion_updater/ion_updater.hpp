#ifndef PHARE_ION_UPDATER_HPP
#define PHARE_ION_UPDATER_HPP

#ifndef PHARE_UPDATER_IMPL
#define PHARE_UPDATER_IMPL 1 // 1 == new impl with less allocation
#endif

#if PHARE_UPDATER_IMPL == 0
#include "core/numerics/ion_updater/ion_updater/ion_updater_impl0.hpp"
#elif PHARE_UPDATER_IMPL == 1
#include "core/numerics/ion_updater/ion_updater/ion_updater_impl1.hpp"
#else
#error
#endif

namespace PHARE::core
{

#if PHARE_UPDATER_IMPL == 0
template<typename... Args>
using IonUpdater = IonUpdater0<Args...>;

#elif PHARE_UPDATER_IMPL == 1
template<typename... Args>
using IonUpdater = IonUpdater1<Args...>;

#endif

} // namespace PHARE::core

#endif // ION_UPDATER_HPP
