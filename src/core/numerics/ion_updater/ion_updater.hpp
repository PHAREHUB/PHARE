#ifndef PHARE_ION_UPDATER_HPP
#define PHARE_ION_UPDATER_HPP


#include "core/numerics/ion_updater/ion_updater/ion_updater_impl0.hpp"
#include "core/numerics/ion_updater/ion_updater/ion_updater_impl1.hpp"
#include "core/numerics/ion_updater/ion_updater/ion_updater_impl2.hpp"

#ifndef PHARE_UPDATER_IMPL
#define PHARE_UPDATER_IMPL 0
#endif


namespace PHARE::core
{

template<typename... Args>
using IonUpdater = std::conditional_t<
    PHARE_UPDATER_IMPL == 0, IonUpdater0<Args...>,
    std::conditional_t<PHARE_UPDATER_IMPL == 1, IonUpdater1<Args...>, IonUpdater2<Args...>>>;


} // namespace PHARE::core


#endif // ION_UPDATER_HPP
