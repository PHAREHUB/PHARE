#ifndef PHARE_ION_UPDATER_HPP
#define PHARE_ION_UPDATER_HPP


#include "core/logger.hpp"

#include "core/numerics/ion_updater/ion_updater/ion_updater_impl0.hpp"
#include "core/numerics/ion_updater/ion_updater/ion_updater_impl1.hpp"
#include "core/numerics/ion_updater/ion_updater/ion_updater_impl2.hpp"

#include <stdexcept>

#ifndef PHARE_UPDATER_IMPL
#define PHARE_UPDATER_IMPL 0
#endif


namespace PHARE::core::detail
{

template<typename Impl_t>
auto reset_updater_impl(Impl_t& impl) -> decltype(impl.reset(), void())
{
    impl.reset();
}

auto reset_updater_impl(auto&&...) {}

} // namespace PHARE::core::detail

namespace PHARE::core
{


template<typename Impl_t>
class IonUpdaterProxy
{
public:
    using Interpolator_t = Impl_t::Interpolator_t;

    IonUpdaterProxy(auto&&... args)
        : impl{args...}
    {
    }

    void updatePopulations(auto&&... args);

    void updateIons(auto& ions);

    void reset() { detail::reset_updater_impl(impl); }

private:
    Impl_t impl;
};


template<typename Impl_t>
void IonUpdaterProxy<Impl_t>::updatePopulations(auto&&... args)
{
    PHARE_LOG_SCOPE(3, "IonUpdater::updatePopulations");

    impl.updatePopulations(args...);
}


template<typename Impl_t>
void IonUpdaterProxy<Impl_t>::updateIons(auto& ions)
{
    ions.computeChargeDensity();
    ions.computeBulkVelocity();
}

template<typename... Args>
struct IonUpdaterImplResolver
{
    auto static constexpr updater_impl()
    {
        if constexpr (PHARE_UPDATER_IMPL == 0)
            return static_cast<IonUpdater0<Args...>*>(0);
        if constexpr (PHARE_UPDATER_IMPL == 1)
            return static_cast<IonUpdater1<Args...>*>(0);
        else
            throw std::runtime_error("NO IMPL FOR VALUE");
    }

    using value_type = std::decay_t<decltype(*updater_impl())>;
};

template<typename... Args>
using IonUpdaterImpl_t = IonUpdaterImplResolver<Args...>::value_type;

template<typename... Args>
using IonUpdater = IonUpdaterProxy<IonUpdaterImpl_t<Args...>>;


} // namespace PHARE::core


#endif // ION_UPDATER_HPP
