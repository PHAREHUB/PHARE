#ifndef PHARE_ION_UPDATERS_HPP
#define PHARE_ION_UPDATERS_HPP

#include "core/def/phare_config.hpp"
#include "core/numerics/ion_updater/ion_updater.hpp"          // IWYU pragma: keep
#include "core/numerics/ion_updater/ion_updater_multi_ts.hpp" // IWYU pragma: keep
// #include "core/numerics/ion_updater/ion_updater_multi_pc.hpp" // IWYU pragma: keep


namespace PHARE::core
{


template<typename Ions, typename Electromag, typename GridLayout>
struct IonUpdaterImplResolverFns
{
    using ParticleArray_t = Ions::particle_array_type;

    auto static constexpr updater_impl()
    {
        using enum LayoutMode;
        if constexpr (any_in(ParticleArray_t::layout_mode, AoSMapped))
            return _as_nullptr_<IonUpdater<Ions, Electromag, GridLayout>*>();
#if PHARE_HAVE_MKN_GPU
        else if constexpr (is_tiled(ParticleArray_t::layout_mode))
            return _as_nullptr_<mkn::IonUpdaterMultiTS<Ions, Electromag, GridLayout>*>();
#endif // PHARE_HAVE_MKN_GPU
        else
            static_assert(dependent_false_v<Ions>);
    }

    template<typename T>
    auto static constexpr _as_nullptr_()
    {
        return T{nullptr};
    }
};

template<typename Ions, typename Electromag, typename GridLayout>
struct IonUpdaterImplResolver : IonUpdaterImplResolverFns<Ions, Electromag, GridLayout>
{
    using Super        = IonUpdaterImplResolverFns<Ions, Electromag, GridLayout>;
    using IonUpdater_t = std::decay_t<decltype(*Super::updater_impl())>;
};


} // namespace PHARE::core


#endif // ION_UPDATER_HPP
