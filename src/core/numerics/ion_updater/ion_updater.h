#ifndef PHARE_ION_UPDATER_H
#define PHARE_ION_UPDATER_H

namespace PHARE
{
namespace core
{
    enum class UpdaterMode { moments_only = 1, particles_and_moments = 2 };

    class IonUpdater
    {
    public:
        template<typename Ions, typename Electromag, typename GridLayout>
        void update(Ions& ions, Electromag const& em, GridLayout const& layout,
                    UpdaterMode = UpdaterMode::particles_and_moments);
    };


    template<typename Ions, typename Electromag, typename GridLayout>
    void IonUpdater::update(Ions& ions, Electromag const& em, GridLayout const& layout,
                            UpdaterMode mode)
    {
    }



} // namespace core
} // namespace PHARE


#endif // ION_UPDATER_H
