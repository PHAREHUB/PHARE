#ifndef PHARE_HYBRID_MESSENGER_STRATEGY_HPP
#define PHARE_HYBRID_MESSENGER_STRATEGY_HPP

#include "phare_mpi.hpp" // IWYU pragma: keep

#include "amr/messengers/messenger_info.hpp"

#include "amr/messengers/scheduler.hpp"

#include <SAMRAI/hier/PatchLevel.h>
#include <SAMRAI/hier/PatchHierarchy.h>


#include <string>
#include <cstdint>
#include <utility>


namespace PHARE
{
namespace amr
{
    //! hybrid communications that are not named after a single resource
    enum class HybridSchedule : std::uint8_t { IonPopBorderSum, IonBorderMax, IonGhostParticles };

    //! the scheduler key and RefinerType of each HybridSchedule
    template<HybridSchedule schedule>
    struct HybridScheduleTraits;

    template<>
    struct HybridScheduleTraits<HybridSchedule::IonPopBorderSum>
    {
        static constexpr auto type = RefinerType::PatchIonPopBorderSum;
        static inline std::string const key{"PHARE_ionPopBorderSum"};
    };

    template<>
    struct HybridScheduleTraits<HybridSchedule::IonBorderMax>
    {
        static constexpr auto type = RefinerType::PatchIonPopBorderMax;
        static inline std::string const key{"PHARE_ionBorderMax"};
    };

    template<>
    struct HybridScheduleTraits<HybridSchedule::IonGhostParticles>
    {
        static constexpr auto type = RefinerType::ExteriorGhostParticles;
        static inline std::string const key{"PHARE_ionGhostParticles"};
    };

    template<typename HybridModel>
    class HybridMessengerStrategy
    {
        using IonsT          = decltype(std::declval<HybridModel>().state.ions);
        using VecFieldT      = decltype(std::declval<HybridModel>().state.electromag.E);
        using IPhysicalModel = typename HybridModel::Interface;
        using amr_types      = HybridModel::amr_types;
        using patch_t        = amr_types::patch_t;
        using level_t        = amr_types::level_t;
        using hierarchy_t    = amr_types::hierarchy_t;
        using level_ptr      = amr_types::level_ptr;
        using hierarchy_ptr  = amr_types::hierarchy_ptr;

    public:
        /**
         * @brief allocate HybridMessengerStrategy internal resources on a given patch for a given
         * allocation time. The method is virtual and is implemented by a concrete
         * HybridMessengerStrategy
         */
        virtual void allocate(patch_t& patch, double const allocateTime) const = 0;


        /**
         * @brief setup prepares the HybridMessengerStrategy to communicate specific data between
         * two levels. The method is called by a HybridMessenger::allocate() and its concrete
         * implementation is found in concrete strategies
         */
        virtual void
        registerQuantities(std::unique_ptr<IMessengerInfo> fromCoarserInfo,
                           [[maybe_unused]] std::unique_ptr<IMessengerInfo> fromFinerInfo)
            = 0;

        virtual void registerQuantities(IPhysicalModel& coarseModel, IPhysicalModel& fineModel) {}


        virtual std::unique_ptr<IMessengerInfo> emptyInfoFromCoarser() = 0;
        virtual std::unique_ptr<IMessengerInfo> emptyInfoFromFiner()   = 0;


        virtual void registerLevel(hierarchy_ptr const& hierarchy, int const levelNumber) = 0;

        virtual void registerLevel(IPhysicalModel& model, hierarchy_ptr const& hierarchy,
                                   int const levelNumber)
        {
            registerLevel(hierarchy, levelNumber);
        }

        virtual void regrid(hierarchy_ptr const& hierarchy, int const levelNumber,
                            level_ptr const& oldLevel, IPhysicalModel& model,
                            double const initDataTime)
            = 0;

        // domain filling methods
        // used during level creation
        virtual void initLevel(IPhysicalModel& model, level_t& level, double const initDataTime)
            = 0;

        // ghost filling

        virtual void fillMagneticGhosts(VecFieldT& B, level_t const& level, double const fillTime)
            = 0;

        virtual void fillElectricGhosts(VecFieldT& E, level_t const& level, double const fillTime)
            = 0;


        virtual void fillIonGhostParticles(IonsT& ions, level_t& level, double const fillTime) = 0;


        virtual void fillIonPopMomentGhosts(IonsT& ions, level_t& level, double const afterPushTime)
            = 0;


        // virtual void fillIonMomentGhosts(IonsT& ions, level_t& level,
        //                                  double const afterPushTime)
        //     = 0;

        virtual std::string fineModelName() const = 0;

        virtual std::string coarseModelName() const = 0;

        virtual void firstStep(IPhysicalModel& model, level_t& level, double const currentTime,
                               double const prevCoarserTime, double const newCoarserTime)
            = 0;

        virtual void lastStep(IPhysicalModel& model, level_t& level) = 0;



        virtual void prepareStep(IPhysicalModel& model, level_t& level, double currentTime) = 0;


        virtual void fillRootGhosts(IPhysicalModel& model, level_t& level,
                                    double const initDataTime)
            = 0;


        virtual void synchronize(level_t& level) = 0;

        virtual void reflux(int const coarserLevelNumber, int const fineLevelNumber,
                            double const syncTime)
            = 0;

        virtual void postSynchronize(IPhysicalModel& model, level_t& level, double const time) = 0;

        virtual void fillIonPopBorders(IonsT& ions, level_t& level, double const fillTime) = 0;
        virtual void fillIonBorders(IonsT& ions, level_t& level, double const fillTime)    = 0;


        std::string name() const { return stratname_; }

        virtual ~HybridMessengerStrategy() = default;


        template<auto rtype>
        void fill(std::string const& dst, SAMRAI::hier::PatchLevel const& level, double time)
        {
            scheduler_.template fill<rtype>(dst, level, time);
        }

        auto& scheduler() { return scheduler_; }

    protected:
        explicit HybridMessengerStrategy(std::string stratName)
            : stratname_{stratName}
        {
        }

        std::string stratname_;
        RefinerScheduler scheduler_;
    };
} // namespace amr
} // namespace PHARE

#endif
