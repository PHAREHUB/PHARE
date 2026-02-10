#ifndef PHARE_HYBRID_MESSENGER_STRATEGY_HPP
#define PHARE_HYBRID_MESSENGER_STRATEGY_HPP

#include "core/def/phare_mpi.hpp" // IWYU pragma: keep

#include "amr/messengers/messenger_info.hpp"

#include <SAMRAI/hier/PatchLevel.h>
#include <SAMRAI/hier/PatchHierarchy.h>


#include <utility>


namespace PHARE
{
namespace amr
{
    template<typename HybridModel>
    class HybridMessengerStrategy
    {
        using IonsT          = decltype(std::declval<HybridModel>().state.ions);
        using VecFieldT      = decltype(std::declval<HybridModel>().state.electromag.E);
        using IPhysicalModel = typename HybridModel::Interface;

    public:
        /**
         * @brief allocate HybridMessengerStrategy internal resources on a given patch for a given
         * allocation time. The method is virtual and is implemented by a concrete
         * HybridMessengerStrategy
         */
        virtual void allocate(SAMRAI::hier::Patch& patch, double const allocateTime) const = 0;


        /**
         * @brief setup prepares the HybridMessengerStrategy to communicate specific data between
         * two levels. The method is called by a HybridMessenger::allocate() and its concrete
         * implementation is found in concrete strategies
         */
        virtual void
        registerQuantities(std::unique_ptr<IMessengerInfo> fromCoarserInfo,
                           [[maybe_unused]] std::unique_ptr<IMessengerInfo> fromFinerInfo)
            = 0;



        virtual std::unique_ptr<IMessengerInfo> emptyInfoFromCoarser() = 0;
        virtual std::unique_ptr<IMessengerInfo> emptyInfoFromFiner()   = 0;


        virtual void registerLevel(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                                   int const levelNumber)
            = 0;

        virtual void regrid(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                            int const levelNumber,
                            std::shared_ptr<SAMRAI::hier::PatchLevel> const& oldLevel,
                            IPhysicalModel& model, double const initDataTime)
            = 0;

        // domain filling methods
        // used during level creation
        virtual void initLevel(IPhysicalModel& model, SAMRAI::hier::PatchLevel& level,
                               double const initDataTime)
            = 0;

        // ghost filling

        virtual void fillMagneticGhosts(VecFieldT& B, SAMRAI::hier::PatchLevel const& level,
                                        double const fillTime)
            = 0;

        virtual void fillElectricGhosts(VecFieldT& E, SAMRAI::hier::PatchLevel const& level,
                                        double const fillTime)
            = 0;


        virtual void fillCurrentGhosts(VecFieldT& J, SAMRAI::hier::PatchLevel const& level,
                                       double const fillTime)
            = 0;


        virtual void fillIonGhostParticles(IonsT& ions, SAMRAI::hier::PatchLevel& level,
                                           double const fillTime)
            = 0;


        virtual void fillIonPopMomentGhosts(IonsT& ions, SAMRAI::hier::PatchLevel& level,
                                            double const afterPushTime)
            = 0;


        // virtual void fillIonMomentGhosts(IonsT& ions, SAMRAI::hier::PatchLevel& level,
        //                                  double const afterPushTime)
        //     = 0;

        virtual std::string fineModelName() const = 0;

        virtual std::string coarseModelName() const = 0;

        virtual void firstStep(IPhysicalModel& model, SAMRAI::hier::PatchLevel& level,
                               std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                               double const currentTime, double const prevCoarserTime,
                               double const newCoarserTime)
            = 0;

        virtual void lastStep(IPhysicalModel& model, SAMRAI::hier::PatchLevel& level) = 0;



        virtual void prepareStep(IPhysicalModel& model, SAMRAI::hier::PatchLevel& level,
                                 double currentTime)
            = 0;


        virtual void fillRootGhosts(IPhysicalModel& model, SAMRAI::hier::PatchLevel& level,
                                    double const initDataTime)
            = 0;


        virtual void synchronize(SAMRAI::hier::PatchLevel& level) = 0;

        virtual void reflux(int const coarserLevelNumber, int const fineLevelNumber,
                            double const syncTime)
            = 0;

        virtual void postSynchronize(IPhysicalModel& model, SAMRAI::hier::PatchLevel& level,
                                     double const time)
            = 0;

        virtual void fillFluxBorders(IonsT& ions, SAMRAI::hier::PatchLevel& level,
                                     double const fillTime)
            = 0;
        virtual void fillDensityBorders(IonsT& ions, SAMRAI::hier::PatchLevel& level,
                                        double const fillTime)
            = 0;
        virtual void fillIonBorders(IonsT& ions, SAMRAI::hier::PatchLevel& level,
                                    double const fillTime)
            = 0;


        std::string name() const { return stratname_; }

        virtual ~HybridMessengerStrategy() = default;


    protected:
        explicit HybridMessengerStrategy(std::string stratName)
            : stratname_{stratName}
        {
        }

        std::string stratname_;
    };
} // namespace amr
} // namespace PHARE

#endif
