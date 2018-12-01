
#ifndef PHARE_HYBRID_MESSENGER_STRATEGY_H
#define PHARE_HYBRID_MESSENGER_STRATEGY_H

#include "evolution/messengers/messenger_info.h"
#include "physical_models/physical_model.h"


#include <SAMRAI/hier/PatchHierarchy.h>
#include <SAMRAI/hier/PatchLevel.h>


#include <utility>


namespace PHARE
{
template<typename HybridModel>
class HybridMessengerStrategy
{
    using IonsT     = decltype(std::declval<HybridModel>().state.ions);
    using VecFieldT = decltype(std::declval<HybridModel>().state.electromag.E);

public:
    /**
     * @brief allocate HybridMessengerStrategy internal resources on a given patch for a given
     * allocation time. The method is virtual and is implemented by a concrete
     * HybridMessengerStrategy
     */
    virtual void allocate(SAMRAI::hier::Patch& patch, double const allocateTime) const = 0;


    /**
     * @brief setup prepares the HybridMessengerStrategy to communicate specific data between two
     * levels. The method is called by a HybridMessenger::allocate() and its concrete
     * implementation is found in concrete strategies
     */
    virtual void registerQuantities(std::unique_ptr<IMessengerInfo> fromCoarserInfo,
                                    [[maybe_unused]] std::unique_ptr<IMessengerInfo> fromFinerInfo)
        = 0;



    virtual std::unique_ptr<IMessengerInfo> emptyInfoFromCoarser() = 0;
    virtual std::unique_ptr<IMessengerInfo> emptyInfoFromFiner()   = 0;


    virtual void registerLevel(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                               int const levelNumber)
        = 0;

    virtual void
    regrid(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy, const int levelNumber,
           std::shared_ptr<SAMRAI::hier::PatchLevel> const& oldLevel, double const initDataTime)
        = 0;

    // domain filling methods
    // used during level creation
    virtual void initLevel(int const levelNumber, double const initDataTime) const = 0;

    // ghost filling

    // used during solver advance
    virtual void fillMagneticGhosts(VecFieldT& B, int const levelNumber, double const fillTime) = 0;
    virtual void fillElectricGhosts(VecFieldT& E, int const levelNumber, double const fillTime) = 0;
    virtual void fillIonGhostParticles(IonsT& ions, int const levelNumber, double const fillTime)
        = 0;
    virtual void fillIonMomentGhosts(IonsT& ions, int const levelNumber, double const fillTime) = 0;



    // synchronization/coarsening methods
    virtual void syncMagnetic(VecFieldT& B)  = 0;
    virtual void syncElectric(VecFieldT& E)  = 0;
    virtual void syncIonMoments(IonsT& ions) = 0;


    virtual std::string fineModelName() const = 0;

    virtual std::string coarseModelName() const = 0;


    virtual void lastStep(IPhysicalModel& model, SAMRAI::hier::PatchLevel& level) = 0;


    std::string name() const { return stratname_; }

    virtual ~HybridMessengerStrategy() = default;


protected:
    explicit HybridMessengerStrategy(std::string stratName)
        : stratname_{stratName}
    {
    }

    std::string stratname_;
};
} // namespace PHARE

#endif
