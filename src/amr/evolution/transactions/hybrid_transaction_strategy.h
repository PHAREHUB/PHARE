
#ifndef PHARE_HYBRID_TRANSACTION_STRATEGY_H
#define PHARE_HYBRID_TRANSACTION_STRATEGY_H

#include "evolution/transactions/transaction_info.h"
#include "physical_models/physical_model.h"


#include <SAMRAI/hier/PatchHierarchy.h>
#include <SAMRAI/hier/PatchLevel.h>


#include <utility>


namespace PHARE
{
template<typename HybridModel>
class HybridTransactionStrategy
{
    using IonsT     = decltype(std::declval<HybridModel>().state.ions);
    using VecFieldT = decltype(std::declval<HybridModel>().state.electromag.E);

public:
    virtual void setup(std::unique_ptr<ITransactionInfo> fromCoarserInfo,
                       [[maybe_unused]] std::unique_ptr<ITransactionInfo> fromFinerInfo)
        = 0;

    virtual std::unique_ptr<ITransactionInfo> emptyInfoFromCoarser() = 0;
    virtual std::unique_ptr<ITransactionInfo> emptyInfoFromFiner()   = 0;

    virtual void allocate(PhysicalModel const& model, SAMRAI::hier::Patch& patch,
                          double const allocateTime) const = 0;

    virtual void update(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                        int const levelNumber)
        = 0;

    virtual void
    regrid(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy, const int levelNumber,
           std::shared_ptr<SAMRAI::hier::PatchLevel> const& oldLevel, double const initDataTime)
        = 0;

    virtual void initLevel(int const levelNumber, double const initDataTime) const = 0;

    // ghost filling

    // used during solver advance
    virtual void fillMagneticGhosts(VecFieldT& B, int const levelNumber, double const fillTime) = 0;
    virtual void fillElectricGhosts(VecFieldT& E, int const levelNumber, double const fillTime) = 0;
    virtual void fillIonGhostParticles(IonsT& ions, int const levelNumber, double const fillTime)
        = 0;
    virtual void fillIonMomentGhosts(IonsT& ions, int const levelNumber, double const fillTime) = 0;


    // domain filling methods
    // used during level creation
    virtual void initialize(HybridModel const& destModel, PhysicalModel const& srcModel) = 0;

    // synchronization/coarsening methods
    virtual void syncMagnetic(VecFieldT& B)  = 0;
    virtual void syncElectric(VecFieldT& E)  = 0;
    virtual void syncIonMoments(IonsT& ions) = 0;


    virtual std::string fineModelName() const = 0;

    virtual std::string coarseModelName() const = 0;

    std::string name() const { return stratname_; }

    virtual ~HybridTransactionStrategy() = default;


protected:
    explicit HybridTransactionStrategy(std::string stratName)
        : stratname_{stratName}
    {
    }

    std::string stratname_;
};
} // namespace PHARE

#endif
