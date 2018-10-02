
#ifndef PHARE_MHD_HYBRID_TRANSACTION_STRATEGY_H
#define PHARE_MHD_HYBRID_TRANSACTION_STRATEGY_H

#include "evolution/transactions/hybrid_transaction_info.h"
#include "evolution/transactions/hybrid_transaction_strategy.h"
#include "evolution/transactions/mhd_transaction_info.h"

#include <string>

namespace PHARE
{
template<typename MHDModel, typename HybridModel>
class MHDHybridTransactionStrategy : public HybridTransactionStrategy<HybridModel>
{
    using IonsT     = decltype(std::declval<HybridModel>().state.ions);
    using VecFieldT = decltype(std::declval<HybridModel>().state.electromag.E);

public:
    static const std::string stratName;

    MHDHybridTransactionStrategy(
        std::shared_ptr<typename MHDModel::resources_manager_type> mhdResourcesManager,
        std::shared_ptr<typename HybridModel::resources_manager_type> hybridResourcesManager,
        int const firstLevel)
        : HybridTransactionStrategy<HybridModel>{stratName}
        , mhdResourcesManager_{std::move(mhdResourcesManager)}
        , hybridResourcesManager_{std::move(hybridResourcesManager)}
        , firstLevel_{firstLevel}
    {
        hybridResourcesManager_->registerResources(EM_old_);
    }


    virtual std::unique_ptr<ITransactionInfo> emptyInfoFromCoarser() override
    {
        return std::make_unique<MHDTransactionInfo>();
    }

    virtual std::unique_ptr<ITransactionInfo> emptyInfoFromFiner() override
    {
        return std::make_unique<HybridTransactionInfo>();
    }

    virtual void allocate(PhysicalModel const& model, SAMRAI::hier::Patch& patch,
                          double const allocateTime) const override
    {
        auto& hybModel = dynamic_cast<HybridModel const&>(model);

        // hybModel.resourcesManager->allocate(EM_old_.E, patch, allocateTime);
        // hybModel.resourcesManager->allocate(EM_old_.B, patch, allocateTime);
        hybModel.resourcesManager->allocate(EM_old_, patch, allocateTime);
    }


    virtual void update(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                        int const levelNumber) override
    {
    }


    virtual void regrid(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                        const int levelNumber,
                        std::shared_ptr<SAMRAI::hier::PatchLevel> const& oldLevel,
                        double const initDataTime) override
    {
        //
    }




    virtual void setup(std::unique_ptr<ITransactionInfo> fromCoarserInfo,
                       [[maybe_unused]] std::unique_ptr<ITransactionInfo> fromFinerInfo) override
    {
    }

    virtual std::string fineModelName() const override { return HybridModel::model_name; }

    virtual std::string coarseModelName() const override { return MHDModel::model_name; }


    virtual void initialize(HybridModel const& destModel, PhysicalModel const& srcModel) override
    {
        auto const& mhdSrcModel = dynamic_cast<MHDModel const&>(srcModel);
    }


    virtual void initLevel(int const levelNumber, double const initDataTime) const override {}

    virtual ~MHDHybridTransactionStrategy() = default;


    virtual void fillMagneticGhosts(VecFieldT& B, int const levelNumber,
                                    double const fillTime) override
    {
    }
    virtual void fillElectricGhosts(VecFieldT& E, int const levelNumber,
                                    double const fillTime) override
    {
    }
    virtual void fillIonGhostParticles(IonsT& ions, int const levelNumber,
                                       double const fillTime) override
    {
    }
    virtual void fillIonMomentGhosts(IonsT& ions, int const levelNumber,
                                     double const fillTime) override
    {
    }



    // synchronization/coarsening methods
    virtual void syncMagnetic(VecFieldT& B) override {}
    virtual void syncElectric(VecFieldT& E) override {}
    virtual void syncIonMoments(IonsT& ions) override {}



private:
    using Electromag = decltype(std::declval<HybridModel>().state.electromag);

    std::shared_ptr<typename MHDModel::resources_manager_type> mhdResourcesManager_;
    std::shared_ptr<typename HybridModel::resources_manager_type> hybridResourcesManager_;
    int const firstLevel_;
    Electromag EM_old_{"EM_old"};
};

template<typename MHDModel, typename HybridModel>
const std::string MHDHybridTransactionStrategy<MHDModel, HybridModel>::stratName
    = "MHDModel-HybridModel";



} // namespace PHARE


#endif
