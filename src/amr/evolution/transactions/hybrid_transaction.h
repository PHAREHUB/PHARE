
#ifndef PHARE_HYBRID_TRANSACTION_H
#define PHARE_HYBRID_TRANSACTION_H




#include "evolution/transactions/hybrid_transaction_strategy.h"
#include "evolution/transactions/mhd_transaction.h"
#include "evolution/transactions/transaction.h"
#include "evolution/transactions/transaction_info.h"
#include "hybrid/hybrid_quantities.h"
#include "physical_models/hybrid_model.h"
#include "physical_models/physical_model.h"




namespace PHARE
{
enum class HybridTransactionStrategyType { HybridHybrid, MHDHybrid };



template<typename HybridModel>
class HybridTransaction : public ITransaction
{
private:
    using IonsT     = decltype(std::declval<HybridModel>().state.ions);
    using VecFieldT = decltype(std::declval<HybridModel>().state.electromag.E);

    using stratT = HybridTransactionStrategy<HybridModel>;

public:
    explicit HybridTransaction(std::unique_ptr<stratT> strat)
        : strat_{std::move(strat)}
    {
    }


    virtual std::unique_ptr<ITransactionInfo> emptyInfoFromCoarser() override
    {
        return strat_->emptyInfoFromCoarser();
    }

    virtual std::unique_ptr<ITransactionInfo> emptyInfoFromFiner() override
    {
        return strat_->emptyInfoFromFiner();
    }

    virtual void initLevel(int const levelNumber, double const initDataTime) const override
    {
        strat_->initLevel(levelNumber, initDataTime);
    }
    virtual std::string fineModelName() const override { return strat_->fineModelName(); }

    virtual std::string coarseModelName() const override { return strat_->coarseModelName(); }



    virtual void allocate(PhysicalModel const& model, SAMRAI::hier::Patch& patch,
                          double const allocateTime) const override
    {
        strat_->allocate(model, patch, allocateTime);
    }


    virtual void setup(std::unique_ptr<ITransactionInfo> fromCoarserInfo,
                       std::unique_ptr<ITransactionInfo> fromFinerInfo) override
    {
        strat_->setup(std::move(fromCoarserInfo), std::move(fromFinerInfo));


        // utilise les strings de la Info pour get les Ids et mettre
        // les IDs la ou il faut dans les registerRefine()
    }

    virtual void update(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                        int const levelNumber) override
    {
        strat_->update(hierarchy, levelNumber);
    }


    virtual void regrid(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                        const int levelNumber,
                        std::shared_ptr<SAMRAI::hier::PatchLevel> const& oldLevel,
                        double const initDataTime) override
    {
        strat_->regrid(hierarchy, levelNumber, oldLevel, initDataTime);
    }



    virtual void initialize(PhysicalModel const& destModel, PhysicalModel const& srcModel) override
    {
        auto const& hybDestModel = dynamic_cast<HybridModel const&>(destModel);
        strat_->initialize(hybDestModel, srcModel);
    }



    virtual std::string name() override
    {
        if (strat_ != nullptr)
        {
            return strat_->name();
        }
        else
        {
            throw std::runtime_error("invalid strat");
        }
    }



    void fillMagneticGhosts(VecFieldT& B, int const levelNumber, double const fillTime)
    {
        strat_->fillMagneticGhosts(B, levelNumber, fillTime);
    }


    void fillElectricGhosts(VecFieldT& E, int const levelNumber, double const fillTime)
    {
        strat_->fillElectricGhosts(E, levelNumber, fillTime);
    }



    void fillIonGhostParticles(IonsT& ions, int const levelNumber, double const fillTime)
    {
        strat_->fillIonGhostParticles(ions, levelNumber, fillTime);
    }



    void fillIonMomentGhosts(IonsT& ions, int const levelNumber, double const fillTime)
    {
        strat_->fillIonMomentGhosts(ions, levelNumber, fillTime);
    }



    // synchronization/coarsening methods
    void syncMagnetic(VecFieldT& B) { strat_->syncMagnetic(B); }
    void syncElectric(VecFieldT& E) { strat_->syncElectric(E); }
    void syncIonMoments(IonsT& ions) { strat_->syncIonMoments(ions); }



    virtual ~HybridTransaction() = default;

private:
    const std::unique_ptr<stratT> strat_;
};




} // namespace PHARE
#endif
