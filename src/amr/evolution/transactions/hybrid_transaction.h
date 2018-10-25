
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


    /**
     * @brief allocate calls the abstract HybridTransactionStrategy to perform the allocation of its
     * internal resources
     */
    virtual void allocate(SAMRAI::hier::Patch& patch, double const allocateTime) const override
    {
        strat_->allocate(patch, allocateTime);
    }



    /**
     * @brief setup prepares a HybridTransaction to communicates variables given in the two
     * information structures. Since the HybridTransaction does not know which concrete strategy it
     * is using, it cannot directly use these informations. It thus passes them to its strategy.
     *
     * @param fromCoarserInfo see ITransaction
     * @param fromFinerInfo see ITransaction
     */
    virtual void setup(std::unique_ptr<ITransactionInfo> fromCoarserInfo,
                       std::unique_ptr<ITransactionInfo> fromFinerInfo) override
    {
        strat_->setup(std::move(fromCoarserInfo), std::move(fromFinerInfo));
    }




    /**
     * @brief setLevel is the hybrid concrete transaction of the virtual method
     * ITransaction::setLevel() it calls the HybridTransactionStrategy::setLevel.
     */
    virtual void setLevel(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                          int const levelNumber) override
    {
        strat_->setLevel(hierarchy, levelNumber);
    }



    /**
     * @brief regrid performs regridding communications for a hybrid transaction.
     * Since the transaction does not know what the coarser model is, it forwards the call to its
     * strategy.
     */
    virtual void regrid(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                        const int levelNumber,
                        std::shared_ptr<SAMRAI::hier::PatchLevel> const& oldLevel,
                        double const initDataTime) final
    {
        strat_->regrid(hierarchy, levelNumber, oldLevel, initDataTime);
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


    virtual void firstStep(PhysicalModel const& model) final {}


    virtual void lastStep(PhysicalModel const& model) final {}




    virtual std::string fineModelName() const override { return strat_->fineModelName(); }

    virtual std::string coarseModelName() const override { return strat_->coarseModelName(); }




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
