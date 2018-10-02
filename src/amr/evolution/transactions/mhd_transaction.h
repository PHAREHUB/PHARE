
#ifndef PHARE_MHD_TRANSACTION_H
#define PHARE_MHD_TRANSACTION_H

#include <memory>
#include <string>

#include <SAMRAI/hier/CoarsenOperator.h>
#include <SAMRAI/hier/RefineOperator.h>

#include "evolution/transactions/transaction.h"
#include "evolution/transactions/transaction_info.h"
#include "hybrid/hybrid_quantities.h"
#include "physical_models/mhd_model.h"
#include "physical_models/physical_model.h"

namespace PHARE
{
template<typename MHDModel>
class MHDTransaction : public ITransaction
{
public:
    MHDTransaction(std::shared_ptr<typename MHDModel::resources_manager_type> resourcesManager,
                   int const firstLevel)
        : resourcesManager_{std::move(resourcesManager)}
        , firstLevel_{firstLevel}
    {
    }

    static const std::string stratName;

    virtual std::string fineModelName() const override { return MHDModel::model_name; }

    virtual std::string coarseModelName() const override { return MHDModel::model_name; }

    virtual void allocate(PhysicalModel const& model, SAMRAI::hier::Patch& patch,
                          double const allocateTime) const override
    {
    }

    virtual void initLevel(int const levelNumber, double const initDataTime) const override {}

    virtual std::unique_ptr<ITransactionInfo> emptyInfoFromCoarser() override
    {
        return std::make_unique<MHDTransactionInfo>();
    }

    virtual std::unique_ptr<ITransactionInfo> emptyInfoFromFiner() override
    {
        return std::make_unique<MHDTransactionInfo>();
    }


    virtual void setup(std::unique_ptr<ITransactionInfo> fromCoarserInfo,
                       [[maybe_unused]] std::unique_ptr<ITransactionInfo> fromFinerInfo) override
    {
        std::unique_ptr<MHDTransactionInfo> mhdInfo{
            dynamic_cast<MHDTransactionInfo*>(fromCoarserInfo.release())};
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
    }



    virtual void initialize(PhysicalModel const& destModel, PhysicalModel const& srcModel) override
    {
    }

    virtual std::string name() override { return stratName; }

    virtual ~MHDTransaction() = default;


private:
    std::shared_ptr<typename MHDModel::resources_manager_type> resourcesManager_;
    int const firstLevel_;
};


template<typename MHDModel>
const std::string MHDTransaction<MHDModel>::stratName = "MHDModel-MHDModel";

} // namespace PHARE
#endif
