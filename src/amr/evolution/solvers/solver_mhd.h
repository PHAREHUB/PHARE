
#ifndef PHARE_SOLVER_MHD_H
#define PHARE_SOLVER_MHD_H

#include "evolution/solvers/solver.h"
#include "evolution/transactions/mhd_transaction_info.h"

namespace PHARE
{
template<typename MHDModel>
class SolverMHD : public ISolver
{
public:
    SolverMHD()
        : ISolver{"MHDSolver"}
    {
    }


    virtual ~SolverMHD() = default;

    virtual std::string modelName() const override { return MHDModel::model_name; }


    virtual void fillTransactionInfo(std::unique_ptr<ITransactionInfo> const& info) const override
    {
        //
    }


    virtual void registerResources(PhysicalModel& model) override {}

    // TODO make this a resourcesUser
    virtual void allocate(PhysicalModel& model, SAMRAI::hier::Patch& patch,
                          double const allocateTime) const override
    {
    }

    virtual void advanceLevel(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                              int const levelNumber, PhysicalModel& model,
                              ITransaction& fromCoarser, const double currentTime,
                              const double newTime) override
    {
    }
};
} // namespace PHARE


#endif
