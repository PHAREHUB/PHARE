
#ifndef PHARE_SOLVER_MHD_H
#define PHARE_SOLVER_MHD_H

#include "messengers/mhd_messenger_info.h"
#include "physical_models/physical_model.h"
#include "solvers/solver.h"

namespace PHARE
{
namespace solver
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


        virtual void fillMessengerInfo(
            std::unique_ptr<amr::IMessengerInfo> const& info) const override
        {
            //
        }


        virtual void registerResources(IPhysicalModel& model) override {}

        // TODO make this a resourcesUser
        virtual void allocate(IPhysicalModel& model, SAMRAI::hier::Patch& patch,
                              double const allocateTime) const override
        {
        }

        virtual void advanceLevel(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                                  int const levelNumber, IPhysicalModel& model,
                                  amr::IMessenger<IPhysicalModel>& fromCoarser,
                                  const double currentTime, const double newTime) override
        {
        }
    };
} // namespace solver
} // namespace PHARE


#endif
