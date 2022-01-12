
#ifndef PHARE_SOLVER_MHD_HPP
#define PHARE_SOLVER_MHD_HPP

#include "amr/messengers/mhd_messenger_info.hpp"
#include "amr/physical_models/physical_model.hpp"
#include "amr/solvers/solver.hpp"

namespace PHARE
{
namespace solver
{
    template<typename MHDModel, typename AMR_Types>
    class SolverMHD : public ISolver<AMR_Types>
    {
    public:
        SolverMHD()
            : ISolver<AMR_Types>{"MHDSolver"}
        {
        }


        virtual ~SolverMHD() = default;

        virtual std::string modelName() const override { return MHDModel::model_name; }


        virtual void
        fillMessengerInfo(std::unique_ptr<amr::IMessengerInfo> const& /*info*/) const override
        {
        }


        virtual void registerResources(IPhysicalModel<AMR_Types>& /*model*/) override {}

        // TODO make this a resourcesUser
        virtual void allocate(IPhysicalModel<AMR_Types>& /*model*/, SAMRAI::hier::Patch& /*patch*/,
                              double const /*allocateTime*/) const override
        {
        }

        virtual void
        advanceLevel(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& /*hierarchy*/,
                     int const /*levelNumber*/, IPhysicalModel<AMR_Types>& /*model*/,
                     amr::IMessenger<IPhysicalModel<AMR_Types>>& /*fromCoarser*/,
                     const double /*currentTime*/, const double /*newTime*/) override
        {
        }
    };
} // namespace solver
} // namespace PHARE


#endif
