
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
        using patch_t          = AMR_Types::patch_t;
        using level_t          = AMR_Types::level_t;
        using hierarchy_t      = AMR_Types::hierarchy_t;
        using IPhysicalModel_t = IPhysicalModel<AMR_Types>;

        SolverMHD()
            : ISolver<AMR_Types>{"MHDSolver"}
        {
        }


        virtual ~SolverMHD() = default;

        std::string modelName() const override { return MHDModel::model_name; }


        void fillMessengerInfo(std::unique_ptr<amr::IMessengerInfo> const& /*info*/) const override
        {
        }


        void registerResources(IPhysicalModel_t& /*model*/) override {}

        // TODO make this a resourcesUser
        void allocate(IPhysicalModel_t& /*model*/, patch_t& /*patch*/,
                      double const /*allocateTime*/) const override
        {
        }

        void prepareStep(IPhysicalModel_t& model, SAMRAI::hier::PatchLevel& level,
                         double const currentTime) override
        {
        }

        void accumulateFluxSum(IPhysicalModel_t& model, SAMRAI::hier::PatchLevel& level,
                               double const coef) override
        {
        }

        void resetFluxSum(IPhysicalModel_t& model, SAMRAI::hier::PatchLevel& level) override {}

        virtual void reflux(IPhysicalModel_t& model, SAMRAI::hier::PatchLevel& level,
                            amr::IMessenger<IPhysicalModel_t>& messenger,
                            double const time) override
        {
        }

        void advanceLevel(hierarchy_t const& /*hierarchy*/, int const /*levelNumber*/,
                          IPhysicalModel_t& /*view*/,
                          amr::IMessenger<IPhysicalModel_t>& /*fromCoarser*/,
                          double const /*currentTime*/, double const /*newTime*/) override
        {
        }
    };
} // namespace solver
} // namespace PHARE


#endif
