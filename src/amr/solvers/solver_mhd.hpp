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
        using patch_t     = typename AMR_Types::patch_t;
        using level_t     = typename AMR_Types::level_t;
        using hierarchy_t = typename AMR_Types::hierarchy_t;

        SolverMHD()
            : ISolver<AMR_Types>{"MHDSolver"}
        {
        }


        virtual ~SolverMHD() = default;

        std::string modelName() const override { return MHDModel::model_name; }


        void fillMessengerInfo(std::unique_ptr<amr::IMessengerInfo> const& /*info*/) const override
        {
        }


        void registerResources(IPhysicalModel<AMR_Types>& /*model*/) override {}

        // TODO make this a resourcesUser
        void allocate(IPhysicalModel<AMR_Types>& /*model*/, patch_t& /*patch*/,
                      double const /*allocateTime*/) const override
        {
        }

        void advanceLevel(hierarchy_t const& /*hierarchy*/, int const /*levelNumber*/,
                          IPhysicalModel<AMR_Types>& /*model*/,
                          amr::IMessenger<IPhysicalModel<AMR_Types>>& /*fromCoarser*/,
                          const double /*currentTime*/, const double /*newTime*/) override
        {
        }

        // TORM
        // std::shared_ptr<ISolverModelView> make_view(level_t&, IPhysicalModel<AMR_Types>&)
        // override
        // {
        //     throw std::runtime_error("Not implemented in mhd solver");
        //     return nullptr;
        // }
    };
} // namespace solver
} // namespace PHARE


#endif
