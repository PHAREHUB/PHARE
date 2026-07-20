#ifndef PHARE_CORE_NUMERICS_TIME_INTEGRATOR_COMPUTE_FLUXES_HPP
#define PHARE_CORE_NUMERICS_TIME_INTEGRATOR_COMPUTE_FLUXES_HPP

#include "initializer/data_provider.hpp"
#include "core/inner_boundary/inner_bc_context.hpp"
#include "core/inner_boundary/inner_boundary_defs.hpp"
#include "core/numerics/godunov_fluxes/godunov_utils.hpp"
#include "amr/resources_manager/amr_utils.hpp"
#include "amr/solvers/solver_mhd_field_evolvers.hpp"

namespace PHARE::solver
{
template<typename FVMethodStrategy, typename MHDModel>
class ComputeFluxes
{
    using level_t                     = MHDModel::level_t;
    using gridlayout_type             = MHDModel::gridlayout_type;
    using state_type                  = MHDModel::state_type;
    using resources_manager_type      = MHDModel::resources_manager_type;
    using inner_boundary_manager_type = MHDModel::inner_boundary_manager_type;
    using Dispatchers_t               = Dispatchers<MHDModel>;

    using Ampere_t = Dispatchers_t::Ampere_t;

    using FVMethod_t     = Dispatchers_t::template FVMethod_t<FVMethodStrategy>;
    using FVMethodInfo_t = FVMethod_t::info_type;

    constexpr static auto Hall             = FVMethod_t::Hall;
    constexpr static auto Resistivity      = FVMethod_t::Resistivity;
    constexpr static auto HyperResistivity = FVMethod_t::HyperResistivity;

    template<typename T>
    using Rec = FVMethod_t::template Rec<T>;

    using ConstrainedTransport_t
        = Dispatchers_t::template ConstrainedTransport_t<Rec, Hall, Resistivity, HyperResistivity>;
    using ConstrainedTransportInfo_t = ConstrainedTransport_t::info_type;

    using ToPrimitiveConverter_t    = Dispatchers_t::ToPrimitiveConverter_t;
    using ToConservativeConverter_t = Dispatchers_t::ToConservativeConverter_t;

    using VecField    = MHDModel::vecfield_type;
    using Equations_t = FVMethod_t::Equations_t;


public:
    ComputeFluxes(PHARE::initializer::PHAREDict const& dict)
        : fVMethodInfo_{FVMethodInfo_t::FROM(dict["fv_method"])}
        , constrainedTransportInfo_{ConstrainedTransportInfo_t::FROM(dict["constrained_transport"])}
        , to_primitive_gamma_{dict["to_primitive"]["heat_capacity_ratio"]}
        , to_conservative_gamma_{dict["to_conservative"]["heat_capacity_ratio"]}
    {
    }

    void operator()(MHDModel& model, auto& state, auto& fluxes, auto& bc, level_t& level,
                    double const newTime)
    {
        ToPrimitiveConverter_t{level, model}(state, to_primitive_gamma_, newTime);

        if constexpr (Hall || Resistivity || HyperResistivity)
        {
            Ampere_t{level, model}(state.B, state.J);
            TimeSetter{level, model, newTime}(state.B, state.J);
        }

        FVMethod_t{level, model, fVMethodInfo_}(fvm_, ct_, state, fluxes, newTime);

        // unecessary if we decide to store both primitive and conservative variables
        ToConservativeConverter_t{level, model}(state, to_conservative_gamma_, newTime);

        ConstrainedTransport_t{level, model, constrainedTransportInfo_}(ct_, state);

        if (model.hasInnerBoundary())
        {
            resources_manager_type& rm       = *model.resourcesManager;
            inner_boundary_manager_type& ibm = *model.innerBoundaryManager;
            core::InnerBCContext<state_type> ctx{state, state, newTime};
            amr::visitLevel<gridlayout_type>(
                level, rm,
                [&](auto& layout, auto&&, auto&&) {
                    ibm.applyBC(state.E, layout, ctx);

                    // Pin E to 0 in inactive cells (deep inside the body) so the large CT
                    // values there neither pollute diagnostics nor leak into the cut-cell
                    // Faraday stencil that consumes E next. Done per component centering.
                    auto& meshData = ibm.getMeshData();
                    auto zeroE     = [&](auto component) {
                        auto centering = layout.centering(state.E(component).physicalQuantity());
                        auto& status   = meshData.getStatusFieldFromCentering(centering);
                        layout.evalOnBox(state.E(component), [&](auto&... args) {
                            auto idx = core::MeshIndex<gridlayout_type::dimension>{args...};
                            if (status(idx) == core::toDouble(core::ElemStatus::Inactive))
                                state.E(component)(idx) = 0.0;
                        });
                    };
                    zeroE(core::Component::X);
                    zeroE(core::Component::Y);
                    zeroE(core::Component::Z);
                },
                ibm, state);
        }
    }

    void registerResources(MHDModel& model)
    {
        model.resourcesManager->registerResources(fvm_);
        model.resourcesManager->registerResources(ct_);
    }

    void allocate(MHDModel& model, auto& patch, double const allocateTime) const
    {
        model.resourcesManager->allocate(fvm_, patch, allocateTime);
        model.resourcesManager->allocate(ct_, patch, allocateTime);
    }

private:
    FVMethodInfo_t fVMethodInfo_;
    ConstrainedTransportInfo_t constrainedTransportInfo_;

    // Ampere_t ampere_;
    core::GodunovState<VecField, Equations_t> fvm_{};
    core::UpwindConstrainedTransportState<VecField, Hall, Resistivity> ct_{};
    // ToPrimitiveConverter_t to_primitive_;
    // ToConservativeConverter_t to_conservative_;
    double to_primitive_gamma_;
    double to_conservative_gamma_;
};
} // namespace PHARE::solver

#endif
