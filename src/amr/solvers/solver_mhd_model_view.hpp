#ifndef PHARE_SOLVER_SOLVER_MHD_MODEL_VIEW_HPP
#define PHARE_SOLVER_SOLVER_MHD_MODEL_VIEW_HPP

#include "amr/solvers/solver.hpp"
#include "amr/solvers/solver_ppc_model_view.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/numerics/constrained_transport/constrained_transport.hpp"
#include "core/numerics/godunov_fluxes/godunov_fluxes.hpp"
#include "core/numerics/finite_volume_euler/finite_volume_euler.hpp"
#include "core/numerics/primite_conservative_converter/to_conservative_converter.hpp"
#include "core/numerics/primite_conservative_converter/to_primitive_converter.hpp"

namespace PHARE::solver
{
template<typename GridLayout>
class GodunovFluxesTransformer
{
    using core_type = PHARE::core::GodunovFluxes<GridLayout>;

public:
    template<typename Layout, typename Field, typename VecField, typename... Fluxes>
    void operator()(Layout const& layouts, Field const& rho, VecField const& V, VecField const& B,
                    Field const& P, VecField const& J, Fluxes&... fluxes)
    {
        assert_equal_sizes(rho, V, B, P, J, fluxes...);
        for (std::size_t i = 0; i < layouts.size(); ++i)
        {
            auto _ = core::SetLayout(layouts[i], godunov_);
            godunov_(*rho[i], *V[i], *B[i], *P[i], *J[i], *fluxes[i]...);
        }
    }

    core_type godunov_;
};

template<typename GridLayout>
class FiniteVolumeEulerTransformer
{
    using core_type = PHARE::core::FiniteVolumeEuler<GridLayout>;

public:
    template<typename Layout, typename Field, typename... Fluxes>
    void operator()(Layout const& layouts, Field const& U, Field& Unew, double const& dt,
                    const Fluxes&... fluxes)
    {
        assert_equal_sizes(U, Unew, fluxes...);
        for (std::size_t i = 0; i < layouts.size(); ++i)
        {
            auto _ = core::SetLayout(layouts[i], finite_volume_euler_);
            finite_volume_euler_(*U[i], *Unew[i], dt, *fluxes[i]...);
        }
    }

    core_type finite_volume_euler_;
};

template<typename GridLayout>
class ConstrainedTransportTransformer
{
    using core_type = PHARE::core::ConstrainedTransport<GridLayout>;

public:
    template<typename Layout, typename VecField, typename... Fluxes>
    void operator()(Layout const& layouts, VecField& E, const Fluxes&... fluxes)
    {
        assert_equal_sizes(E, fluxes...);
        for (std::size_t i = 0; i < layouts.size(); ++i)
        {
            auto _ = core::SetLayout(layouts[i], constrained_transport_);
            constrained_transport_(*E[i], *fluxes[i]...);
        }
    }

    core_type constrained_transport_;
};

template<typename GridLayout>
class ToConservativeTransformer
{
    using core_type = PHARE::core::ToConservativeConverter<GridLayout>;

public:
    template<typename Layout, typename Field, typename VecField>
    void operator()(Layout const& layouts, const Field& rho, const VecField& V, const VecField& B,
                    const Field& P, VecField& rhoV, Field& Etot)
    {
        assert_equal_sizes(rho, V, B, P, rhoV, Etot);
        for (std::size_t i = 0; i < layouts.size(); ++i)
        {
            auto _ = core::SetLayout(layouts[i], to_conservative_);
            to_conservative_.vToRhoV(*rho[i], *V[i], *rhoV[i]);
            to_conservative_.eosPToEtot(*rho[i], *V[i], *B[i], *P[i], *Etot[i]);
        }
    }

    core_type to_conservative_;
};

template<typename GridLayout>
class ToPrimitiveTransformer
{
    using core_type = PHARE::core::ToPrimitiveConverter<GridLayout>;

public:
    template<typename Layout, typename Field, typename VecField>
    void operator()(Layout const& layouts, const Field& rho, const VecField& rhoV,
                    const VecField& B, const Field& Etot, VecField& V, Field& P)
    {
        assert_equal_sizes(rho, rhoV, B, Etot, V, P);
        for (std::size_t i = 0; i < layouts.size(); ++i)
        {
            auto _ = core::SetLayout(layouts[i], to_primtitve_);
            to_primtitve_.rhoVToV(*rho[i], *rhoV[i], *V[i]);
            to_primtitve_.eosEtotToP(*rho[i], *rhoV[i], *B[i], *Etot[i], *P[i]);
        }
    }

    core_type to_primtitve_;
};


template<typename MHDModel_>
class MHDModelView : public ISolverModelView
{
public:
    using MHDModel_t                = MHDModel_;
    using GridLayout                = typename MHDModel_t::gridlayout_type;
    using GodunovFluxes_t           = GodunovFluxesTransformer<GridLayout>;
    using Ampere_t                  = AmpereTransformer<GridLayout>;
    using FiniteVolumeEuler_t       = FiniteVolumeEulerTransformer<GridLayout>;
    using ConstrainedTransport_t    = ConstrainedTransportTransformer<GridLayout>;
    using Faraday_t                 = FaradayTransformer<GridLayout>;
    using ToPrimitiveConverter_t    = ToPrimitiveTransformer<GridLayout>;
    using ToConservativeConverter_t = ToConservativeTransformer<GridLayout>;
};

}; // namespace PHARE::solver

#endif // PHARE_SOLVER_SOLVER_MHD_MODEL_VIEW_HPP
