#ifndef PHARE_MHD_STATE_HPP
#define PHARE_MHD_STATE_HPP

#include "core/def.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/data/vecfield/vecfield_component.hpp"
#include "core/numerics/primite_conservative_converter/to_conservative_converter.hpp"
#include "core/data/field/initializers/field_user_initializer.hpp"
#include "core/data/vecfield/vecfield_initializer.hpp"
#include "core/mhd/mhd_quantities.hpp"
#include "core/models/physical_state.hpp"

#include "initializer/data_provider.hpp"

#include <stdexcept>

namespace PHARE
{
namespace core
{
    // Fill Bout = curl(A) over the ghost box, where A = (Ax, Ay, Az) is a full 3D vector potential.
    // A is sampled at E (edge) centring into the A scratch vecfield, then the discrete `deriv` (the
    // same operator Faraday uses) yields a face-centred B with discrete div B = 0 to machine
    // precision (div . curl = 0), unlike a component-wise init sampled at face centres. A is the
    // caller's E buffer (overwritten before its real use by the constrained transport).
    //   Bx = dAz/dy - dAy/dz   (2D: dAz/dy ; 1D: 0)
    //   By = dAx/dz - dAz/dx   (1D/2D: -dAz/dx)
    //   Bz = dAy/dx - dAx/dy   (1D: dAy/dx)
    // The per-dimension terms mirror Faraday's curl(E); dropping a purely-out-of-plane A (only Az,
    // 2D) recovers the legacy (dAz/dy, -dAz/dx, 0).
    template<typename GridLayout, typename VecField>
    void initBFromPotential(VecFieldInitializer<GridLayout::dimension>& aInit, VecField& Bout,
                            VecField& A, GridLayout const& layout)
    {
        aInit.initialize(A, layout);

        auto const& Ax = A(Component::X);
        auto const& Ay = A(Component::Y);
        auto const& Az = A(Component::Z);

        auto& Bx = Bout(Component::X);
        auto& By = Bout(Component::Y);
        auto& Bz = Bout(Component::Z);

        // Fill the full ghost box so B0 ghosts (read by the flux reconstruction and refreshed each
        // substep by updateExternalFields) and B1 ghosts are not left NaN. A is filled over its own
        // ghost box, so the discrete deriv is exact on the domain; the outermost ghost layer is
        // approximate but finite.
        layout.evalOnGhostBox(Bx, [&](auto&... args) {
            if constexpr (GridLayout::dimension == 1)
                Bx(args...) = 0.0;
            else if constexpr (GridLayout::dimension == 2)
                Bx(args...) = layout.template deriv<Direction::Y>(Az, {args...});
            else
                Bx(args...) = layout.template deriv<Direction::Y>(Az, {args...})
                              - layout.template deriv<Direction::Z>(Ay, {args...});
        });
        layout.evalOnGhostBox(By, [&](auto&... args) {
            if constexpr (GridLayout::dimension == 3)
                By(args...) = layout.template deriv<Direction::Z>(Ax, {args...})
                              - layout.template deriv<Direction::X>(Az, {args...});
            else
                By(args...) = -layout.template deriv<Direction::X>(Az, {args...});
        });
        layout.evalOnGhostBox(Bz, [&](auto&... args) {
            if constexpr (GridLayout::dimension == 1)
                Bz(args...) = layout.template deriv<Direction::X>(Ay, {args...});
            else
                Bz(args...) = layout.template deriv<Direction::X>(Ay, {args...})
                              - layout.template deriv<Direction::Y>(Ax, {args...});
        });
    }

    // Time-dependent overload: A = A(x,t) is sampled at the given time before taking the discrete
    // curl. Applied to both A0(t) -> B0 = curl(A0) and dA0/dt(t) -> dB0/dt = curl(dA0/dt), so the
    // time derivative of B0 uses exactly the same discrete curl as B0 itself (hence discretely
    // div-free and consistent). The body is identical to the static overload aside from the timed
    // sampling of A.
    template<typename GridLayout, typename VecField>
    void initBFromPotential(SpaceTimeVecFieldInitializer<GridLayout::dimension>& aInit,
                            VecField& Bout, VecField& A, GridLayout const& layout, double time)
    {
        aInit.initialize(A, layout, time);

        auto const& Ax = A(Component::X);
        auto const& Ay = A(Component::Y);
        auto const& Az = A(Component::Z);

        auto& Bx = Bout(Component::X);
        auto& By = Bout(Component::Y);
        auto& Bz = Bout(Component::Z);

        layout.evalOnGhostBox(Bx, [&](auto&... args) {
            if constexpr (GridLayout::dimension == 1)
                Bx(args...) = 0.0;
            else if constexpr (GridLayout::dimension == 2)
                Bx(args...) = layout.template deriv<Direction::Y>(Az, {args...});
            else
                Bx(args...) = layout.template deriv<Direction::Y>(Az, {args...})
                              - layout.template deriv<Direction::Z>(Ay, {args...});
        });
        layout.evalOnGhostBox(By, [&](auto&... args) {
            if constexpr (GridLayout::dimension == 3)
                By(args...) = layout.template deriv<Direction::Z>(Ax, {args...})
                              - layout.template deriv<Direction::X>(Az, {args...});
            else
                By(args...) = -layout.template deriv<Direction::X>(Az, {args...});
        });
        layout.evalOnGhostBox(Bz, [&](auto&... args) {
            if constexpr (GridLayout::dimension == 1)
                Bz(args...) = layout.template deriv<Direction::X>(Ay, {args...});
            else
                Bz(args...) = layout.template deriv<Direction::X>(Ay, {args...})
                              - layout.template deriv<Direction::Y>(Ax, {args...});
        });
    }

    // The dynamic MHD state of the B = B0 + B1 split. It holds the evolved perturbation field B1
    // (the user prescribes the total field B; B1 = B - B0 at init) and the conserved perturbation
    // energy Etot1 (kinetic + thermal + 1/2|B1|^2). The static background B0 is NOT stored here:
    // it is a single face-centered vector field held by the model and shared by every RK stage.
    template<typename VecFieldT>
    class MHDState : public IPhysicalState
    {
    public:
        using vecfield_type = VecFieldT;
        using field_type    = typename VecFieldT::field_type;

        static constexpr auto dimension = VecFieldT::dimension;

        //-------------------------------------------------------------------------
        //                  start the ResourcesUser interface
        //-------------------------------------------------------------------------

        NO_DISCARD bool isUsable() const
        {
            return rho.isUsable() and V.isUsable() and B1.isUsable() and B.isUsable()
                   and P.isUsable() and rhoV.isUsable() and Etot1.isUsable() and J.isUsable()
                   and E.isUsable();
        }

        NO_DISCARD bool isSettable() const
        {
            return rho.isSettable() and V.isSettable() and B1.isSettable() and B.isSettable()
                   and P.isSettable() and rhoV.isSettable() and Etot1.isSettable()
                   and J.isSettable() and E.isSettable();
        }

        NO_DISCARD auto getCompileTimeResourcesViewList() const
        {
            return std::forward_as_tuple(rho, V, B1, B, P, rhoV, Etot1, J, E);
        }

        NO_DISCARD auto getCompileTimeResourcesViewList()
        {
            return std::forward_as_tuple(rho, V, B1, B, P, rhoV, Etot1, J, E);
        }

        //-------------------------------------------------------------------------
        //                  ends the ResourcesUser interface
        //-------------------------------------------------------------------------

        MHDState(PHARE::initializer::PHAREDict const& dict)
            : rho{dict["name"].template to<std::string>() + "_" + "rho", MHDQuantity::Scalar::rho}
            , V{dict["name"].template to<std::string>() + "_" + "V", MHDQuantity::Vector::V}
            , B1{dict["name"].template to<std::string>() + "_" + "B1", MHDQuantity::Vector::B1}
            , B{dict["name"].template to<std::string>() + "_" + "B", MHDQuantity::Vector::B}
            , P{dict["name"].template to<std::string>() + "_" + "P", MHDQuantity::Scalar::P}


            , rhoV{dict["name"].template to<std::string>() + "_" + "rhoV",
                   MHDQuantity::Vector::rhoV}
            , Etot1{dict["name"].template to<std::string>() + "_" + "Etot1",
                    MHDQuantity::Scalar::Etot1}


            , E{dict["name"].template to<std::string>() + "_" + "E", MHDQuantity::Vector::E}
            , J{dict["name"].template to<std::string>() + "_" + "J", MHDQuantity::Vector::J}


            , rhoinit_{dict["density"]["initializer"]
                           .template to<initializer::InitFunction<dimension>>()}
            , Vinit_{dict["velocity"]["initializer"]}
            , totalBInit_{dict["magnetic"]["initializer"]}
            , Pinit_{dict["pressure"]["initializer"]
                         .template to<initializer::InitFunction<dimension>>()}
            , gamma_{dict["to_conservative_init"]["heat_capacity_ratio"].template to<double>()}
            , b0FromPotential_{cppdict::get_value(dict, "b0_init_mode", std::string{"components"})
                               == "potential"}
            , b1FromPotential_{cppdict::get_value(dict, "b1_init_mode", std::string{"components"})
                               == "potential"}
        {
            // Vector-potential init: the potential is read only when B1 is in "potential" mode, so
            // dicts that predate this feature (e.g. C++ unit-test dicts) need not provide the
            // potential key. b0FromPotential_ is needed only to skip the B1 = total - B0
            // subtraction below (the Python total cannot fold a potential B0).
            if (b1FromPotential_)
                a1Init_ = VecFieldInitializer<dimension>{
                    dict["perturbed_magnetic"]["initializer"]["potential"]};
        }

        MHDState(std::string name)
            : rho{name + "_" + "rho", MHDQuantity::Scalar::rho}
            , V{name + "_" + "V", MHDQuantity::Vector::V}
            , B1{name + "_" + "B1", MHDQuantity::Vector::B1}
            , B{name + "_" + "B", MHDQuantity::Vector::B}
            , P{name + "_" + "P", MHDQuantity::Scalar::P}


            , rhoV{name + "_" + "rhoV", MHDQuantity::Vector::rhoV}
            , Etot1{name + "_" + "Etot1", MHDQuantity::Scalar::Etot1}


            , E{name + "_" + "E", MHDQuantity::Vector::E}
            , J{name + "_" + "J", MHDQuantity::Vector::J}

            , gamma_{}
        {
        }

        // Initialize the dynamic state. The user prescribes the TOTAL field B: it is initialized
        // into B1, from which the (already-evaluated) background B0 is subtracted so B1 holds the
        // perturbation. The conserved state is then formed (Etot1 uses B1 only).
        template<typename GridLayout, typename VecField>
        void initialize(GridLayout const& layout, VecField const& B0)
        {
            FieldUserFunctionInitializer::initialize(rho, layout, rhoinit_);
            Vinit_.initialize(V, layout);
            FieldUserFunctionInitializer::initialize(P, layout, Pinit_);

            if (b1FromPotential_)
            {
                // B1 = curl(A1), built directly as the perturbation (no subtraction). E is the
                // (edge-centred) A scratch buffer, cleared below before its real use.
                initBFromPotential(a1Init_, B1, E, layout);
            }
            else
            {
                totalBInit_.initialize(B1, layout);
                // In component mode dict["magnetic"] holds the total field B0 + B1, so subtract the
                // analytic B0 to recover B1. When B0 itself comes from a vector potential the Python
                // side could not fold it into the total, so dict["magnetic"] already equals the
                // perturbation and no subtraction is needed.
                if (!b0FromPotential_)
                    for (auto const& component : {Component::X, Component::Y, Component::Z})
                    {
                        auto& B1c       = B1(component);
                        auto const& B0c = B0(component);
                        layout.evalOnGhostBox(
                            B1c, [&](auto&... args) mutable { B1c(args...) -= B0c(args...); });
                    }
            }

            // The B1 potential init used E as the A scratch buffer; clear E so t=0 diagnostics
            // and the first read see 0 (constrained transport recomputes E before its real use).
            if (b1FromPotential_)
                for (auto const& component : {Component::X, Component::Y, Component::Z})
                {
                    auto& Ec = E(component);
                    layout.evalOnGhostBox(Ec, [&](auto&... args) mutable { Ec(args...) = 0.0; });
                }

            ToConservativeConverter{layout, gamma_}(
                rho, V, B1, P, rhoV, Etot1); // initial to conservative conversion because we
                                             // store conservative quantities on the grid
        }

        field_type rho;
        VecFieldT V;
        VecFieldT B1; // the evolved perturbation field (the real conserved magnetic variable)
        VecFieldT B;  // total-field scratch B = B0 + B1, recomputed on demand (tagging); NOT evolved
        field_type P;

        VecFieldT rhoV;
        field_type Etot1;

        VecFieldT E;
        VecFieldT J;

    private:
        initializer::InitFunction<dimension> rhoinit_;
        VecFieldInitializer<dimension> Vinit_;
        VecFieldInitializer<dimension> totalBInit_;
        initializer::InitFunction<dimension> Pinit_;

        double const gamma_;

        // Vector-potential init: build B1 = curl(A1) from a full 3D vector potential so that
        // div B = 0 discretely. b0FromPotential_ only governs the B1 = total - B0 subtraction.
        // Defaults keep the legacy component-wise init.
        bool b0FromPotential_ = false;
        bool b1FromPotential_ = false;
        VecFieldInitializer<dimension> a1Init_;
    };
} // namespace core
} // namespace PHARE

#endif // PHARE_MHD_STATE_HPP
