#ifndef PHARE_ELECTRONS_HPP
#define PHARE_ELECTRONS_HPP

#include "core/def.hpp"
#include "core/models/quantities/hybrid_quantities.hpp"
#include "core/utilities/variants.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/data/vecfield/vecfield_component.hpp"
#include "core/data/field/initializers/field_user_initializer.hpp"

#include "initializer/data_provider.hpp"



#include <array>
#include <cmath>
#include <memory>
#include <optional>
#include <algorithm>
#include <sstream>
#include <stdexcept>
#include <variant>

namespace PHARE::core
{

template<typename Ions>
class StandardHybridElectronFluxComputer
{
public:
    using VecField   = Ions::vecfield_type;
    using Field      = Ions::field_type;
    using GridLayout = Ions::gridlayout_type;

    StandardHybridElectronFluxComputer(Ions& ions, VecField& J)
        : ions_{ions}
        , J_{J}
        , Ve_{"StandardHybridElectronFluxComputer_Ve", HybridQuantity::Vector::V}
    {
    }

    //-------------------------------------------------------------------------
    //                  start the ResourcesUser interface
    //-------------------------------------------------------------------------

    NO_DISCARD bool isUsable() const { return core::isUsable(ions_, J_, Ve_); }

    NO_DISCARD bool isSettable() const { return core::isSettable(ions_, J_, Ve_); }

    NO_DISCARD auto getCompileTimeResourcesViewList() const
    {
        return std::forward_as_tuple(Ve_, ions_, J_);
    }

    NO_DISCARD auto getCompileTimeResourcesViewList()
    {
        return std::forward_as_tuple(Ve_, ions_, J_);
    }

    //-------------------------------------------------------------------------
    //                  ends the ResourcesUser interface
    //-------------------------------------------------------------------------


    NO_DISCARD Field const& density() const { return ions_.chargeDensity(); }

    NO_DISCARD Field& density() { return ions_.chargeDensity(); }

    NO_DISCARD VecField& velocity() { return Ve_; }

    void computeDensity() {}

    void computeBulkVelocity(GridLayout const& layout)
    {
        auto const& Jx = J_(Component::X);
        auto const& Jy = J_(Component::Y);
        auto const& Jz = J_(Component::Z);
        // Ve = (sum_s q_s F_s - J) / Ne, with Ne the total ion charge density.
        // The ion bulk velocity is mass weighted and cannot be used here, see #1466
        auto const& Fx = ions_.chargeFlux()(Component::X);
        auto const& Fy = ions_.chargeFlux()(Component::Y);
        auto const& Fz = ions_.chargeFlux()(Component::Z);
        auto const& Ne = ions_.chargeDensity();

        auto& Vex = Ve_(Component::X);
        auto& Vey = Ve_(Component::Y);
        auto& Vez = Ve_(Component::Z);

        // from Ni because all components are defined on primal
        layout.evalOnBox(Ne, [&](auto const&... args) {
            auto arr = std::array{args...};

            auto const JxOnVx = GridLayout::template project<GridLayout::JxToMoments>(Jx, arr);
            auto const JyOnVy = GridLayout::template project<GridLayout::JyToMoments>(Jy, arr);
            auto const JzOnVz = GridLayout::template project<GridLayout::JzToMoments>(Jz, arr);

            Vex(arr) = (Fx(arr) - JxOnVx) / Ne(arr);
            Vey(arr) = (Fy(arr) - JyOnVy) / Ne(arr);
            Vez(arr) = (Fz(arr) - JzOnVz) / Ne(arr);
        });
    }

    auto& getIons() const { return ions_; }

private:
    Ions ions_;
    VecField J_;
    VecField Ve_;
};



template<typename FluxComputer>
class ElectronPressureClosure
{
    using GridLayout = FluxComputer::GridLayout;
    using VecField   = FluxComputer::VecField;
    using Field      = FluxComputer::Field;

public:
    ElectronPressureClosure(PHARE::initializer::PHAREDict const& dict, FluxComputer const& flux)
        : flux_{flux}
        , Pe_{"Pe", HybridQuantity::Scalar::P}
    {
    }

    virtual ~ElectronPressureClosure() {}

    //-------------------------------------------------------------------------
    //                  start the ResourcesUser interface
    //-------------------------------------------------------------------------

    NO_DISCARD virtual bool isUsable() const { return core::isUsable(Pe_, flux_); }

    NO_DISCARD virtual bool isSettable() const { return core::isSettable(Pe_, flux_); }

    struct PressureProperty // needed?
    {
        std::string name;
        typename HybridQuantity::Scalar qty;
    };

    using PressureProperties = std::vector<PressureProperty>;


    NO_DISCARD auto getCompileTimeResourcesViewList() const
    {
        return std::forward_as_tuple(flux_, Pe_);
    }

    NO_DISCARD auto getCompileTimeResourcesViewList() { return std::forward_as_tuple(flux_, Pe_); }

    //-------------------------------------------------------------------------
    //                  ends the ResourcesUser interface
    //-------------------------------------------------------------------------


    void virtual initialize(GridLayout const& layout) = 0;

    NO_DISCARD Field& pressure() { return Pe_; }

    NO_DISCARD Field const& pressure() const { return Pe_; }

    // called once per time step, before the first computePressure of that step,
    // for closures that evolve the pressure in time (stages of a step all advance
    // from the state saved here, like Faraday advances from B at time n)
    void virtual beginStep(GridLayout const& /*layout*/) {}

    // name of the field saved by beginStep whose ghost nodes computePressure reads,
    // if any: the model fills them between the two calls
    NO_DISCARD virtual std::optional<std::string> stepStateName() const { return std::nullopt; }

    void virtual computePressure(GridLayout const& /*layout*/, double const /* dt */) = 0;


    NO_DISCARD auto& getRunTimeResourcesViewList() { return resources; }
    NO_DISCARD auto& getRunTimeResourcesViewList() const { return resources; }

    NO_DISCARD auto& B() const { return get_as_ref_or_throw<VecField const>(resources); }
    NO_DISCARD auto& B() { return get_as_ref_or_throw<VecField>(resources); }

    NO_DISCARD auto& Te() const
    {
        auto const& field = get_as_ref_or_throw<Field const>(resources);
        assert(field.name() == "Te");
        return field;
    }

protected:
    FluxComputer flux_;
    Field Pe_;

    using Resources = std::variant<VecField, Field>;
    std::vector<Resources> resources;
};



template<typename FluxComputer>
class IsothermalElectronPressureClosure : public ElectronPressureClosure<FluxComputer>
{
    using GridLayout = FluxComputer::GridLayout;
    using VecField   = FluxComputer::VecField;
    using Field      = FluxComputer::Field;

    using Super = ElectronPressureClosure<FluxComputer>;
    using Super::getCompileTimeResourcesViewList;
    using Super::isSettable;
    using Super::isUsable;

public:
    using field_type = Field;

    IsothermalElectronPressureClosure(PHARE::initializer::PHAREDict const& dict,
                                      FluxComputer const& flux)
        : Super{dict, flux}
        , T0_{dict["pressure_closure"]["Te"]}
    {
    }

    void initialize(GridLayout const& layout) override {}

    void computePressure(GridLayout const& /*layout*/, double const /* dt */) override
    {
        static_assert(Field::is_contiguous, "Error - assumes Field date is contiguous");

        if (!this->Pe_.isUsable())
            throw std::runtime_error("Error - isothermal closure is not usable");

        auto const& Ne_ = this->flux_.density();
        std::transform(std::begin(Ne_), std::end(Ne_), std::begin(this->Pe_),
                       [this](auto n) { return n * T0_; });
    }

private:
    double const T0_ = 0;
};



template<typename FluxComputer>
class PolytropicElectronPressureClosure : public ElectronPressureClosure<FluxComputer>
{
    using GridLayout = FluxComputer::GridLayout;
    using VecField   = FluxComputer::VecField;
    using Field      = FluxComputer::Field;

    using Super = ElectronPressureClosure<FluxComputer>;
    using Super::getCompileTimeResourcesViewList;
    using Super::isSettable;
    using Super::isUsable;

    static constexpr std::size_t dim = VecField::dimension;


public:
    using field_type = Field;
    using Super::resources;

    PolytropicElectronPressureClosure(PHARE::initializer::PHAREDict const& dict,
                                      FluxComputer const& flux)
        : Super{dict, flux}
        , gamma_{dict["pressure_closure"]["Gamma"]}
        , Pe_init_{dict["pressure_closure"]["Pe"].template to<initializer::InitFunction<dim>>()}
    {
        resources.emplace_back(Te_);
        resources.emplace_back(sN_);
    }

    void initialize(GridLayout const& layout) override
    {
        FieldUserFunctionInitializer::initialize(this->Pe_, layout, Pe_init_);
    }

    // Pe_ and the density are both at time n here, so is sN_. Domain only: on
    // ghost nodes the density can be incomplete (next to a level border, level
    // ghost particles only cover part of the stencil), so Pe/N^gamma would be
    // wrong there. The ghosts of sN_ are filled by the model afterwards, from
    // neighbor patches or interpolated from the coarser level's sN_, see
    // stepStateName(). The initial density is checked against min_density in
    // pyphare. If the deposited density still drops to 0 at a node during the run
    // (particles leaving it), s, Pe and Te are left as they were there: the node keeps
    // its entropy for when it refills, and neighbors do not upwind from a zero.
    // Vacuum is not supported otherwise (Ve = (sum_s q_s F_s - J)/n and Ohm diverge there).
    void beginStep(GridLayout const& layout) override
    {
        if (!this->Pe_.isUsable())
            throw std::runtime_error("Error - polytropic closure not usable");

        auto const& N_ = this->flux_.density();
        auto& sN       = get_from_variants(resources, sN_);

        layout.evalOnBox(sN, [&](auto const&... ijk) {
            auto const n = N_(ijk...);
            if (n > 0)
                sN(ijk...) = this->Pe_(ijk...) / std::pow(n, gamma_);
        });
    }

    NO_DISCARD std::optional<std::string> stepStateName() const override { return sN_.name(); }

    // the polytropic law dPe/dt = -(Ve.grad)Pe - gamma*Pe*div(Ve) together with the
    // electron continuity (ne = ni, div J = 0) is equivalent to the transport of
    // s = Pe/N^gamma: ds/dt + (Ve.grad)s = 0, and then Pe = s*N^gamma. This avoids
    // div(Ve), which contains div(J/N) and is noisy, and the compression is carried
    // by the ion density, consistent with the particles. With a uniform initial s
    // (Pe = P0*(N/N0)^gamma) the transport is inert and this reduces to the
    // algebraic polytropic law; a non uniform s is carried by the electron flow.
    //
    // every stage of a step advances s from sN_ (time n) over the full dt with the
    // stage's velocity and density, so calling it once per predictor/corrector does
    // not accumulate: the corrector's result is the new state.
    void computePressure(GridLayout const& layout, double const dt) override
    {
        if (!this->Pe_.isUsable())
            throw std::runtime_error("Error - polytropic closure not usable");

        auto const& N_ = this->flux_.density();
        auto const& V_ = this->flux_.velocity();

        auto& Te       = get_from_variants(resources, Te_);
        auto const& sN = get_from_variants(resources, sN_);

        auto const inverseMeshSize = layout.inverseMeshSize();

        // checked on the whole patch before anything is written, so that a throw
        // leaves Pe_ as it was
        double maxCourant = 0;
        std::array<std::uint32_t, dim> worst{};
        layout.evalOnBox(this->Pe_, [&](auto const&... ijk) {
            auto const index = std::array{ijk...};
            if (N_(index) <= 0)
                return;
            if (auto const c = courant_(V_, index, inverseMeshSize, dt); c > maxCourant)
            {
                maxCourant = c;
                worst      = index;
            }
        });
        if (maxCourant > 1)
        {
            std::stringstream msg;
            msg << "Error - polytropic closure: sum over directions of |Ve|*dt/dx is " << maxCourant
                << " > 1 on level " << layout.levelNumber() << " at AMR node "
                << layout.localToAMR(Point{worst}) << " (n = " << N_(worst) << ", Ve =";
            for (std::size_t d = 0; d < dim; ++d)
                msg << " " << V_[d](worst);
            msg << "), the electron pressure advection is unstable, reduce the time step";
            throw std::runtime_error(msg.str());
        }

        // domain only: Ohm reads Pe_ on domain nodes only, its ghosts are not used
        layout.evalOnBox(this->Pe_, [&](auto const&... ijk) {
            auto const index = std::array{ijk...};
            auto const n     = N_(index);
            if (n <= 0) // see beginStep
                return;

            auto const s     = sN(index) - dt * upwindAdvection_(V_, sN, index, inverseMeshSize);
            this->Pe_(index) = s * std::pow(n, gamma_);
            Te(index)        = this->Pe_(index) / n;
        });
    }

private:
    double const gamma_ = 5. / 3.;
    initializer::InitFunction<dim> Pe_init_;
    Field Te_{"Te", HybridQuantity::Scalar::P};

    // s = Pe/N^gamma at the start of the current time step, set by beginStep. Must
    // persist across the stages of a step, so it cannot share a scratch resource
    // like "PHARE_sumField", which the messenger reuses while moving ions.
    Field sN_{"Se_n", HybridQuantity::Scalar::P};


    // the update below is unsplit: s_new = (1 - sum_d c_d) s_i + sum_d c_d s_upwind_d
    // with c_d = |Ve_d|*dt/dx_d. It is a convex combination, hence stable, monotone
    // and positivity preserving, only if sum_d c_d <= 1 (not max_d c_d <= 1).
    auto courant_(VecField const& Ve, auto const& index,
                  std::array<double, dim> const& inverseMeshSize, double const dt) const
    {
        double courant = 0;
        for (std::size_t d = 0; d < dim; ++d)
            courant += std::abs(Ve[d](index)) * dt * inverseMeshSize[d];
        return courant;
    }

    // (Ve.grad)s, each direction differenced on the side the flow comes from. Ve and
    // s are both primal on Yee, so the stencil only involves neighbor nodes.
    auto upwindAdvection_(VecField const& Ve, Field const& s, auto const& index,
                          std::array<double, dim> const& inverseMeshSize) const
    {
        double advection = 0;

        for (std::size_t d = 0; d < dim; ++d)
        {
            auto const v  = Ve[d](index);
            auto neighbor = index;
            if (v > 0)
            {
                --neighbor[d];
                advection += v * inverseMeshSize[d] * (s(index) - s(neighbor));
            }
            else
            {
                ++neighbor[d];
                advection += v * inverseMeshSize[d] * (s(neighbor) - s(index));
            }
        }

        return advection;
    }
};



template<typename FluxComputer>
class CGLElectronPressureClosure : public ElectronPressureClosure<FluxComputer>
{
    using GridLayout = FluxComputer::GridLayout;
    using VecField   = FluxComputer::VecField;
    using Field      = FluxComputer::Field;

    using Super = ElectronPressureClosure<FluxComputer>;
    using Super::getCompileTimeResourcesViewList;
    using Super::isSettable;
    using Super::isUsable;
    using Super::resources;

    static constexpr std::size_t dim = VecField::dimension;


public:
    using field_type = Field;

    CGLElectronPressureClosure(PHARE::initializer::PHAREDict const& dict, FluxComputer const& flux,
                               VecField const& B)
        : Super{dict, flux}
    {
        resources.emplace_back(B);
    }

    // CGL needs a gyrotropic pressure tensor (P_par, P_perp) and its divergence in
    // Ohm's law, which the scalar Pe_ interface cannot represent yet: only the
    // construction and B plumbing exist, using it must fail rather than run
    void initialize(GridLayout const& /*layout*/) override
    {
        throw std::runtime_error("Error - CGL electron pressure closure is not implemented yet");
    }

    void computePressure(GridLayout const& /*layout*/, double const /*dt*/) override
    {
        throw std::runtime_error("Error - CGL electron pressure closure is not implemented yet");
    }
};



template<typename FluxComputer>
std::unique_ptr<ElectronPressureClosure<FluxComputer>>
ElectronPressureClosureFactory(PHARE::initializer::PHAREDict const& dict, FluxComputer& flux,
                               typename FluxComputer::VecField const& B)
{
    if (dict["pressure_closure"]["name"].template to<std::string>() == "isothermal")
    {
        return std::make_unique<IsothermalElectronPressureClosure<FluxComputer>>(dict, flux);
    }
    else if (dict["pressure_closure"]["name"].template to<std::string>() == "polytropic")
    {
        return std::make_unique<PolytropicElectronPressureClosure<FluxComputer>>(dict, flux);
    }
    else if (dict["pressure_closure"]["name"].template to<std::string>() == "CGL")
    {
        return std::make_unique<CGLElectronPressureClosure<FluxComputer>>(dict, flux, B);
    }
    return nullptr;
}



template<typename FluxComputer>
class ElectronMomentModel
{
    using GridLayout = FluxComputer::GridLayout;
    using VecField   = FluxComputer::VecField;
    using Field      = FluxComputer::Field;

public:
    ElectronMomentModel(PHARE::initializer::PHAREDict const& dict, FluxComputer& flux,
                        VecField const& B)
        : dict_{dict}
        , fluxComput_{flux}
        , B_{B}
        , pressureClosure_{ElectronPressureClosureFactory<FluxComputer>(dict, flux, B)}
    {
        assert(pressureClosure_);
    }

    ElectronMomentModel(ElectronMomentModel const& self)
        : dict_{self.dict_}
        , fluxComput_{self.fluxComput_}
        , B_{self.B_}
        , pressureClosure_{ElectronPressureClosureFactory<FluxComputer>(dict_, fluxComput_, B_)}
    {
        *pressureClosure_ = *self.pressureClosure_;
    }

    //-------------------------------------------------------------------------
    //                  start the ResourcesUser interface
    //-------------------------------------------------------------------------

    NO_DISCARD bool isUsable() const
    {
        return fluxComput_.isUsable() /*and B_.isUsable()*/ and pressureClosure_->isUsable();
    }

    NO_DISCARD bool isSettable() const
    {
        return fluxComput_.isSettable() /*and B_.isSettable()*/ and pressureClosure_->isSettable();
    }

    NO_DISCARD auto getCompileTimeResourcesViewList() const
    {
        return std::forward_as_tuple(fluxComput_, *pressureClosure_);
    }

    NO_DISCARD auto getCompileTimeResourcesViewList()
    {
        return std::forward_as_tuple(fluxComput_, *pressureClosure_);
    }

    //-------------------------------------------------------------------------
    //                  ends the ResourcesUser interface
    //-------------------------------------------------------------------------

    void initialize(GridLayout const& layout) { pressureClosure_->initialize(layout); }

    NO_DISCARD Field const& density() const { return fluxComput_.density(); }
    NO_DISCARD VecField const& velocity() const { return fluxComput_.velocity(); }
    NO_DISCARD Field const& pressure() const { return pressureClosure_->pressure(); }

    NO_DISCARD Field& density() { return fluxComput_.density(); }
    NO_DISCARD VecField& velocity() { return fluxComput_.velocity(); }
    NO_DISCARD Field& pressure() { return pressureClosure_->pressure(); }

    void computeDensity() { fluxComput_.computeDensity(); }
    void computeBulkVelocity(GridLayout const& layout) { fluxComput_.computeBulkVelocity(layout); }
    void beginStep(GridLayout const& layout) { pressureClosure_->beginStep(layout); }
    NO_DISCARD auto stepStateName() const { return pressureClosure_->stepStateName(); }
    void computePressure(GridLayout const& layout, double const dt)
    {
        pressureClosure_->computePressure(layout, dt);
    }

private:
    initializer::PHAREDict dict_;
    FluxComputer fluxComput_;
    VecField B_; // needed for CGL
    std::unique_ptr<ElectronPressureClosure<FluxComputer>> pressureClosure_;
};



template<typename FluxComputer>
class Electrons
{
    using GridLayout = FluxComputer::GridLayout;
    using VecField   = FluxComputer::VecField;
    using Field      = FluxComputer::Field;

public:
    Electrons(initializer::PHAREDict const& dict, FluxComputer flux, VecField const& B)
        : dict_{dict}
        , momentModel_{dict, flux, B}
    {
    }

    Electrons(Electrons const& that) = default;

    void initialize(GridLayout const& layout) { momentModel_.initialize(layout); }

    // exposed separately so callers needing a level-wide step between
    // sub-steps (e.g. ghost filling) can drive them instead of update() below.
    void computeDensity() { momentModel_.computeDensity(); }
    void computeBulkVelocity(GridLayout const& layout) { momentModel_.computeBulkVelocity(layout); }
    void beginStep(GridLayout const& layout) { momentModel_.beginStep(layout); }
    NO_DISCARD auto stepStateName() const { return momentModel_.stepStateName(); }
    void computePressure(GridLayout const& layout, double const dt)
    {
        momentModel_.computePressure(layout, dt);
    }

    // convenience for closures without a step state (e.g. isothermal), used by unit
    // tests. A closure with a step state (polytropic) reads its ghost nodes, which
    // must be filled between beginStep and computePressure, even on a single patch:
    // use the split calls as SolverPPC::update_electrons does.
    void update(GridLayout const& layout, auto const dt)
    {
        if (stepStateName())
            throw std::runtime_error(
                "Error - Electrons::update cannot be used with a closure that has a step "
                "state ("
                + *stepStateName()
                + "): its ghosts must be filled between "
                  "beginStep and computePressure");

        if (isUsable())
        {
            computeDensity();
            computeBulkVelocity(layout);
            beginStep(layout);
            computePressure(layout, dt);
        }
        else
            throw std::runtime_error("Error - Electron  is not usable");
    }

    //-------------------------------------------------------------------------
    //                  start the ResourcesUser interface
    //-------------------------------------------------------------------------

    NO_DISCARD bool isUsable() const { return momentModel_.isUsable(); }

    NO_DISCARD bool isSettable() const { return momentModel_.isSettable(); }

    NO_DISCARD auto getCompileTimeResourcesViewList() const
    {
        return std::forward_as_tuple(momentModel_);
    }

    NO_DISCARD auto getCompileTimeResourcesViewList()
    {
        return std::forward_as_tuple(momentModel_);
    }

    //-------------------------------------------------------------------------
    //                  ends the ResourcesUser interface
    //-------------------------------------------------------------------------

    NO_DISCARD Field const& density() const { return momentModel_.density(); }
    NO_DISCARD VecField const& velocity() const { return momentModel_.velocity(); }
    NO_DISCARD Field const& pressure() const { return momentModel_.pressure(); }

    NO_DISCARD Field& density() { return momentModel_.density(); }
    NO_DISCARD VecField& velocity() { return momentModel_.velocity(); }
    NO_DISCARD Field& pressure() { return momentModel_.pressure(); }

private:
    initializer::PHAREDict dict_;
    ElectronMomentModel<FluxComputer> momentModel_;
};

} // namespace PHARE::core


#endif
