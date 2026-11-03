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
#include <vector>
#include <utility>

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

    NO_DISCARD VecField const& velocity() const { return Ve_; }
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

    // a closure refers to the flux computer of the model that owns it: a copy of the
    // closure is made for another model, see clone
    ElectronPressureClosure(ElectronPressureClosure const& that, FluxComputer const& flux)
        : flux_{flux}
        , Pe_{that.Pe_}
        , resources{that.resources}
    {
    }

    ElectronPressureClosure(ElectronPressureClosure const&)            = delete;
    ElectronPressureClosure& operator=(ElectronPressureClosure const&) = delete;

    // a copy of this closure, bound to the given flux computer
    NO_DISCARD virtual std::unique_ptr<ElectronPressureClosure> clone(FluxComputer const& flux) const
        = 0;

    virtual ~ElectronPressureClosure() {}

    //-------------------------------------------------------------------------
    //                  start the ResourcesUser interface
    //-------------------------------------------------------------------------

    // the flux computer is a resource of the model, not of its closure; the fields a
    // closure keeps for itself (resources) are, like Pe_
    NO_DISCARD virtual bool isUsable() const
    {
        return core::isUsable(Pe_) and std::all_of(std::begin(resources), std::end(resources),
                                                   [](auto const& var) {
                                                       return std::visit(
                                                           [](auto const& r) {
                                                               return r.isUsable();
                                                           },
                                                           var);
                                                   });
    }

    NO_DISCARD virtual bool isSettable() const
    {
        return core::isSettable(Pe_)
               and std::all_of(std::begin(resources), std::end(resources), [](auto const& var) {
                      return std::visit([](auto const& r) { return r.isSettable(); }, var);
                  });
    }

    struct PressureProperty // needed?
    {
        std::string name;
        typename HybridQuantity::Scalar qty;
    };

    using PressureProperties = std::vector<PressureProperty>;


    NO_DISCARD auto getCompileTimeResourcesViewList() const { return std::forward_as_tuple(Pe_); }

    NO_DISCARD auto getCompileTimeResourcesViewList() { return std::forward_as_tuple(Pe_); }

    //-------------------------------------------------------------------------
    //                  ends the ResourcesUser interface
    //-------------------------------------------------------------------------


    void virtual initialize(GridLayout const& layout) = 0;

    NO_DISCARD Field& pressure() { return Pe_; }

    NO_DISCARD Field const& pressure() const { return Pe_; }

    // sets the pressure at the initial time, once the ion density is deposited on the
    // layout's domain nodes; for closures whose pressure depends on their own state
    // and the density (see the polytropic one), the initial state is the user's
    void virtual seedPressure(GridLayout const& /*layout*/) {}

    // sets the pressure from the state saved by beginStep, which is then set by the
    // caller, on a level that has just been regridded, once its density is deposited
    void virtual pressureFromEntropy(GridLayout const& /*layout*/) {}

    // called once per time step, before the first computePressure of that step,
    // for closures that evolve the pressure in time (stages of a step all advance
    // from the state saved here, like Faraday advances from B at time n)
    void virtual beginStep(GridLayout const& /*layout*/) {}

    // name of the field saved by beginStep whose ghost nodes computePressure reads,
    // if any: the model fills them between the two calls
    NO_DISCARD virtual std::optional<std::string> stepStateName() const { return std::nullopt; }

    // sub-cycling of a closure that transports its step state: each stage of a step
    // calls startSubcycles, then subcycle and a ghost fill of the field named by
    // subcycleStateName (nbrSubcycles - 1) times, then computePressure, all over
    // dt / nbrSubcycles. nbrSubcycles is the max of subcycles over the level, so that
    // every patch takes the same sub-steps between the collective ghost fills.
    NO_DISCARD virtual std::size_t subcycles(GridLayout const& /*layout*/,
                                             double const /*dt*/) const
    {
        return 1;
    }
    NO_DISCARD virtual std::optional<std::string> subcycleStateName() const { return std::nullopt; }
    void virtual startSubcycles(GridLayout const& /*layout*/) {}
    void virtual subcycle(GridLayout const& /*layout*/, double const /*dt*/) {}

    void virtual computePressure(GridLayout const& /*layout*/, double const /* dt */) = 0;


    NO_DISCARD auto& getRunTimeResourcesViewList() { return resources; }
    NO_DISCARD auto& getRunTimeResourcesViewList() const { return resources; }

protected:
    FluxComputer const& flux_;
    Field Pe_;

    using Resources = std::variant<Field>;
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
        , Te_{dict["pressure_closure"]["Te"]}
    {
    }

    IsothermalElectronPressureClosure(IsothermalElectronPressureClosure const& that,
                                      FluxComputer const& flux)
        : Super{that, flux}
        , Te_{that.Te_}
    {
    }

    NO_DISCARD std::unique_ptr<Super> clone(FluxComputer const& flux) const override
    {
        return std::make_unique<IsothermalElectronPressureClosure>(*this, flux);
    }

    void initialize(GridLayout const& layout) override {}

    // the pressure only depends on the density: it is set on a new level, or a regridded
    // one, as soon as its density is known
    void seedPressure(GridLayout const& layout) override { computePressure(layout, 0.0); }
    void pressureFromEntropy(GridLayout const& layout) override { computePressure(layout, 0.0); }

    void computePressure(GridLayout const& /*layout*/, double const /* dt */) override
    {
        static_assert(Field::is_contiguous, "Error - assumes Field date is contiguous");

        if (!this->Pe_.isUsable())
            throw std::runtime_error("Error - isothermal closure is not usable");

        auto const& Ne_ = this->flux_.density();
        std::transform(std::begin(Ne_), std::end(Ne_), std::begin(this->Pe_),
                       [this](auto n) { return n * Te_; });
    }

private:
    double const Te_ = 0;
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
        , gamma_{dict["pressure_closure"]["gamma"]}
        , Pe_init_{dict["pressure_closure"]["Pe"].template to<initializer::InitFunction<dim>>()}
        , s_init_{dict["pressure_closure"]["s"].template to<initializer::InitFunction<dim>>()}
    {
        resources.emplace_back(sOld_);
        resources.emplace_back(s_);
    }

    PolytropicElectronPressureClosure(PolytropicElectronPressureClosure const& that,
                                      FluxComputer const& flux)
        : Super{that, flux}
        , gamma_{that.gamma_}
        , Pe_init_{that.Pe_init_}
        , s_init_{that.s_init_}
        , sOld_{that.sOld_}
        , s_{that.s_}
    {
    }

    NO_DISCARD std::unique_ptr<Super> clone(FluxComputer const& flux) const override
    {
        return std::make_unique<PolytropicElectronPressureClosure>(*this, flux);
    }

    void initialize(GridLayout const& layout) override
    {
        FieldUserFunctionInitializer::initialize(this->Pe_, layout, Pe_init_);
    }

    // the initial s is the user's analytic entropy Pe(x)/n(x)^gamma, see pharein, not
    // Pe_ over the deposited density: that one carries the particle noise of the
    // density, which s would then keep for the whole run since it is only transported.
    // So Pe_ = s*N^gamma on the domain, N being the deposited density, and beginStep
    // gets s back from it (same N) up to round off.
    void seedPressure(GridLayout const& layout) override
    {
        FieldUserFunctionInitializer::initialize(get_from_variants(resources, sOld_), layout,
                                                 s_init_);
        pressureFromEntropy(layout);
    }

    // where the density is not positive, Pe_ is left as it is
    void pressureFromEntropy(GridLayout const& layout) override
    {
        auto const& sOld = get_from_variants(resources, sOld_);
        auto const& N_   = this->flux_.density();

        layout.evalOnBox(this->Pe_, [&](auto const&... ijk) {
            auto const n = N_(ijk...);
            if (n > 0)
                this->Pe_(ijk...) = sOld(ijk...) * std::pow(n, gamma_);
        });
    }

    // Pe_ and the density are both at time n here, so is sOld_. Domain only: on
    // ghost nodes the density can be incomplete (next to a level border, level
    // ghost particles only cover part of the stencil), so Pe/N^gamma would be
    // wrong there. The ghosts of sOld_ are filled by the model afterwards, from
    // neighbor patches or interpolated from the coarser level's sOld_, see
    // stepStateName(). The initial density is checked against min_initial_density in
    // pyphare. If the deposited density still drops to 0 at a node during the run
    // (particles leaving it), s and Pe are left as they were there: the node keeps
    // its entropy for when it refills, and neighbors do not upwind from a zero.
    // Vacuum is not supported otherwise (Ve = (sum_s q_s F_s - J)/n and Ohm diverge there).
    void beginStep(GridLayout const& layout) override
    {
        if (!this->Pe_.isUsable())
            throw std::runtime_error("Error - polytropic closure not usable");

        auto const& N_ = this->flux_.density();
        auto& sOld     = get_from_variants(resources, sOld_);

        layout.evalOnBox(sOld, [&](auto const&... ijk) {
            auto const n = N_(ijk...);
            if (n > 0)
                sOld(ijk...) = this->Pe_(ijk...) / std::pow(n, gamma_);
        });
    }

    NO_DISCARD std::optional<std::string> stepStateName() const override { return sOld_.name(); }

    // the polytropic law dPe/dt = -(Ve.grad)Pe - gamma*Pe*div(Ve) together with the
    // electron continuity (ne = ni, div J = 0) is equivalent to the transport of
    // s = Pe/N^gamma: ds/dt + (Ve.grad)s = 0, and then Pe = s*N^gamma. This avoids
    // div(Ve), which contains div(J/N) and is noisy, and the compression is carried
    // by the ion density, consistent with the particles. With a uniform initial s
    // (Pe = P0*(N/N0)^gamma) the transport is inert and this reduces to the
    // algebraic polytropic law; a non uniform s is carried by the electron flow.
    //
    // every stage of a step advances s from sOld_ (time n) over the full dt with the
    // stage's velocity and density, so calling it once per predictor/corrector does
    // not accumulate: the corrector's result is the new state.
    //
    // the upwind update is stable only if the electron flow crosses at most one cell
    // per update (see courant_), and Ve contains J/n, which can exceed that in thin
    // or rarefied current layers at a time step fine for the ions and the fields. A
    // stage is then sub-cycled: s_ starts from sOld_ and takes ceil(C) sub-steps of
    // dt/ceil(C), its ghosts filled in between, the last one also setting Pe. Beyond
    // maxSubcycledCourant_ the time step is refused rather than sub-cycled further: such a
    // flow is most likely a nearly empty node, see min_initial_density.
    NO_DISCARD std::size_t subcycles(GridLayout const& layout, double const dt) const override
    {
        if (!this->Pe_.isUsable())
            throw std::runtime_error("Error - polytropic closure not usable");

        auto const [courant, worst] = patchMaxCourant_(layout, dt);
        if (courant > maxSubcycledCourant_)
        {
            auto const& N_ = this->flux_.density();
            auto const& V_ = this->flux_.velocity();

            std::stringstream msg;
            msg << "Error - polytropic closure: sum over directions of |Ve|*dt/dx is " << courant
                << " > " << maxSubcycledCourant_ << " on level " << layout.levelNumber()
                << " at AMR node " << layout.localToAMR(Point{worst}) << " (n = " << N_(worst)
                << ", Ve =";
            for (std::size_t d = 0; d < dim; ++d)
                msg << " " << V_[d](worst);
            msg << "), too many sub-steps for the electron pressure advection, reduce the time "
                   "step";
            throw std::runtime_error(msg.str());
        }

        // the tolerance keeps an integer Courant number, up to round off, from taking
        // one more sub-step than needed
        return std::max<std::size_t>(1, std::ceil(courant * (1 - courantTolerance_)));
    }

    NO_DISCARD std::optional<std::string> subcycleStateName() const override { return s_.name(); }

    // ghost nodes included: they hold sOld_'s neighbor and coarser level values at the start
    // of the step, which the sub-steps' ghost fill refreshes at the time of each sub-step,
    // from neighbor patches and time interpolated from the coarser level.
    void startSubcycles(GridLayout const& /*layout*/) override
    {
        auto const& sOld = get_from_variants(resources, sOld_);
        auto& s          = get_from_variants(resources, s_);
        std::copy(std::begin(sOld), std::end(sOld), std::begin(s));
    }

    void subcycle(GridLayout const& layout, double const dt) override
    {
        auto& s       = get_from_variants(resources, s_);
        auto const s1 = advance_(layout, dt);

        std::size_t i = 0;
        layout.evalOnBox(s, [&](auto const&... ijk) { s(ijk...) = s1[i++]; });
    }

    // the last sub-step (the whole stage without sub-cycling) of s_, then
    // Pe = s*N^gamma
    void computePressure(GridLayout const& layout, double const dt) override
    {
        auto const& N_ = this->flux_.density();
        auto& sWork    = get_from_variants(resources, s_);
        auto const s1  = advance_(layout, dt);

        // domain only: Ohm reads Pe_ on domain nodes only, its ghosts are not used.
        // s_ is left at the end of the stage, i.e. time n+1 after the last one: finer
        // levels time interpolate their ghosts between sOld_ (time n) and s_
        std::size_t i = 0;
        layout.evalOnBox(this->Pe_, [&](auto const&... ijk) {
            auto const index = std::array{ijk...};
            auto const s     = s1[i++];
            sWork(index)     = s;
            auto const n     = N_(index);
            if (n <= 0) // see beginStep
                return;

            this->Pe_(index) = s * std::pow(n, gamma_);
        });
    }

private:
    // most sub-steps per stage, see subcycles
    static constexpr double maxSubcycledCourant_ = 10;
    static constexpr double courantTolerance_    = 1e-12;

    // max over the patch domain of courant_, and where it is reached
    auto patchMaxCourant_(GridLayout const& layout, double const dt) const
    {
        auto const& N_             = this->flux_.density();
        auto const& V_             = this->flux_.velocity();
        auto const inverseMeshSize = layout.inverseMeshSize();

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
        return std::make_pair(maxCourant, worst);
    }

    // s_ advanced over dt on the domain nodes, in evalOnBox order, from s_ with its
    // ghosts filled.
    // Returned rather than written since the stencil reads s_'s neighbors. The
    // Courant number is checked before, so that a throw leaves s_ and Pe_ as they were
    std::vector<double> advance_(GridLayout const& layout, double const dt) const
    {
        if (!this->Pe_.isUsable())
            throw std::runtime_error("Error - polytropic closure not usable");

        auto const& N_ = this->flux_.density();
        auto const& V_ = this->flux_.velocity();
        auto const& s  = get_from_variants(resources, s_);

        if (auto const courant = patchMaxCourant_(layout, dt).first;
            courant > 1 + courantTolerance_)
            throw std::runtime_error(
                "Error - polytropic closure: sum over directions of |Ve|*dt/dx is "
                + std::to_string(courant) + " > 1 for one update, see subcycles");

        auto const inverseMeshSize = layout.inverseMeshSize();

        std::vector<double> s1;
        layout.evalOnBox(s, [&](auto const&... ijk) {
            auto const index = std::array{ijk...};
            if (N_(index) <= 0) // see beginStep
                s1.push_back(s(index));
            else
                s1.push_back(s(index) - dt * upwindAdvection_(V_, s, index, inverseMeshSize));
        });
        return s1;
    }

    double const gamma_ = 5. / 3.;
    initializer::InitFunction<dim> Pe_init_;
    initializer::InitFunction<dim> s_init_;

    // s = Pe/N^gamma at the start of the current time step, set by beginStep. Must
    // persist across the stages of a step, so it cannot share a scratch resource
    // like "PHARE_sumField", which the messenger reuses while moving ions.
    Field sOld_{"Se_old", HybridQuantity::Scalar::P};

    // s within a stage, advanced from sOld_ by the sub-steps, see subcycles
    Field s_{"Se", HybridQuantity::Scalar::P};


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
std::unique_ptr<ElectronPressureClosure<FluxComputer>>
ElectronPressureClosureFactory(PHARE::initializer::PHAREDict const& dict, FluxComputer const& flux)
{
    auto const name = dict["pressure_closure"]["name"].template to<std::string>();

    if (name == "isothermal")
        return std::make_unique<IsothermalElectronPressureClosure<FluxComputer>>(dict, flux);
    if (name == "polytropic")
        return std::make_unique<PolytropicElectronPressureClosure<FluxComputer>>(dict, flux);

    throw std::runtime_error("Error - unknown electron pressure closure \"" + name + "\"");
}



template<typename FluxComputer>
class ElectronMomentModel
{
    using GridLayout = FluxComputer::GridLayout;
    using VecField   = FluxComputer::VecField;
    using Field      = FluxComputer::Field;

public:
    ElectronMomentModel(PHARE::initializer::PHAREDict const& dict, FluxComputer& flux)
        : dict_{dict}
        , fluxComput_{flux}
        , pressureClosure_{ElectronPressureClosureFactory<FluxComputer>(dict, fluxComput_)}
    {
    }

    ElectronMomentModel(ElectronMomentModel const& self)
        : dict_{self.dict_}
        , fluxComput_{self.fluxComput_}
        , pressureClosure_{self.pressureClosure_->clone(fluxComput_)}
    {
    }

    //-------------------------------------------------------------------------
    //                  start the ResourcesUser interface
    //-------------------------------------------------------------------------

    NO_DISCARD bool isUsable() const
    {
        return fluxComput_.isUsable() and pressureClosure_->isUsable();
    }

    NO_DISCARD bool isSettable() const
    {
        return fluxComput_.isSettable() and pressureClosure_->isSettable();
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
    void seedPressure(GridLayout const& layout) { pressureClosure_->seedPressure(layout); }
    void pressureFromEntropy(GridLayout const& layout)
    {
        pressureClosure_->pressureFromEntropy(layout);
    }
    void beginStep(GridLayout const& layout) { pressureClosure_->beginStep(layout); }
    NO_DISCARD auto stepStateName() const { return pressureClosure_->stepStateName(); }
    NO_DISCARD auto subcycles(GridLayout const& layout, double const dt) const
    {
        return pressureClosure_->subcycles(layout, dt);
    }
    NO_DISCARD auto subcycleStateName() const { return pressureClosure_->subcycleStateName(); }
    void startSubcycles(GridLayout const& layout) { pressureClosure_->startSubcycles(layout); }
    void subcycle(GridLayout const& layout, double const dt)
    {
        pressureClosure_->subcycle(layout, dt);
    }
    void computePressure(GridLayout const& layout, double const dt)
    {
        pressureClosure_->computePressure(layout, dt);
    }

private:
    initializer::PHAREDict dict_;
    FluxComputer fluxComput_;
    std::unique_ptr<ElectronPressureClosure<FluxComputer>> pressureClosure_;
};



template<typename FluxComputer>
class Electrons
{
    using GridLayout = FluxComputer::GridLayout;
    using VecField   = FluxComputer::VecField;
    using Field      = FluxComputer::Field;

public:
    Electrons(initializer::PHAREDict const& dict, FluxComputer flux)
        : dict_{dict}
        , momentModel_{dict, flux}
    {
    }

    Electrons(Electrons const& that) = default;

    void initialize(GridLayout const& layout) { momentModel_.initialize(layout); }

    // exposed separately so callers needing a level-wide step between
    // sub-steps (e.g. ghost filling) can drive them instead of update() below.
    void computeDensity() { momentModel_.computeDensity(); }
    void computeBulkVelocity(GridLayout const& layout) { momentModel_.computeBulkVelocity(layout); }
    void seedPressure(GridLayout const& layout) { momentModel_.seedPressure(layout); }
    void pressureFromEntropy(GridLayout const& layout) { momentModel_.pressureFromEntropy(layout); }
    void beginStep(GridLayout const& layout) { momentModel_.beginStep(layout); }
    NO_DISCARD auto stepStateName() const { return momentModel_.stepStateName(); }
    NO_DISCARD auto subcycles(GridLayout const& layout, double const dt) const
    {
        return momentModel_.subcycles(layout, dt);
    }
    NO_DISCARD auto subcycleStateName() const { return momentModel_.subcycleStateName(); }
    void startSubcycles(GridLayout const& layout) { momentModel_.startSubcycles(layout); }
    void subcycle(GridLayout const& layout, double const dt) { momentModel_.subcycle(layout, dt); }
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
