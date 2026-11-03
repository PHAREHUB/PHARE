#ifndef PHARE_ELECTRONS_HPP
#define PHARE_ELECTRONS_HPP

#include "core/def.hpp"
#include "core/models/quantities/hybrid_quantities.hpp"
#include "core/utilities/variants.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/data/vecfield/vecfield_component.hpp"

#include "initializer/data_provider.hpp"



#include <memory>
#include <algorithm>
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


    NO_DISCARD Field& pressure() { return Pe_; }

    NO_DISCARD Field const& pressure() const { return Pe_; }

    void virtual computePressure(GridLayout const& /*layout*/) = 0;


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

    void computePressure(GridLayout const& /*layout*/) override
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
std::unique_ptr<ElectronPressureClosure<FluxComputer>>
ElectronPressureClosureFactory(PHARE::initializer::PHAREDict const& dict, FluxComputer const& flux)
{
    auto const name = dict["pressure_closure"]["name"].template to<std::string>();

    if (name == "isothermal")
        return std::make_unique<IsothermalElectronPressureClosure<FluxComputer>>(dict, flux);

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

    NO_DISCARD Field const& density() const { return fluxComput_.density(); }
    NO_DISCARD VecField const& velocity() const { return fluxComput_.velocity(); }
    NO_DISCARD Field const& pressure() const { return pressureClosure_->pressure(); }

    NO_DISCARD Field& density() { return fluxComput_.density(); }
    NO_DISCARD VecField& velocity() { return fluxComput_.velocity(); }
    NO_DISCARD Field& pressure() { return pressureClosure_->pressure(); }

    void computeDensity() { fluxComput_.computeDensity(); }
    void computeBulkVelocity(GridLayout const& layout) { fluxComput_.computeBulkVelocity(layout); }
    void computePressure(GridLayout const& layout)
    {
        pressureClosure_->computePressure(layout);
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

    // exposed separately so callers needing a level-wide step between
    // sub-steps (e.g. ghost filling) can drive them instead of update() below.
    void computeDensity() { momentModel_.computeDensity(); }
    void computeBulkVelocity(GridLayout const& layout) { momentModel_.computeBulkVelocity(layout); }
    void computePressure(GridLayout const& layout)
    {
        momentModel_.computePressure(layout);
    }

    void update(GridLayout const& layout)
    {
        if (isUsable())
        {
            computeDensity();
            computeBulkVelocity(layout);
            computePressure(layout);
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
