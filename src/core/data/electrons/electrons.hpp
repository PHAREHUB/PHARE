#ifndef PHARE_ELECTRONS_HPP
#define PHARE_ELECTRONS_HPP

#include "core/def.hpp"
#include "core/utilities/variants.hpp"
#include "core/hybrid/hybrid_quantities.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/data/vecfield/vecfield_component.hpp"
#include "core/data/field/initializers/field_user_initializer.hpp"
#include "core/hybrid/hybrid_quantities.hpp"
#include "core/data/vecfield/vecfield_component.hpp"

#include "initializer/data_provider.hpp"



#include <memory>
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


    NO_DISCARD Field const& density() const
    {
        if (!isUsable())
            throw std::runtime_error("Error, cannot return density because "
                                     "StandardHybridElectronFluxComputer is not usable");

        return ions_.chargeDensity();
    }

    NO_DISCARD Field& density()
    {
        if (!isUsable())
            throw std::runtime_error("Error, cannot return density because "
                                     "StandardHybridElectronFluxComputer is not usable");

        return ions_.chargeDensity();
    }

    NO_DISCARD VecField& velocity()
    {
        if (!isUsable())
            throw std::runtime_error("Error, cannot return velocity because "
                                     "StandardHybridElectronFluxComputer is not usable");

        return Ve_;
    }

    void computeDensity() {}

    void computeBulkVelocity(GridLayout const& layout)
    {
        auto const& Jx  = J_(Component::X);
        auto const& Jy  = J_(Component::Y);
        auto const& Jz  = J_(Component::Z);
        auto const& Vix = ions_.velocity()(Component::X);
        auto const& Viy = ions_.velocity()(Component::Y);
        auto const& Viz = ions_.velocity()(Component::Z);
        auto const& Ne  = ions_.chargeDensity();

        auto& Vex = Ve_(Component::X);
        auto& Vey = Ve_(Component::Y);
        auto& Vez = Ve_(Component::Z);

        // from Ni because all components are defined on primal
        layout.evalOnBox(Ne, [&](auto const&... args) {
            auto arr = std::array{args...};

            auto const JxOnVx = GridLayout::template project<GridLayout::JxToMoments>(Jx, arr);
            auto const JyOnVy = GridLayout::template project<GridLayout::JyToMoments>(Jy, arr);
            auto const JzOnVz = GridLayout::template project<GridLayout::JzToMoments>(Jz, arr);

            Vex(arr) = Vix(arr) - JxOnVx / Ne(arr);
            Vey(arr) = Viy(arr) - JyOnVy / Ne(arr);
            Vez(arr) = Viz(arr) - JzOnVz / Ne(arr);
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

    NO_DISCARD Field& pressure()
    {
        if (!Pe_.isUsable())
            throw std::runtime_error("Error - ! isothermal closure pressure not usable");
        return Pe_;
    }

    NO_DISCARD Field const& pressure() const
    {
        if (!Pe_.isUsable())
            throw std::runtime_error("Error - !! isothermal closure pressure not usable");
        return Pe_;
    }

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
        , T0_{dict["pressure_closure"]["Te"].template to<double>()}
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
        , gamma_{dict["pressure_closure"]["Gamma"].template to<double>()}
        , Pe_init_{dict["pressure_closure"]["Pe"].template to<initializer::InitFunction<dim>>()}
    {
        resources.emplace_back(Te_);
    }

    void initialize(GridLayout const& layout) override
    {
        FieldUserFunctionInitializer::initialize(this->Pe_, layout, Pe_init_);
    }

    void computePressure(GridLayout const& layout, double const dt) override
    {
        static_assert(Field::is_contiguous, "Error - assumes Field date is contiguous");

        if (!this->Pe_.isUsable())
            throw std::runtime_error("Error - polytropic closure not usable");

        auto const& N_ = this->flux_.density();
        auto const& V_ = this->flux_.velocity();

        auto& Te = get_from_variants(resources, Te_);

        std::transform(this->Pe_.begin(), this->Pe_.end(), N_.begin(), Te.begin(),
                       [](auto p, auto n) { return p / n; });

        this->dt_ = dt;

        layout.evalOnBox(Te, [&](auto&... ijk) mutable { P_Eq_(layout, V_, Te, ijk...); });


        // std::transform(std::begin(Ne_), std::end(Ne_), std::begin(this->Pe_),
        //                [this](auto n) { return n * 0.1; });  // TODO utiliser autre chose que
        //                transform et aller chercher Tnew
    }

private:
    double const gamma_ = 5. / 3.;
    double dt_;
    initializer::InitFunction<dim> Pe_init_;
    Field Te_{"Te", HybridQuantity::Scalar::P};


    template<typename Field, typename VecField>
    void P_Eq_(GridLayout const& layout, VecField const& Ve, Field const& Te, auto&... ijk) const
    {
        Te(ijk...)
            = Te(ijk...)
              + dt_ * (advection_(layout, Ve, Te, ijk...) + compression_(layout, Ve, Te, ijk...));
    }

    template<typename Field, typename VecField>
    auto advection_(GridLayout const& layout, VecField const& Ve, Field const& Te,
                    auto&... ijk) const
    {
        auto l_ = layout; // TODO so that it is non const

        // Both Ve and Te are primal on Yee, that is with the same centering
        // Hence, using 'derivOnSameCentering' guarantee that all the terms involved
        // are also on the same centering, that is primal on Yee
        auto const& Vx = Ve(Component::X);
        auto gradT_X   = l_.template derivOnSameCentering<Direction::X>(Te, {ijk...});

        if constexpr (dim == 1)
            return Vx(ijk...) * gradT_X;

        else
        {
            auto const& Vy = Ve(Component::Y);
            auto gradT_Y   = l_.template derivOnSameCentering<Direction::Y>(Te, {ijk...});

            if constexpr (dim == 2)
                return Vx(ijk...) * gradT_X + Vy(ijk...) * gradT_Y;
            else
            {
                auto const& Vz = Ve(Component::Z);
                auto gradT_Z   = l_.template derivOnSameCentering<Direction::Z>(Te, {ijk...});

                return Vx(ijk...) * gradT_X + Vy(ijk...) * gradT_Y + Vz(ijk...) * gradT_Z;
            }
        }
    }



    template<typename Field, typename VecField>
    auto compression_(GridLayout const& layout, VecField const& Ve, Field const& Te,
                      auto&... ijk) const
    {
        auto l_ = layout; // TODO so that it is non const

        // Both Ve and Te are primal on Yee, that is with the same centering
        // Hence, using 'derivOnSameCentering' guarantee that all the terms involved
        // are also on the same centering, that is primal on Yee
        auto const& Vx = Ve(Component::X);
        auto gradV_X   = l_.template derivOnSameCentering<Direction::X>(Vx, {ijk...});

        if constexpr (dim == 1)
            return (gamma_ - 1) * Te(ijk...) * gradV_X;

        else
        {
            auto const& Vy = Ve(Component::Y);
            auto gradV_Y   = l_.template derivOnSameCentering<Direction::Y>(Vx, {ijk...});

            if constexpr (dim == 2)
                return (gamma_ - 1) * Te(ijk...) * (gradV_X + gradV_Y);
            else
            {
                auto const& Vz = Ve(Component::Z);
                auto gradV_Z   = l_.template derivOnSameCentering<Direction::Z>(Vz, {ijk...});

                return (gamma_ - 1) * Te(ijk...) * (gradV_X + gradV_Y + gradV_Z);
            }
        }
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
        , gamma_{dict["pressure_closure"]["Gamma"].template to<double>()}
        , Pe_init_{dict["pressure_closure"]["Pe"].template to<initializer::InitFunction<dim>>()}
    {
        resources.emplace_back(B);
    }

    void initialize(GridLayout const& layout) override
    {
        FieldUserFunctionInitializer::initialize(this->Pe_, layout, Pe_init_);
    }

    void computePressure(GridLayout const& layout, double const dt) override
    {
        static_assert(Field::is_contiguous, "Error - assumes Field date is contiguous");

        if (!this->Pe_.isUsable())
            throw std::runtime_error("Error - CGL closure not usable");

        auto const& Ne_ = this->flux_.density();
        std::transform(std::begin(Ne_), std::end(Ne_), std::begin(this->Pe_),
                       [](auto n) { return n * 0.1; });
    }

private:
    double const gamma_ = 5. / 3.;
    initializer::InitFunction<dim> Pe_init_;
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
        // return fluxComput_.isUsable() and pressureClosure_->isUsable();
        return fluxComput_.isUsable() /*and B_.isUsable()*/ and pressureClosure_->isUsable();
    }

    // NO_DISCARD bool isSettable() const { return fluxComput_.isSettable(); }  // TODO
    // pressureClosure needs also to be settable ?
    NO_DISCARD bool isSettable() const
    {
        return fluxComput_.isSettable() /*and B_.isSettable()*/ and pressureClosure_->isSettable();
    }

    NO_DISCARD auto getCompileTimeResourcesViewList() const
    {
        // return std::forward_as_tuple(fluxComput_, *pressureClosure_);
        return std::forward_as_tuple(fluxComput_, *pressureClosure_);
    }

    NO_DISCARD auto getCompileTimeResourcesViewList()
    {
        // return std::forward_as_tuple(fluxComput_, *pressureClosure_);
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
    void computePressure(GridLayout const& layout, double const dt)
    {
        pressureClosure_->computePressure(layout, dt);
    }

private:
    initializer::PHAREDict dict_;
    FluxComputer fluxComput_;
    VecField B_; // NOT USED HERE!
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

    void update(GridLayout const& layout, auto const dt)
    {
        if (isUsable())
        {
            momentModel_.computeDensity();
            momentModel_.computeBulkVelocity(layout);
            momentModel_.computePressure(layout, dt);
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
