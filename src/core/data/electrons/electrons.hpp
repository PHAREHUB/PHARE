#ifndef PHARE_ELECTRONS_HPP
#define PHARE_ELECTRONS_HPP

#include "core/def.hpp"
#include "core/models/quantities/hybrid_quantities.hpp"
#include "core/data/vecfield/vecfield_component.hpp"
#include "core/numerics/pressure_term/pressure_term.hpp"

#include "initializer/data_provider.hpp"


namespace PHARE::core
{


template<typename Ions>
class StandardHybridElectronFluxComputer
{
public:
    using VecField   = Ions::vecfield_type;
    using Field      = Ions::field_type;
    using GridLayout = Ions::gridlayout_type;
    using Scalar     = GridLayout::Scalar;

    StandardHybridElectronFluxComputer(Ions& ions, VecField& J)
        : ions_{ions}
        , J_{J}
        , Ve_{"StandardHybridElectronFluxComputer_Ve", HybridQuantity::Vector::V}
    {
    }

    //-------------------------------------------------------------------------
    //                  start the ResourcesUser interface
    //-------------------------------------------------------------------------


    NO_DISCARD bool isUsable() const { return ions_.isUsable() && J_.isUsable() && Ve_.isUsable(); }

    NO_DISCARD bool isSettable() const
    {
        return Ve_.isSettable() && ions_.isSettable() && J_.isSettable();
    }

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

        // from Ni because all components defined on primal
        layout.evalOnBox(Ne, [&](auto const&... args) {
            auto arr = std::array{args...};

            auto const JxOnVx = GridLayout::template project<Scalar::Jx, Scalar::Vx>(Jx, arr);
            auto const JyOnVy = GridLayout::template project<Scalar::Jy, Scalar::Vy>(Jy, arr);
            auto const JzOnVz = GridLayout::template project<Scalar::Jz, Scalar::Vz>(Jz, arr);

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



template<typename Ions>
class IsothermalElectronPressureClosure
{
    using GridLayout = Ions::gridlayout_type;
    using VecField   = Ions::vecfield_type;
    using Field      = Ions::field_type;

public:
    using field_type = Field;

    IsothermalElectronPressureClosure(PHARE::initializer::PHAREDict const& dict, Ions const& ions)
        : ions_{ions}
        , Te_{dict["Te"].template to<double>()}
        , Pe_{"Pe", HybridQuantity::Scalar::P}
        , pressureTerm_{"pressureTerm", HybridQuantity::Vector::E}
        , isothermalPressureTerm_{Te_}
    {
    }


    //-------------------------------------------------------------------------
    //                  start the ResourcesUser interface
    //-------------------------------------------------------------------------


    NO_DISCARD bool isUsable() const
    {
        return Pe_.isUsable() and pressureTerm_.isUsable() and ions_.isUsable();
    }

    NO_DISCARD bool isSettable() const { return Pe_.isSettable() and pressureTerm_.isSettable(); }

    struct PressureProperty
    {
        std::string name;
        HybridQuantity::Scalar qty;
    };

    using PressureProperties = std::vector<PressureProperty>;

    NO_DISCARD auto getCompileTimeResourcesViewList() const
    {
        return std::forward_as_tuple(ions_, Pe_, pressureTerm_);
    }

    NO_DISCARD auto getCompileTimeResourcesViewList()
    {
        return std::forward_as_tuple(ions_, Pe_, pressureTerm_);
    }

    //-------------------------------------------------------------------------
    //                  ends the ResourcesUser interface
    //-------------------------------------------------------------------------



    // Ohm does not read Pe: its pressure term comes from computePressureTerm, which only needs n
    // and Te here. Pe is kept for now for further internal use by non-isothermal pressure laws, and
    // for the upcoming electron pressure diagnostics (#1548).
    NO_DISCARD Field& pressure() { return Pe_; }
    NO_DISCARD Field const& pressure() const { return Pe_; }

    void computePressure(GridLayout const& /*layout*/)
    {
        static_assert(Field::is_contiguous, "Error - assumes Field date is contiguous");

        if (!Pe_.isUsable())
            throw std::runtime_error("Error - isothermal closure pressure not usable");

        auto const& Ne_ = ions_.chargeDensity();
        std::transform(std::begin(Ne_), std::end(Ne_), std::begin(Pe_),
                       [this](auto n) { return n * Te_; });
    }

    // Ohm's electron pressure term on the E edges of this patch (see pressure_term.hpp).
    // Electrons::update calls it after computePressure; code that runs the electron steps one by
    // one must call it too.
    void computePressureTerm(GridLayout const& layout)
    {
        isothermalPressureTerm_(layout, ions_.chargeDensity(), pressureTerm_);
    }

    NO_DISCARD VecField const& pressureTerm() const { return pressureTerm_; }

private:
    Ions ions_;
    double const Te_ = 0;
    Field Pe_;
    VecField pressureTerm_;
    IsothermalPressureTerm<GridLayout> isothermalPressureTerm_;
};




template<typename Ions>
class ElectronMomentModel
{
    using VecField     = Ions::vecfield_type;
    using Field        = Ions::field_type;
    using GridLayout   = Ions::gridlayout_type;
    using FluxComputer = StandardHybridElectronFluxComputer<Ions>;

public:
    ElectronMomentModel(PHARE::initializer::PHAREDict const& dict, Ions& ions, VecField& J)
        : fluxComput_{ions, J}
        , pressureClosure_{dict["pressure_closure"], fluxComput_.getIons()}
    {
    }

    //-------------------------------------------------------------------------
    //                  start the ResourcesUser interface
    //-------------------------------------------------------------------------

    NO_DISCARD bool isUsable() const
    {
        return fluxComput_.isUsable() and pressureClosure_.isUsable();
    }

    NO_DISCARD bool isSettable() const { return fluxComput_.isSettable(); }

    NO_DISCARD auto getCompileTimeResourcesViewList() const
    {
        return std::forward_as_tuple(fluxComput_, pressureClosure_);
    }

    NO_DISCARD auto getCompileTimeResourcesViewList()
    {
        return std::forward_as_tuple(fluxComput_, pressureClosure_);
    }

    //-------------------------------------------------------------------------
    //                  ends the ResourcesUser interface
    //-------------------------------------------------------------------------


    NO_DISCARD Field const& density() const { return fluxComput_.density(); }

    NO_DISCARD Field& density() { return fluxComput_.density(); }


    NO_DISCARD VecField const& velocity() const { return fluxComput_.velocity(); }

    NO_DISCARD VecField& velocity() { return fluxComput_.velocity(); }


    NO_DISCARD Field const& pressure() const { return pressureClosure_.pressure(); }
    NO_DISCARD Field& pressure() { return pressureClosure_.pressure(); }



    void computeDensity() { fluxComput_.computeDensity(); }
    void computeBulkVelocity(GridLayout const& layout) { fluxComput_.computeBulkVelocity(layout); }
    void computePressure(GridLayout const& layout) { pressureClosure_.computePressure(layout); }

    void computePressureTerm(GridLayout const& layout)
    {
        pressureClosure_.computePressureTerm(layout);
    }

    NO_DISCARD VecField const& pressureTerm() const { return pressureClosure_.pressureTerm(); }

private:
    FluxComputer fluxComput_;
    IsothermalElectronPressureClosure<Ions> pressureClosure_;
};



template<typename Ions>
class Electrons
{
    using VecField   = Ions::vecfield_type;
    using Field      = Ions::field_type;
    using GridLayout = Ions::gridlayout_type;

public:
    Electrons(PHARE::initializer::PHAREDict const& dict, Ions& ions, VecField& J)
        : momentModel_{dict, ions, J}
    {
    }


    void update(GridLayout const& layout)
    {
        if (!isUsable())
            throw std::runtime_error("Error - Electron is not usable");

        momentModel_.computeDensity();
        momentModel_.computeBulkVelocity(layout);
        momentModel_.computePressure(layout);
        momentModel_.computePressureTerm(layout);
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

    void computePressureTerm(GridLayout const& layout) { momentModel_.computePressureTerm(layout); }
    NO_DISCARD VecField const& pressureTerm() const { return momentModel_.pressureTerm(); }

private:
    ElectronMomentModel<Ions> momentModel_;
};

} // namespace PHARE::core


#endif
