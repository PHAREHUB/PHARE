#ifndef PHARE_ELECTRONS_HPP
#define PHARE_ELECTRONS_HPP

#include "core/def.hpp"
#include "core/models/quantities/hybrid_quantities.hpp"
#include "core/data/vecfield/vecfield_component.hpp"
#include "core/utilities/types.hpp"
#include "core/utilities/point/point.hpp"
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

    static constexpr auto dimension = GridLayout::dimension;

    StandardHybridElectronFluxComputer(Ions& ions, VecField& J)
        : ions_{ions}
        , J_{J}
        , Ge_{"StandardHybridElectronFluxComputer_flux", HybridQuantity::Vector::J}
    {
    }

    //-------------------------------------------------------------------------
    //                  start the ResourcesUser interface
    //-------------------------------------------------------------------------


    NO_DISCARD bool isUsable() const { return ions_.isUsable() && J_.isUsable() && Ge_.isUsable(); }

    NO_DISCARD bool isSettable() const
    {
        return Ge_.isSettable() && ions_.isSettable() && J_.isSettable();
    }

    NO_DISCARD auto getCompileTimeResourcesViewList() const
    {
        return std::forward_as_tuple(Ge_, ions_, J_);
    }

    NO_DISCARD auto getCompileTimeResourcesViewList()
    {
        return std::forward_as_tuple(Ge_, ions_, J_);
    }

    //-------------------------------------------------------------------------
    //                  ends the ResourcesUser interface
    //-------------------------------------------------------------------------




    NO_DISCARD Field const& density() const { return ions_.chargeDensity(); }

    NO_DISCARD Field& density() { return ions_.chargeDensity(); }


    // electron flux G = n Ve on the E nodes (where J and E live), what Ohm's ideal term reads
    NO_DISCARD VecField const& flux() const { return Ge_; }
    NO_DISCARD VecField& flux() { return Ge_; }


    void computeDensity() {}

    // G = n Ve = Ji - J, with Ji = sum_s q_s fluxE_s the ion current deposited on the E nodes
    // (Ions::chargeFluxE) and J the total current: J = Ji + Je and Je = -n Ve. It is charge
    // weighted, not mass weighted (#1466). Computed one node beyond the domain because Ohm reads
    // it there (#1577); Ji and J are valid there (Ampere fills J on the ghost box shrunk by 1).
    void computeFlux(GridLayout const& layout)
    {
        auto const& Ji = ions_.chargeFluxE();
        auto const one = Point<std::uint32_t, dimension>{ConstArray<std::uint32_t, dimension>(1)};
        for (auto const c : {Component::X, Component::Y, Component::Z})
        {
            auto& Gc        = Ge_(c);
            auto const& Jic = Ji(c);
            auto const& Jc  = J_(c);
            layout.evalOnBiggerBox(
                Gc, one, [&](auto const&... ijk) { Gc(ijk...) = Jic(ijk...) - Jc(ijk...); });
        }
    }

    // Ve = G / n on the moment nodes, G averaged from the E nodes, for code that needs the
    // electron velocity (e.g. a non-isothermal closure). Nothing in the time loop does, so it is
    // not stored: the caller owns Ve (HybridQuantity::Vector::V) and calls this after
    // computeFlux. Domain nodes only.
    void computeVelocity(GridLayout const& layout, VecField& Ve) const
    {
        auto const& Ne              = ions_.chargeDensity();
        auto const& [Gx, Gy, Gz]    = Ge_();
        auto const& [Vex, Vey, Vez] = Ve(); // a tuple of references: Vex is writable
        layout.evalOnBox(Ne, [&](auto const&... args) {
            auto const arr = std::array{args...};
            Vex(arr) = GridLayout::template project<Scalar::Jx, Scalar::Vx>(Gx, arr) / Ne(arr);
            Vey(arr) = GridLayout::template project<Scalar::Jy, Scalar::Vy>(Gy, arr) / Ne(arr);
            Vez(arr) = GridLayout::template project<Scalar::Jz, Scalar::Vz>(Gz, arr) / Ne(arr);
        });
    }


    auto& getIons() const { return ions_; }

private:
    Ions ions_;
    VecField J_;
    VecField Ge_;
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


    NO_DISCARD bool isUsable() const { return core::isUsable(Pe_, pressureTerm_, ions_); }

    NO_DISCARD bool isSettable() const { return core::isSettable(Pe_, pressureTerm_); }

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


    NO_DISCARD VecField const& flux() const { return fluxComput_.flux(); }

    NO_DISCARD VecField& flux() { return fluxComput_.flux(); }

    void computeVelocity(GridLayout const& layout, VecField& Ve) const
    {
        fluxComput_.computeVelocity(layout, Ve);
    }


    NO_DISCARD Field const& pressure() const { return pressureClosure_.pressure(); }
    NO_DISCARD Field& pressure() { return pressureClosure_.pressure(); }



    void computeDensity() { fluxComput_.computeDensity(); }
    void computeFlux(GridLayout const& layout) { fluxComput_.computeFlux(layout); }
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
        momentModel_.computeFlux(layout);
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
    NO_DISCARD VecField const& flux() const { return momentModel_.flux(); }
    NO_DISCARD Field const& pressure() const { return momentModel_.pressure(); }


    NO_DISCARD Field& density() { return momentModel_.density(); }
    NO_DISCARD VecField& flux() { return momentModel_.flux(); }
    NO_DISCARD Field& pressure() { return momentModel_.pressure(); }

    // electron velocity on the moment nodes into a caller-owned Ve, on demand (not stored)
    void computeVelocity(GridLayout const& layout, VecField& Ve) const
    {
        momentModel_.computeVelocity(layout, Ve);
    }

    void computePressureTerm(GridLayout const& layout) { momentModel_.computePressureTerm(layout); }
    NO_DISCARD VecField const& pressureTerm() const { return momentModel_.pressureTerm(); }

private:
    ElectronMomentModel<Ions> momentModel_;
};

} // namespace PHARE::core


#endif
