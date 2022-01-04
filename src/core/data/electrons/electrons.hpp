#ifndef PHARE_ELECTRONS_H
#define PHARE_ELECTRONS_H

#include "core/hybrid/hybrid_quantities.hpp"
#include "core/data/vecfield/vecfield_component.hpp"
#include "core/data/grid/gridlayout_utils.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/utilities/index/index.hpp"

#include "initializer/data_provider.hpp"
#include <memory>


namespace PHARE::core
{
template<typename Ions>
class StandardHybridElectronFluxComputer
{
public:
    using VecField   = typename Ions::vecfield_type;
    using Field      = typename Ions::field_type;
    using GridLayout = typename Ions::gridlayout_type;

    StandardHybridElectronFluxComputer(Ions& ions, VecField& J)
        : ions_{ions}
        , J_{J}
        , Ve_{"StandardHybridElectronFluxComputer_Ve", HybridQuantity::Vector::V}
    {
    }

    //-------------------------------------------------------------------------
    //                  start the ResourcesUser interface
    //-------------------------------------------------------------------------


    bool isUsable() const { return ions_.isUsable() && J_.isUsable() && Ve_.isUsable(); }

    bool isSettable() const { return Ve_.isSettable() && ions_.isSettable() && J_.isSettable(); }

    auto getCompileTimeResourcesUserList() const { return std::forward_as_tuple(Ve_, ions_, J_); }

    auto getCompileTimeResourcesUserList() { return std::forward_as_tuple(Ve_, ions_, J_); }

    //-------------------------------------------------------------------------
    //                  ends the ResourcesUser interface
    //-------------------------------------------------------------------------




    Field const& density() const
    {
        if (isUsable())
        {
            return ions_.density();
        }
        else
        {
            throw std::runtime_error("Error, cannot return density because "
                                     "StandardHybridElectronFluxComputer is not usable");
        }
    }

    Field& density()
    {
        if (isUsable())
        {
            return ions_.density();
        }
        else
        {
            throw std::runtime_error("Error, cannot return density because "
                                     "StandardHybridElectronFluxComputer is not usable");
        }
    }


    VecField& velocity()
    {
        if (isUsable())
        {
            return Ve_;
        }
        else
        {
            throw std::runtime_error("Error, cannot return velocity because "
                                     "StandardHybridElectronFluxComputer is not usable");
        }
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
        auto const& Ni  = ions_.density();

        auto& Vex = Ve_(Component::X);
        auto& Vey = Ve_(Component::Y);
        auto& Vez = Ve_(Component::Z);

        // from Ni because all components defined on primal
        layout.evalOnBox(Ni, [&](auto const&... args) {
            auto arr = std::array{args...};

            auto const JxOnVx = GridLayout::project(Jx, arr, GridLayout::JxToMoments());
            auto const JyOnVy = GridLayout::project(Jy, arr, GridLayout::JyToMoments());
            auto const JzOnVz = GridLayout::project(Jz, arr, GridLayout::JzToMoments());

            Vex(arr) = Vix(arr) - JxOnVx / Ni(arr);
            Vey(arr) = Viy(arr) - JyOnVy / Ni(arr);
            Vez(arr) = Viz(arr) - JzOnVz / Ni(arr);
        });
    }



private:
    Ions& ions_;
    VecField& J_;
    VecField Ve_;

}; // namespace PHARE::core



template<typename FluxComputer>
class IsothermalElectronPressureClosure
{
    using GridLayout = typename FluxComputer::GridLayout;
    using VecField   = typename FluxComputer::VecField;
    using Field      = typename VecField::field_type;

public:
    IsothermalElectronPressureClosure(PHARE::initializer::PHAREDict const& dict,
                                      FluxComputer const& fc)
        : fluxComputer_{fc}
        , Te_{dict["Te"].template to<double>()}
    {
    }


    //-------------------------------------------------------------------------
    //                  start the ResourcesUser interface
    //-------------------------------------------------------------------------

    using field_type = Field;

    bool isUsable() { return Pe_ != nullptr; }

    bool isSettable() { return Pe_ == nullptr; }

    struct PressureProperty
    {
        std::string name;
        typename HybridQuantity::Scalar qty;
    };

    using PressureProperties = std::vector<PressureProperty>;

    PressureProperties getFieldNamesAndQuantities() const
    {
        return {{{"Pe", HybridQuantity::Scalar::P}}};
    }



    void setBuffer(std::string name, Field* field)
    {
        if (name == "Pe")
            Pe_ = field;
        else
            throw std::runtime_error("Error - Pe is not usable");
    }


    //-------------------------------------------------------------------------
    //                  ends the ResourcesUser interface
    //-------------------------------------------------------------------------



    Field& pressure()
    {
        if (Pe_ != nullptr)
            return *Pe_;
        else
            throw std::runtime_error("Error - isothermal closure pressure not usable");
    }
    Field const& pressure() const { return *Pe_; }

    void computePressure([[maybe_unused]] GridLayout const& layout)
    {
        if (Pe_ != nullptr)
        {
            static_assert(Field::is_contiguous, "Error - assumes Field date is contiguous");

            auto const& Ne_ = fluxComputer_.density();

            std::transform(std::begin(Ne_), std::end(Ne_), std::begin(*Pe_),
                           [this](auto n) { return n * Te_; });
        }
        else
            throw std::runtime_error("Error - isothermal closure pressure not usable");
    }

private:
    FluxComputer const& fluxComputer_;
    const double Te_;
    Field* Pe_;
};




template<typename Ions>
class ElectronMomentModel
{
    using VecField     = typename Ions::vecfield_type;
    using Field        = typename Ions::field_type;
    using GridLayout   = typename Ions::gridlayout_type;
    using FluxComputer = StandardHybridElectronFluxComputer<Ions>;

public:
    ElectronMomentModel(PHARE::initializer::PHAREDict const& dict, Ions& ions, VecField& J)
        : fluxComput_{ions, J}
        , pressureClosure_{dict["pressure_closure"], fluxComput_}
    {
    }

    //-------------------------------------------------------------------------
    //                  start the ResourcesUser interface
    //-------------------------------------------------------------------------

    bool isUsable() const { return fluxComput_.isUsable(); }

    bool isSettable() const { return fluxComput_.isSettable(); }

    auto getCompileTimeResourcesUserList() const
    {
        return std::forward_as_tuple(fluxComput_, pressureClosure_);
    }

    auto getCompileTimeResourcesUserList()
    {
        return std::forward_as_tuple(fluxComput_, pressureClosure_);
    }

    //-------------------------------------------------------------------------
    //                  ends the ResourcesUser interface
    //-------------------------------------------------------------------------


    Field const& density() const { return fluxComput_.density(); }

    Field& density() { return fluxComput_.density(); }


    VecField const& velocity() const { return fluxComput_.velocity(); }

    VecField& velocity() { return fluxComput_.velocity(); }


    Field const& pressure() const { return pressureClosure_.pressure(); }
    Field& pressure() { return pressureClosure_.pressure(); }



    void computeDensity() { fluxComput_.computeDensity(); }
    void computeBulkVelocity(GridLayout const& layout) { fluxComput_.computeBulkVelocity(layout); }
    void computePressure(GridLayout const& layout) { pressureClosure_.computePressure(layout); }

private:
    FluxComputer fluxComput_;
    IsothermalElectronPressureClosure<FluxComputer> pressureClosure_;
};



template<typename Ions>
class Electrons : public LayoutHolder<typename Ions::gridlayout_type>
{
    using VecField   = typename Ions::vecfield_type;
    using Field      = typename Ions::field_type;
    using GridLayout = typename Ions::gridlayout_type;

public:
    Electrons(PHARE::initializer::PHAREDict const& dict, Ions& ions, VecField& J)
        : momentModel_{dict, ions, J}
    {
    }

    void update(GridLayout const& layout)
    {
        if (isUsable())
        {
            momentModel_.computeDensity();
            momentModel_.computeBulkVelocity(layout);
            momentModel_.computePressure(layout);
        }
        else
            throw std::runtime_error("Error - Electron  is not usable");
    }


    //-------------------------------------------------------------------------
    //                  start the ResourcesUser interface
    //-------------------------------------------------------------------------

    bool isUsable() const { return momentModel_.isUsable(); }

    bool isSettable() const { return momentModel_.isSettable(); }

    auto getCompileTimeResourcesUserList() const { return std::forward_as_tuple(momentModel_); }

    auto getCompileTimeResourcesUserList() { return std::forward_as_tuple(momentModel_); }



    //-------------------------------------------------------------------------
    //                  ends the ResourcesUser interface
    //-------------------------------------------------------------------------




    Field const& density() const { return momentModel_.density(); }
    VecField const& velocity() const { return momentModel_.velocity(); }
    Field const& pressure() const { return momentModel_.pressure(); }


    Field& density() { return momentModel_.density(); }
    VecField& velocity() { return momentModel_.velocity(); }
    Field& pressure() { return momentModel_.pressure(); }

private:
    ElectronMomentModel<Ions> momentModel_;
};

} // namespace PHARE::core


#endif
