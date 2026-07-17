#ifndef PHARE_ELECTRONS_HPP
#define PHARE_ELECTRONS_HPP

#include "core/def.hpp"
#include "core/data/grid/grid_tiles.hpp"
#include "core/hybrid/hybrid_quantities.hpp"
#include "core/data/tensorfield/tensorfield.hpp"
#include "core/data/vecfield/vecfield_component.hpp"

#include "initializer/data_provider.hpp"


#include <tuple>


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




    NO_DISCARD Field const& density() const
    {
        if (isUsable())
        {
            return ions_.chargeDensity();
        }
        else
        {
            throw std::runtime_error("Error, cannot return density because "
                                     "StandardHybridElectronFluxComputer is not usable");
        }
    }

    NO_DISCARD Field& density()
    {
        if (isUsable())
        {
            return ions_.chargeDensity();
        }
        else
        {
            throw std::runtime_error("Error, cannot return density because "
                                     "StandardHybridElectronFluxComputer is not usable");
        }
    }


    NO_DISCARD VecField& velocity()
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


    void computeChargeDensity() {}



    void static Vxyz(auto const& layout, auto const& Ne, auto&&... args)
    {
        auto&& [J, Vi, Ve] = std::forward_as_tuple(args...);
        auto const& Jx     = J(Component::X);
        auto const& Jy     = J(Component::Y);
        auto const& Jz     = J(Component::Z);
        auto const& Vix    = Vi(Component::X);
        auto const& Viy    = Vi(Component::Y);
        auto const& Viz    = Vi(Component::Z);
        auto& Vex          = Ve(Component::X);
        auto& Vey          = Ve(Component::Y);
        auto& Vez          = Ve(Component::Z);

        layout.evalGrownOnBox(Ne, 1, [=] _PHARE_ALL_FN_(auto const& ijk) mutable {
            auto const JxOnVx = GridLayout::template project<GridLayout::JxToMoments>(Jx, ijk);
            auto const JyOnVy = GridLayout::template project<GridLayout::JyToMoments>(Jy, ijk);
            auto const JzOnVz = GridLayout::template project<GridLayout::JzToMoments>(Jz, ijk);

            Vex(ijk) = Vix(ijk) - JxOnVx / Ne(ijk);
            Vey(ijk) = Viy(ijk) - JyOnVy / Ne(ijk);
            Vez(ijk) = Viz(ijk) - JzOnVz / Ne(ijk);
        });
    }

    template<typename V_t>
    V_t static tt(auto& vf, auto i)
    {
        return vf.template as<V_t>([&](auto& c) { return c()[i](); });
    }

    void computeBulkVelocity(GridLayout const& layout)
    {
        check_field(ions_.chargeDensity(), layout);
        check_tensor_field(ions_.velocity(), layout);
        check_tensor_field(J_, layout);

        if constexpr (is_field_tile_set_v<Field>)
        {
            using Tile_vt = Field::value_type::value_type;
            using V_t     = basic::TensorField<Tile_vt, 1>;
            for (std::size_t tidx = 0; tidx < J_[0]().size(); ++tidx)
            {
                auto Ve = Ve_.template as<V_t>([&](auto& c) { return c()[tidx](); });
                Vxyz(J_[0]()[tidx].layout(), ions_.chargeDensity()()[tidx](), tt<V_t>(J_, tidx),
                     tt<V_t>(ions_.velocity(), tidx), Ve);
            }
        }
        else
        {
            Vxyz(layout, ions_.chargeDensity(), J_, ions_.velocity(), Ve_);
        }

        check_tensor_field(Ve_, layout);
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
    using GridLayout = typename Ions::gridlayout_type;
    using VecField   = typename Ions::vecfield_type;
    using Field      = typename Ions::field_type;

public:
    using field_type = Field;

    IsothermalElectronPressureClosure(PHARE::initializer::PHAREDict const& dict, Ions const& ions)
        : ions_{ions}
        , Te_{dict["Te"].template to<double>()}
        , Pe_{"Pe", HybridQuantity::Scalar::P}
    {
    }


    //-------------------------------------------------------------------------
    //                  start the ResourcesUser interface
    //-------------------------------------------------------------------------


    NO_DISCARD bool isUsable() const { return Pe_.isUsable() and ions_.isUsable(); }

    NO_DISCARD bool isSettable() const { return Pe_.isSettable(); }

    struct PressureProperty
    {
        std::string name;
        typename HybridQuantity::Scalar qty;
    };

    using PressureProperties = std::vector<PressureProperty>;

    NO_DISCARD auto getCompileTimeResourcesViewList() const
    {
        return std::forward_as_tuple(ions_, Pe_);
    }

    NO_DISCARD auto getCompileTimeResourcesViewList() { return std::forward_as_tuple(ions_, Pe_); }

    //-------------------------------------------------------------------------
    //                  ends the ResourcesUser interface
    //-------------------------------------------------------------------------



    NO_DISCARD Field& pressure()
    {
        if (!Pe_.isUsable())
            throw std::runtime_error("Error - isothermal closure pressure not usable");
        return Pe_;
    }
    NO_DISCARD Field const& pressure() const
    {
        if (!Pe_.isUsable())
            throw std::runtime_error("Error - isothermal closure pressure not usable");
        return Pe_;
    }

    void computePressure(GridLayout const& /*layout*/)
    {
        // static_assert(Field::is_contiguous, "Error - assumes Field date is contiguous");

        if (!Pe_.isUsable())
            throw std::runtime_error("Error - isothermal closure pressure not usable");

        auto const& Ne_ = ions_.chargeDensity();
        transform(Ne_, Pe_, [this](auto n) { return n * Te_; });
    }

private:
    Ions ions_;
    double const Te_ = 0;
    Field Pe_;
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



    void computeChargeDensity() { fluxComput_.computeChargeDensity(); }
    void computeBulkVelocity(GridLayout const& layout) { fluxComput_.computeBulkVelocity(layout); }
    void computePressure(GridLayout const& layout) { pressureClosure_.computePressure(layout); }

private:
    FluxComputer fluxComput_;
    IsothermalElectronPressureClosure<Ions> pressureClosure_;
};



template<typename Ions>
class Electrons
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
            momentModel_.computeChargeDensity();
            momentModel_.computeBulkVelocity(layout);
            momentModel_.computePressure(layout);
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
    ElectronMomentModel<Ions> momentModel_;
};

} // namespace PHARE::core


#endif
