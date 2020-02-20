#ifndef PHARE_ELECTRONS_H
#define PHARE_ELECTRONS_H

#include "core/hybrid/hybrid_quantities.h"
#include "core/data/vecfield/vecfield_component.h"
#include "core/data/grid/gridlayout_utils.h"
#include "core/data/grid/gridlayoutdefs.h"

#include "initializer/data_provider.h"
#include <memory>


namespace PHARE::core
{
template<typename Ions>
class StandardHybridElectronFluxComputer
{
    using VecField   = typename Ions::vecfield_type;
    using Field      = typename Ions::field_type;
    using GridLayout = typename Ions::gridlayout_type;

public:
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

    VecField const& velocity() const
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

    void computeBulkVelocity(GridLayout& layout)
    {
        auto const& Ni = ions_.density();
        auto const& Vi = ions_.velocity();

        auto& Vex = Ve_.getComponent(Component::X);
        auto& Vey = Ve_.getComponent(Component::Y);
        auto& Vez = Ve_.getComponent(Component::Z);

        auto& Vix = Vi.getComponent(Component::X);
        auto& Viy = Vi.getComponent(Component::Y);
        auto& Viz = Vi.getComponent(Component::Z);

        if constexpr (GridLayout::dimension == 1)
        {
            auto ix0 = layout.physicalStartIndex(Vex, Direction::X);
            auto ix1 = layout.physicalEndIndex(Vex, Direction::X);

            for (auto ix = ix0; ix <= ix1; ++ix)
            {
                //                     auto const JxOnVx = GridLayout::project(Jx, index,
                //                     GridLayout::JxToMoment());

                Vex(ix) = Vix(ix) /* - JxOnVx/Ni(ix) */;
            }
        }
        else
        {
            throw std::runtime_error("computeBulkVelocity not implemented if not 1D");
        }
    }



private:
    Ions& ions_;
    VecField& J_;
    VecField Ve_;

}; // namespace PHARE::core




template<typename Ions, typename Electromag>
class ElectronMomentModel
{
    using VecField   = typename Ions::vecfield_type;
    using Field      = typename Ions::field_type;
    using GridLayout = typename Ions::gridlayout_type;

public:
    ElectronMomentModel(Ions& ions, Electromag& electromag, VecField& J)
        : electromag_{electromag}
        , fluxComput_{ions, J}
    {
    }

    //-------------------------------------------------------------------------
    //                  start the ResourcesUser interface
    //-------------------------------------------------------------------------

    bool isUsable() const { return fluxComput_.isUsable() && electromag_.isUsable(); }

    bool isSettable() const { return fluxComput_.isSettable(); }


    auto getCompileTimeResourcesUserList() const { return std::forward_as_tuple(fluxComput_); }

    auto getCompileTimeResourcesUserList() { return std::forward_as_tuple(fluxComput_); }

    //-------------------------------------------------------------------------
    //                  ends the ResourcesUser interface
    //-------------------------------------------------------------------------

    Field& density() { return fluxComput_.density(); }
    Field const& density() const { return fluxComput_.density(); }

    VecField& velocity() { return fluxComput_->velocity(); }
    VecField const& velocity() const { return fluxComput_.velocity(); }

    void computeDensity() { fluxComput_.computeDensity(); }

    Field& pressure() { return *Pe_; }

private:
    Field* Pe_;
    Electromag& electromag_;
    StandardHybridElectronFluxComputer<Ions> fluxComput_;
};



template<typename Ions, typename Electromag>
class Electrons : public LayoutHolder<typename Ions::gridlayout_type>
{
    using VecField = typename Ions::vecfield_type;
    using Field    = typename Ions::field_type;

public:
    Electrons(PHARE::initializer::PHAREDict dict, Ions& ions, Electromag& electromag, VecField& J)
        : momentModel_{ions, electromag, J}
    {
    }

    void update()
    {
        if (isUsable())
        {
            momentModel_.computeDensity();
        }
        else
            throw std::runtime_error("Errror - Electron  is not usable");
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


    Field& density() { return momentModel_.density(); }
    Field& velocity() { return momentModel_.velocity(); }
    Field& pressure() { return momentModel_.pressure(); }

private:
    ElectronMomentModel<Ions, Electromag> momentModel_;
};

} // namespace PHARE::core


#endif
