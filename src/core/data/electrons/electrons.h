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
template<typename VecField, typename GridLayout>
class FluxComputer
{
    using Field = typename VecField::field_type;

public:
    virtual bool isUsable() const                        = 0;
    virtual ~FluxComputer()                              = default;
    virtual Field& density()                             = 0;
    virtual Field const& density() const                 = 0;
    virtual VecField& velocity()                         = 0;
    virtual VecField const& velocity() const             = 0;
    virtual void computeDensity()                        = 0;
    virtual void computeBulkVelocity(GridLayout& layout) = 0;

private:
};


template<typename Ions>
class StandardHybridElectronFluxComputer
    : public FluxComputer<typename Ions::vecfield_type, typename Ions::gridlayout_type>
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
    virtual ~StandardHybridElectronFluxComputer() override = default;

    virtual bool isUsable() const override
    {
        return ions_.isUsable() && J_.isUsable() && Ne_ != nullptr && Ve_.isUsable();
    }
    virtual Field& density() override
    {
        if (isUsable())
        {
            return *Ne_;
        }
        else
        {
            throw std::runtime_error("Error, cannot return density because "
                                     "StandardHybridElectronFluxComputer is not usable");
        }
    }

    virtual Field const& density() const override
    {
        if (isUsable())
        {
            return *Ne_;
        }
        else
        {
            throw std::runtime_error("Error, cannot return density because "
                                     "StandardHybridElectronFluxComputer is not usable");
        }
    }

    virtual VecField const& velocity() const override
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


    virtual VecField& velocity() override
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

    virtual void computeDensity() override
    {
        auto& Ni = ions_.density();
        Ne_      = &Ni;
    }


    virtual void computeBulkVelocity(GridLayout& layout) override
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
    Field* Ne_;
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
        , fluxComput_{std::make_unique<StandardHybridElectronFluxComputer<Ions>>(ions, J)}
    {
    }

    bool isUsable() const { return fluxComput_->isUsable() && electromag_.isUsable(); }

    // Field& density() { return *Ne_; }

    // VecField& velocity() { return Ve_; }

    Field& pressure() { return *Pe_; }

private:
    Field* Pe_;
    Electromag& electromag_;
    std::unique_ptr<FluxComputer<VecField, GridLayout>> fluxComput_;
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
            // momentModel_.fluxComputer->computeDensity(layout_);
        }
    }

    bool isUsable() const { return momentModel_.isUsable(); }

    Field& density() { return momentModel_.density(); }
    Field& velocity() { return momentModel_.velocity(); }
    Field& pressure() { return momentModel_.pressure(); }

private:
    ElectronMomentModel<Ions, Electromag> momentModel_;
};

} // namespace PHARE::core


#endif
