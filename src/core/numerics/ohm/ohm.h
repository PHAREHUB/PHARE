#ifndef PHARE_OHM_H
#define PHARE_OHM_H

#include <cstddef>
#include <iostream>

#include "core/data/grid/gridlayoutdefs.h"
#include "core/data/grid/gridlayout.h"
#include "core/data/grid/gridlayout_utils.h"
#include "core/data/vecfield/vecfield_component.h"

#include "initializer/data_provider.h"


namespace PHARE::core
{
template<typename GridLayout>
struct StandardOhmComputer
{
    constexpr static auto dimension = GridLayout::dimension;

    template<typename VecField, typename ComponentTag>
    auto ideal1D_(VecField const& Ve, VecField const& B, MeshIndex<1> index, ComponentTag) const
    {
        if constexpr (ComponentTag::component == Component::X)
        {
            auto const& Vy = Ve.getComponent(Component::Y);
            auto const& Vz = Ve.getComponent(Component::Z);

            auto const& By = B.getComponent(Component::Y);
            auto const& Bz = B.getComponent(Component::Z);

            auto constexpr momentsToEx = GridLayout::momentsToEx();
            auto const vyOnEx          = GridLayout::project(Vy, index, momentsToEx);
            auto const vzOnEx          = GridLayout::project(Vz, index, momentsToEx);
            auto const byOnEx          = GridLayout::project(By, index, GridLayout::ByToEx());
            auto const bzOnEx          = GridLayout::project(Bz, index, GridLayout::BzToEx());

            return -vyOnEx * bzOnEx + vzOnEx * byOnEx;
        }

        if constexpr (ComponentTag::component == Component::Y)
        {
            auto const& Vx = Ve.getComponent(Component::X);
            auto const& Vz = Ve.getComponent(Component::Z);
            auto const& Bx = B.getComponent(Component::X);
            auto const& Bz = B.getComponent(Component::Z);

            auto constexpr momentsToEy = GridLayout::momentsToEy();
            auto const vxOnEy          = GridLayout::project(Vx, index, momentsToEy);
            auto const vzOnEy          = GridLayout::project(Vz, index, momentsToEy);
            auto const bxOnEy          = GridLayout::project(Bx, index, GridLayout::BxToEy());
            auto const bzOnEy          = GridLayout::project(Bz, index, GridLayout::BzToEy());

            return -vzOnEy * bxOnEy + vxOnEy * bzOnEy;
        }

        if constexpr (ComponentTag::component == Component::Z)
        {
            auto const& Vx = Ve.getComponent(Component::X);
            auto const& Vy = Ve.getComponent(Component::Y);
            auto const& Bx = B.getComponent(Component::X);
            auto const& By = B.getComponent(Component::Y);

            auto constexpr momentsToEz = GridLayout::momentsToEz();
            auto const vxOnEz          = GridLayout::project(Vx, index, momentsToEz);
            auto const vyOnEz          = GridLayout::project(Vy, index, momentsToEz);
            auto const bxOnEz          = GridLayout::project(Bx, index, GridLayout::BxToEz());
            auto const byOnEz          = GridLayout::project(By, index, GridLayout::ByToEz());

            return -vxOnEz * byOnEz + vyOnEz * bxOnEz;
        }
    }


    template<typename VecField, typename ComponentTag>
    auto ideal2D_(VecField const& Ve, VecField const& B, MeshIndex<2> index, ComponentTag) const
    {
        if constexpr (ComponentTag::component == Component::X)
        {
            auto const& Vy = Ve.getComponent(Component::Y);
            auto const& Vz = Ve.getComponent(Component::Z);
            auto const& By = B.getComponent(Component::Y);
            auto const& Bz = B.getComponent(Component::Z);

            auto constexpr momentsToEx = GridLayout::momentsToEx();
            auto const vyOnEx          = GridLayout::project(Vy, index, momentsToEx);
            auto const vzOnEx          = GridLayout::project(Vz, index, momentsToEx);
            auto const byOnEx          = GridLayout::project(By, index, GridLayout::ByToEx());
            auto const bzOnEx          = GridLayout::project(Bz, index, GridLayout::BzToEx());

            return -vyOnEx * bzOnEx + vzOnEx * byOnEx;
        }


        if constexpr (ComponentTag::component == Component::Y)
        {
            auto const& Vx = Ve.getComponent(Component::X);
            auto const& Vz = Ve.getComponent(Component::Z);
            auto const& Bx = B.getComponent(Component::X);
            auto const& Bz = B.getComponent(Component::Z);

            auto constexpr momentsToEy = GridLayout::momentsToEy();
            auto const vxOnEy          = GridLayout::project(Vx, index, momentsToEy);
            auto const vzOnEy          = GridLayout::project(Vz, index, momentsToEy);
            auto const bxOnEy          = GridLayout::project(Bx, index, GridLayout::BxToEy());
            auto const bzOnEy          = GridLayout::project(Bz, index, GridLayout::BzToEy());

            return -vzOnEy * bxOnEy + vxOnEy * bzOnEy;
        }

        if constexpr (ComponentTag::component == Component::Z)
        {
            auto const& Vx = Ve.getComponent(Component::X);
            auto const& Vy = Ve.getComponent(Component::Y);
            auto const& Bx = B.getComponent(Component::X);
            auto const& By = B.getComponent(Component::Y);

            auto constexpr momentsToEz = GridLayout::momentsToEz();
            auto const vxOnEz          = GridLayout::project(Vx, index, momentsToEz);
            auto const vyOnEz          = GridLayout::project(Vy, index, momentsToEz);
            auto const bxOnEz          = GridLayout::project(Bx, index, GridLayout::BxToEz());
            auto const byOnEz          = GridLayout::project(By, index, GridLayout::ByToEz());

            return -vxOnEz * byOnEz + vyOnEz * bxOnEz;
        }
    }



    template<typename VecField, typename ComponentTag>
    auto ideal3D_(VecField const& Ve, VecField const& B, MeshIndex<3> index, ComponentTag) const
    {
        if constexpr (ComponentTag::component == Component::X)
        {
            auto const& Vy = Ve.getComponent(Component::Y);
            auto const& Vz = Ve.getComponent(Component::Z);
            auto const& By = B.getComponent(Component::Y);
            auto const& Bz = B.getComponent(Component::Z);

            auto constexpr momentsToEx = GridLayout::momentsToEx();
            auto const vyOnEx          = GridLayout::project(Vy, index, momentsToEx);
            auto const vzOnEx          = GridLayout::project(Vz, index, momentsToEx);
            auto const byOnEx          = GridLayout::project(By, index, GridLayout::ByToEx());
            auto const bzOnEx          = GridLayout::project(Bz, index, GridLayout::BzToEx());

            return -vyOnEx * bzOnEx + vzOnEx * byOnEx;
        }

        if constexpr (ComponentTag::component == Component::Y)
        {
            auto const& Vx = Ve.getComponent(Component::X);
            auto const& Vz = Ve.getComponent(Component::Z);
            auto const& Bx = B.getComponent(Component::X);
            auto const& Bz = B.getComponent(Component::Z);

            auto constexpr momentsToEy = GridLayout::momentsToEy();
            auto const vxOnEy          = GridLayout::project(Vx, index, momentsToEy);
            auto const vzOnEy          = GridLayout::project(Vz, index, momentsToEy);
            auto const bxOnEy          = GridLayout::project(Bx, index, GridLayout::BxToEy());
            auto const bzOnEy          = GridLayout::project(Bz, index, GridLayout::BzToEy());

            return -vzOnEy * bxOnEy + vxOnEy * bzOnEy;
        }

        if constexpr (ComponentTag::component == Component::Z)
        {
            auto const& Vx = Ve.getComponent(Component::X);
            auto const& Vy = Ve.getComponent(Component::Y);
            auto const& Bx = B.getComponent(Component::X);
            auto const& By = B.getComponent(Component::Y);

            auto constexpr momentsToEz = GridLayout::momentsToEz();
            auto const vxOnEz          = GridLayout::project(Vx, index, momentsToEz);
            auto const vyOnEz          = GridLayout::project(Vy, index, momentsToEz);
            auto const bxOnEz          = GridLayout::project(Bx, index, GridLayout::BxToEz());
            auto const byOnEz          = GridLayout::project(By, index, GridLayout::ByToEz());

            return -vxOnEz * byOnEz + vyOnEz * bxOnEz;
        }
    }


    template<typename VecField, typename ComponentTag>
    auto ideal_(VecField const& Ve, VecField const& B, MeshIndex<dimension> index,
                ComponentTag tag) const
    {
        if constexpr (dimension == 1)
            return ideal1D_(Ve, B, index, tag);
        if constexpr (dimension == 2)
            return ideal2D_(Ve, B, index, tag);
        if constexpr (dimension == 3)
            return ideal3D_(Ve, B, index, tag);
    }

    template<typename Field, typename ComponentTag>
    auto pressure_(Field const& n, Field const& Pe, MeshIndex<Field::dimension> index,
                   ComponentTag) const
    {
        static_assert(ComponentTag::component == Component::X
                      || ComponentTag::component == Component::Y
                      || ComponentTag::component == Component::Z);

        if constexpr (ComponentTag::component == Component::X)
        {
            auto const nOnEx = GridLayout::project(n, index, GridLayout::momentsToEx());

            auto gradPOnEx
                = layout.deriv(Pe, index, DirectionTag<Direction::X>{}); // TODO : issue 3391

            return -gradPOnEx / nOnEx;
        }

        else if constexpr (ComponentTag::component == Component::Y)
        {
            if constexpr (Field::dimension >= 2)
            {
                auto const nOnEy = GridLayout::project(n, index, GridLayout::momentsToEy());

                auto gradPOnEy
                    = layout.deriv(Pe, index, DirectionTag<Direction::Y>{}); // TODO : issue 3391

                return -gradPOnEy / nOnEy;
            }
            else
            {
                return 0.;
            }
        }

        else if constexpr (ComponentTag::component == Component::Z)
        {
            if constexpr (Field::dimension >= 3)
            {
                auto const nOnEz = GridLayout::project(n, index, GridLayout::momentsToEz());

                auto gradPOnEz
                    = layout.deriv(Pe, index, DirectionTag<Direction::Z>{}); // TODO : issue 3391

                return -gradPOnEz / nOnEz;
            }
            else
            {
                return 0.;
            }
        }
    }




    template<typename VecField, typename ComponentTag>
    auto resistive_(VecField const& J, MeshIndex<VecField::dimension> index, ComponentTag) const
    {
        if constexpr (ComponentTag::component == Component::X)
        {
            auto const& Jx    = J.getComponent(Component::X);
            auto const jxOnEx = GridLayout::project(Jx, index, GridLayout::JxToEx());
            return eta * jxOnEx;
        }

        if constexpr (ComponentTag::component == Component::Y)
        {
            auto const& Jy    = J.getComponent(Component::Y);
            auto const jyOnEy = GridLayout::project(Jy, index, GridLayout::JyToEy());
            return eta * jyOnEy;
        }

        if constexpr (ComponentTag::component == Component::Z)
        {
            auto const& Jz    = J.getComponent(Component::Z);
            auto const jzOnEz = GridLayout::project(Jz, index, GridLayout::JzToEz());
            return eta * jzOnEz;
        }
    }




    template<typename VecField, typename ComponentTag>
    auto hyperresistive_(VecField const& J, MeshIndex<VecField::dimension> index,
                         ComponentTag) const
    {
        if constexpr (ComponentTag::component == Component::X)
        {
            auto const& Jx = J.getComponent(Component::X);
            auto lapJx     = layout.laplacian(Jx, index); // TODO : issue 3391
            return -nu * lapJx;
        }

        if constexpr (ComponentTag::component == Component::Y)
        {
            auto const& Jy = J.getComponent(Component::Y);
            auto lapJy     = layout.laplacian(Jy, index); // TODO : issue 3391
            return -nu * lapJy;
        }

        if constexpr (ComponentTag::component == Component::Z)
        {
            auto const& Jz = J.getComponent(Component::Z);
            auto lapJz     = layout.laplacian(Jz, index); // TODO : issue 3391
            return -nu * lapJz;
        }
    }



    template<typename Component, Component Tag, typename OhmPack, typename... IDXs>
    void exyz(OhmPack&& pack, IDXs const&... ijk) const
    {
        auto const& [E, n, Pe, Ve, B, J] = pack;
        auto& Exyz                       = E(Tag);

        Exyz(ijk...) = ideal_(Ve, B, {ijk...}, ComponentTag<Tag>{})
                       + pressure_(n, Pe, {ijk...}, ComponentTag<Tag>{})
                       + resistive_(J, {ijk...}, ComponentTag<Tag>{})
                       + hyperresistive_(J, {ijk...}, ComponentTag<Tag>{});
    }

    double eta;
    double nu;
    GridLayout& layout;
};


template<typename GridLayout, typename Computer = StandardOhmComputer<GridLayout>>
class Ohm : public LayoutHolder<GridLayout>
{
    using LayoutHolder<GridLayout>::layout_;

    template<typename VecField_ref, typename Field, typename VecField_cref>
    struct OhmPack
    {
        VecField_ref& Exyz;
        Field const &n, &Pe;
        VecField_cref const &Ve, &B, &J;
    };

public:
    explicit Ohm(PHARE::initializer::PHAREDict const& dict)
        : eta_{dict["resistivity"].template to<double>()}
        , nu_{dict["hyper_resistivity"].template to<double>()}
    {
    }

    template<typename VecField>
    void operator()(typename VecField::field_type const& n_, VecField const& Ve_,
                    typename VecField::field_type const& Pe_, VecField const& B_,
                    VecField const& J_, VecField& Enew_)
    {
        using Pack
            = OhmPack<decltype(Enew_.as_view()), decltype(n_.as_view()), decltype(Ve_.as_view())>;

        if (!this->hasLayout())
            throw std::runtime_error(
                "Error - Ohm - GridLayout not set, cannot proceed to calculate ohm()");

        auto n    = n_.as_view();
        auto Ve   = Ve_.as_view();
        auto Pe   = Pe_.as_view();
        auto B    = B_.as_view();
        auto J    = J_.as_view();
        auto Enew = Enew_.as_view();

        Computer op{eta_, nu_, *this->layout_};
        layout_->scan(Enew_(Component::X), [=](auto const&... args) mutable {
            op.template exyz<Component, Component::X>(Pack{Enew, n, Pe, Ve, B, J}, args...);
        });
        layout_->scan(Enew_(Component::Y), [=](auto const&... args) mutable {
            op.template exyz<Component, Component::Y>(Pack{Enew, n, Pe, Ve, B, J}, args...);
        });
        layout_->scan(Enew_(Component::Z), [=](auto const&... args) mutable {
            op.template exyz<Component, Component::Z>(Pack{Enew, n, Pe, Ve, B, J}, args...);
        });
    }

private:
    double const eta_;
    double const nu_;
};


} // namespace PHARE::core
#endif
