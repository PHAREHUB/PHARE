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
class Ohm : public LayoutHolder<GridLayout>
{
    constexpr static auto dimension = GridLayout::dimension;
    using LayoutHolder<GridLayout>::layout_;

public:
    explicit Ohm(PHARE::initializer::PHAREDict const& dict)
        : eta_{dict["resistivity"].template to<double>()}
        , nu_{dict["hyper_resistivity"].template to<double>()}
    {
    }

    template<typename VecField, typename Field>
    void operator()(Field const& n, VecField const& Ve, Field const& Pe, VecField const& B,
                    VecField const& J, VecField& Enew)
    {
        using Pack = OhmPack<VecField, Field>;

        if (!this->hasLayout())
            throw std::runtime_error(
                "Error - Ohm - GridLayout not set, cannot proceed to calculate ohm()");

        auto const& [Exnew, Eynew, Eznew] = Enew();

        layout_->evalOnBox(Exnew, [&](auto&... args) mutable {
            this->template E_Eq_<Component, Component::X>(Pack{Enew, n, Pe, Ve, B, J}, args...);
        });
        layout_->evalOnBox(Eynew, [&](auto&... args) mutable {
            this->template E_Eq_<Component, Component::Y>(Pack{Enew, n, Pe, Ve, B, J}, args...);
        });
        layout_->evalOnBox(Eznew, [&](auto&... args) mutable {
            this->template E_Eq_<Component, Component::Z>(Pack{Enew, n, Pe, Ve, B, J}, args...);
        });
    }

    double const eta_;
    double const nu_;

private:
    template<typename VecField, typename Field>
    struct OhmPack
    {
        VecField& Exyz;
        Field const &n, &Pe;
        VecField const &Ve, &B, &J;
    };


    template<typename Component, Component Tag, typename OhmPack, typename... IDXs>
    void E_Eq_(OhmPack&& pack, IDXs const&... ijk) const
    {
        auto const& [E, n, Pe, Ve, B, J] = pack;
        auto& Exyz                       = E(Tag);

        ComponentTag<Tag> tag;
        static_assert(Components::check(tag));

        Exyz(ijk...) = ideal_(Ve, B, {ijk...}, tag) + pressure_(n, Pe, {ijk...}, tag)
                       + resistive_(J, {ijk...}, tag) + hyperresistive_(J, {ijk...}, tag);
    }



    template<typename VecField, typename ComponentTag>
    auto ideal1D_(VecField const& Ve, VecField const& B, MeshIndex<1> index, ComponentTag) const
    {
        if constexpr (ComponentTag::component == Component::X)
        {
            auto const& Vy = Ve(Component::Y);
            auto const& Vz = Ve(Component::Z);

            auto const& By = B(Component::Y);
            auto const& Bz = B(Component::Z);

            auto constexpr momentsToEx = GridLayout::momentsToEx();
            auto const vyOnEx          = GridLayout::project(Vy, index, momentsToEx);
            auto const vzOnEx          = GridLayout::project(Vz, index, momentsToEx);
            auto const byOnEx          = GridLayout::project(By, index, GridLayout::ByToEx());
            auto const bzOnEx          = GridLayout::project(Bz, index, GridLayout::BzToEx());

            return -vyOnEx * bzOnEx + vzOnEx * byOnEx;
        }

        if constexpr (ComponentTag::component == Component::Y)
        {
            auto const& Vx = Ve(Component::X);
            auto const& Vz = Ve(Component::Z);
            auto const& Bx = B(Component::X);
            auto const& Bz = B(Component::Z);

            auto constexpr momentsToEy = GridLayout::momentsToEy();
            auto const vxOnEy          = GridLayout::project(Vx, index, momentsToEy);
            auto const vzOnEy          = GridLayout::project(Vz, index, momentsToEy);
            auto const bxOnEy          = GridLayout::project(Bx, index, GridLayout::BxToEy());
            auto const bzOnEy          = GridLayout::project(Bz, index, GridLayout::BzToEy());

            return -vzOnEy * bxOnEy + vxOnEy * bzOnEy;
        }

        if constexpr (ComponentTag::component == Component::Z)
        {
            auto const& Vx = Ve(Component::X);
            auto const& Vy = Ve(Component::Y);
            auto const& Bx = B(Component::X);
            auto const& By = B(Component::Y);

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
            auto const& Vy = Ve(Component::Y);
            auto const& Vz = Ve(Component::Z);
            auto const& By = B(Component::Y);
            auto const& Bz = B(Component::Z);

            auto constexpr momentsToEx = GridLayout::momentsToEx();
            auto const vyOnEx          = GridLayout::project(Vy, index, momentsToEx);
            auto const vzOnEx          = GridLayout::project(Vz, index, momentsToEx);
            auto const byOnEx          = GridLayout::project(By, index, GridLayout::ByToEx());
            auto const bzOnEx          = GridLayout::project(Bz, index, GridLayout::BzToEx());

            return -vyOnEx * bzOnEx + vzOnEx * byOnEx;
        }


        if constexpr (ComponentTag::component == Component::Y)
        {
            auto const& Vx = Ve(Component::X);
            auto const& Vz = Ve(Component::Z);
            auto const& Bx = B(Component::X);
            auto const& Bz = B(Component::Z);

            auto constexpr momentsToEy = GridLayout::momentsToEy();
            auto const vxOnEy          = GridLayout::project(Vx, index, momentsToEy);
            auto const vzOnEy          = GridLayout::project(Vz, index, momentsToEy);
            auto const bxOnEy          = GridLayout::project(Bx, index, GridLayout::BxToEy());
            auto const bzOnEy          = GridLayout::project(Bz, index, GridLayout::BzToEy());

            return -vzOnEy * bxOnEy + vxOnEy * bzOnEy;
        }

        if constexpr (ComponentTag::component == Component::Z)
        {
            auto const& Vx = Ve(Component::X);
            auto const& Vy = Ve(Component::Y);
            auto const& Bx = B(Component::X);
            auto const& By = B(Component::Y);

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
            auto const& Vy = Ve(Component::Y);
            auto const& Vz = Ve(Component::Z);
            auto const& By = B(Component::Y);
            auto const& Bz = B(Component::Z);

            auto constexpr momentsToEx = GridLayout::momentsToEx();
            auto const vyOnEx          = GridLayout::project(Vy, index, momentsToEx);
            auto const vzOnEx          = GridLayout::project(Vz, index, momentsToEx);
            auto const byOnEx          = GridLayout::project(By, index, GridLayout::ByToEx());
            auto const bzOnEx          = GridLayout::project(Bz, index, GridLayout::BzToEx());

            return -vyOnEx * bzOnEx + vzOnEx * byOnEx;
        }

        if constexpr (ComponentTag::component == Component::Y)
        {
            auto const& Vx = Ve(Component::X);
            auto const& Vz = Ve(Component::Z);
            auto const& Bx = B(Component::X);
            auto const& Bz = B(Component::Z);

            auto constexpr momentsToEy = GridLayout::momentsToEy();
            auto const vxOnEy          = GridLayout::project(Vx, index, momentsToEy);
            auto const vzOnEy          = GridLayout::project(Vz, index, momentsToEy);
            auto const bxOnEy          = GridLayout::project(Bx, index, GridLayout::BxToEy());
            auto const bzOnEy          = GridLayout::project(Bz, index, GridLayout::BzToEy());

            return -vzOnEy * bxOnEy + vxOnEy * bzOnEy;
        }

        if constexpr (ComponentTag::component == Component::Z)
        {
            auto const& Vx = Ve(Component::X);
            auto const& Vy = Ve(Component::Y);
            auto const& Bx = B(Component::X);
            auto const& By = B(Component::Y);

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
                   ComponentTag tag) const
    {
        if constexpr (ComponentTag::component == Component::X)
        {
            auto const nOnEx = GridLayout::project(n, index, GridLayout::momentsToEx());

            auto gradPOnEx
                = layout_->deriv(Pe, index, DirectionTag<Direction::X>{}); // TODO : issue 3391

            return -gradPOnEx / nOnEx;
        }

        else if constexpr (ComponentTag::component == Component::Y)
        {
            if constexpr (Field::dimension >= 2)
            {
                auto const nOnEy = GridLayout::project(n, index, GridLayout::momentsToEy());

                auto gradPOnEy
                    = layout_->deriv(Pe, index, DirectionTag<Direction::Y>{}); // TODO : issue 3391

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
                    = layout_->deriv(Pe, index, DirectionTag<Direction::Z>{}); // TODO : issue 3391

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
            auto const& Jx    = J(Component::X);
            auto const jxOnEx = GridLayout::project(Jx, index, GridLayout::JxToEx());
            return eta_ * jxOnEx;
        }

        if constexpr (ComponentTag::component == Component::Y)
        {
            auto const& Jy    = J(Component::Y);
            auto const jyOnEy = GridLayout::project(Jy, index, GridLayout::JyToEy());
            return eta_ * jyOnEy;
        }

        if constexpr (ComponentTag::component == Component::Z)
        {
            auto const& Jz    = J(Component::Z);
            auto const jzOnEz = GridLayout::project(Jz, index, GridLayout::JzToEz());
            return eta_ * jzOnEz;
        }
    }




    template<typename VecField, typename ComponentTag>
    auto hyperresistive_(VecField const& J, MeshIndex<VecField::dimension> index,
                         ComponentTag) const
    {
        if constexpr (ComponentTag::component == Component::X)
        {
            auto const& Jx = J(Component::X);
            auto lapJx     = layout_->laplacian(Jx, index); // TODO : issue 3391
            return -nu_ * lapJx;
        }

        if constexpr (ComponentTag::component == Component::Y)
        {
            auto const& Jy = J(Component::Y);
            auto lapJy     = layout_->laplacian(Jy, index); // TODO : issue 3391
            return -nu_ * lapJy;
        }

        if constexpr (ComponentTag::component == Component::Z)
        {
            auto const& Jz = J(Component::Z);
            auto lapJz     = layout_->laplacian(Jz, index); // TODO : issue 3391
            return -nu_ * lapJz;
        }
    }
};


} // namespace PHARE::core
#endif
