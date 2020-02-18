#ifndef PHARE_OHM_H
#define PHARE_OHM_H

#include <cstddef>
#include <iostream>

#include "core/data/grid/gridlayoutdefs.h"
#include "core/data/grid/gridlayout.h"
#include "core/data/grid/gridlayout_utils.h"
#include "core/data/vecfield/vecfield_component.h"
#include "core/utilities/index/index.h"

namespace PHARE
{
namespace core
{
    template<typename GridLayout>
    class Ohm : public LayoutHolder<GridLayout>
    {
    private:
        template<typename VecField, typename ComponentTag>
        auto ideal1D_(VecField const& Ve, VecField const& B, MeshIndex<1> index, ComponentTag)
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
        auto ideal2D_(VecField const& Ve, VecField const& B, MeshIndex<2> index, ComponentTag)
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
        auto ideal3D_(VecField const& Ve, VecField const& B, MeshIndex<3> index, ComponentTag)
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




        template<typename Field, typename ComponentTag>
        auto pressure_(Field const& n, Field const& Pe, MeshIndex<Field::dimension> index,
                       ComponentTag)
        {
            if constexpr (ComponentTag::component == Component::X)
            {
                auto const nOnEx = GridLayout::project(n, index, GridLayout::momentsToEx());

                auto gradPOnEx = this->layout_->deriv(
                    Pe, index, DirectionTag<Direction::X>{}); // TODO : issue 3391

                return gradPOnEx / nOnEx;
            }

            if constexpr (ComponentTag::component == Component::Y)
            {
                if constexpr (Field::dimension >= 2)
                {
                    auto const nOnEy = GridLayout::project(n, index, GridLayout::momentsToEy());

                    auto gradPOnEy = this->layout_->deriv(
                        Pe, index, DirectionTag<Direction::Y>{}); // TODO : issue 3391

                    return gradPOnEy / nOnEy;
                }
                else
                {
                    return 0.;
                }
            }

            if constexpr (ComponentTag::component == Component::Z)
            {
                if constexpr (Field::dimension >= 3)
                {
                    auto const nOnEz = GridLayout::project(n, index, GridLayout::momentsToEz());

                    auto gradPOnEz = this->layout_->deriv(
                        Pe, index, DirectionTag<Direction::Z>{}); // TODO : issue 3391

                    return gradPOnEz / nOnEz;
                }
                else
                {
                    return 0.;
                }
            }
        }




        template<typename VecField, typename ComponentTag>
        auto resistive_(VecField const& J, MeshIndex<VecField::dimension> index, ComponentTag)
        {
            auto const eta = 1.0; // TODO : eta should comme from input file

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
        auto hyperresistive_(VecField const& J, MeshIndex<VecField::dimension> index, ComponentTag)
        {
            auto const nu = 0.001; // TODO : nu should comme from input file

            if constexpr (ComponentTag::component == Component::X)
            {
                auto const& Jx = J.getComponent(Component::X);
                auto lapJx     = this->layout_->laplacian(Jx, index); // TODO : issue 3391
                return -nu * lapJx;
            }

            if constexpr (ComponentTag::component == Component::Y)
            {
                auto const& Jy = J.getComponent(Component::Y);
                auto lapJy     = this->layout_->laplacian(Jy, index); // TODO : issue 3391
                return -nu * lapJy;
            }

            if constexpr (ComponentTag::component == Component::Z)
            {
                auto const& Jz = J.getComponent(Component::Z);
                auto lapJz     = this->layout_->laplacian(Jz, index); // TODO : issue 3391
                return -nu * lapJz;
            }
        }




        template<typename VecField, std::enable_if_t<VecField::dimension == 1, int> = 0>
        void compute_(typename VecField::field_type const& n, VecField const& Ve,
                      typename VecField::field_type const& Pe, VecField const& B, VecField const& J,
                      VecField& Enew)
        {
            auto& Ex = Enew.getComponent(Component::X);
            auto& Ey = Enew.getComponent(Component::Y);
            auto& Ez = Enew.getComponent(Component::Z);

            auto ix0 = this->layout_->physicalStartIndex(Ex, Direction::X);
            auto ix1 = this->layout_->physicalEndIndex(Ex, Direction::X);

            for (auto ix = ix0; ix <= ix1; ++ix)
            {
                Ex(ix) = ideal1D_(Ve, B, {ix}, ComponentTag<Component::X>{})
                         + pressure_(n, Pe, {ix}, ComponentTag<Component::X>{})
                         + resistive_(J, {ix}, ComponentTag<Component::X>{})
                         + hyperresistive_(J, {ix}, ComponentTag<Component::X>{});
            }

            ix0 = this->layout_->physicalStartIndex(Ey, Direction::X);
            ix1 = this->layout_->physicalEndIndex(Ey, Direction::X);

            for (auto ix = ix0; ix <= ix1; ++ix)
            {
                Ey(ix) = ideal1D_(Ve, B, {ix}, ComponentTag<Component::Y>{})
                         + pressure_(n, Pe, {ix}, ComponentTag<Component::Y>{})
                         + resistive_(J, {ix}, ComponentTag<Component::Y>{})
                         + hyperresistive_(J, {ix}, ComponentTag<Component::Y>{});
            }

            ix0 = this->layout_->physicalStartIndex(Ez, Direction::X);
            ix1 = this->layout_->physicalEndIndex(Ez, Direction::X);

            for (auto ix = ix0; ix <= ix1; ++ix)
            {
                Ez(ix) = ideal1D_(Ve, B, {ix}, ComponentTag<Component::Z>{})
                         + pressure_(n, Pe, {ix}, ComponentTag<Component::Z>{})
                         + resistive_(J, {ix}, ComponentTag<Component::Z>{})
                         + hyperresistive_(J, {ix}, ComponentTag<Component::Z>{});
            }
        }



        template<typename VecField, std::enable_if_t<VecField::dimension == 2, int> = 0>
        void compute_(typename VecField::field_type const& n, VecField const& Ve,
                      typename VecField::field_type const& Pe, VecField const& B, VecField const& J,
                      VecField& Enew)
        {
            auto& Ex = Enew.getComponent(Component::X);
            auto& Ey = Enew.getComponent(Component::Y);
            auto& Ez = Enew.getComponent(Component::Z);

            auto ix0 = this->layout_->physicalStartIndex(Ex, Direction::X);
            auto ix1 = this->layout_->physicalEndIndex(Ex, Direction::X);
            auto iy0 = this->layout_->physicalStartIndex(Ex, Direction::Y);
            auto iy1 = this->layout_->physicalEndIndex(Ex, Direction::Y);

            for (auto ix = ix0; ix <= ix1; ++ix)
            {
                for (auto iy = iy0; iy <= iy1; ++iy)
                {
                    Ex(ix, iy) = ideal2D_(Ve, B, {ix, iy}, ComponentTag<Component::X>{})
                                 + pressure_(n, Pe, {ix, iy}, ComponentTag<Component::X>{})
                                 + resistive_(J, {ix, iy}, ComponentTag<Component::X>{})
                                 + hyperresistive_(J, {ix, iy}, ComponentTag<Component::X>{});
                }
            }

            ix0 = this->layout_->physicalStartIndex(Ey, Direction::X);
            ix1 = this->layout_->physicalEndIndex(Ey, Direction::X);
            iy0 = this->layout_->physicalStartIndex(Ey, Direction::Y);
            iy1 = this->layout_->physicalEndIndex(Ey, Direction::Y);

            for (auto ix = ix0; ix <= ix1; ++ix)
            {
                for (auto iy = iy0; iy <= iy1; ++iy)
                {
                    Ey(ix, iy) = ideal2D_(Ve, B, {ix, iy}, ComponentTag<Component::Y>{})
                                 + pressure_(n, Pe, {ix, iy}, ComponentTag<Component::Y>{})
                                 + resistive_(J, {ix, iy}, ComponentTag<Component::Y>{})
                                 + hyperresistive_(J, {ix, iy}, ComponentTag<Component::Y>{});
                }
            }

            ix0 = this->layout_->physicalStartIndex(Ez, Direction::X);
            ix1 = this->layout_->physicalEndIndex(Ez, Direction::X);
            iy0 = this->layout_->physicalStartIndex(Ez, Direction::Y);
            iy1 = this->layout_->physicalEndIndex(Ez, Direction::Y);

            for (auto ix = ix0; ix <= ix1; ++ix)
            {
                for (auto iy = iy0; iy <= iy1; ++iy)
                {
                    Ez(ix, iy) = ideal2D_(Ve, B, {ix, iy}, ComponentTag<Component::Z>{})
                                 + pressure_(n, Pe, {ix, iy}, ComponentTag<Component::Z>{})
                                 + resistive_(J, {ix, iy}, ComponentTag<Component::Z>{})
                                 + hyperresistive_(J, {ix, iy}, ComponentTag<Component::Z>{});
                }
            }
        }


        template<typename VecField, std::enable_if_t<VecField::dimension == 3, int> = 0>
        void compute_(typename VecField::field_type const& n, VecField const& Ve,
                      typename VecField::field_type const& Pe, VecField const& B, VecField const& J,
                      VecField& Enew)
        {
            auto& Ex = Enew.getComponent(Component::X);
            auto& Ey = Enew.getComponent(Component::Y);
            auto& Ez = Enew.getComponent(Component::Z);

            auto ix0 = this->layout_->physicalStartIndex(Ex, Direction::X);
            auto ix1 = this->layout_->physicalEndIndex(Ex, Direction::X);
            auto iy0 = this->layout_->physicalStartIndex(Ex, Direction::Y);
            auto iy1 = this->layout_->physicalEndIndex(Ex, Direction::Y);
            auto iz0 = this->layout_->physicalStartIndex(Ex, Direction::Z);
            auto iz1 = this->layout_->physicalEndIndex(Ex, Direction::Z);

            for (auto ix = ix0; ix <= ix1; ++ix)
            {
                for (auto iy = iy0; iy <= iy1; ++iy)
                {
                    for (auto iz = iz0; iz <= iz1; ++iz)
                    {
                        Ex(ix, iy, iz)
                            = ideal3D_(Ve, B, {ix, iy, iz}, ComponentTag<Component::X>{})
                              + pressure_(n, Pe, {ix, iy, iz}, ComponentTag<Component::X>{})
                              + resistive_(J, {ix, iy, iz}, ComponentTag<Component::X>{})
                              + hyperresistive_(J, {ix, iy, iz}, ComponentTag<Component::X>{});
                    }
                }
            }

            ix0 = this->layout_->physicalStartIndex(Ey, Direction::X);
            ix1 = this->layout_->physicalEndIndex(Ey, Direction::X);
            iy0 = this->layout_->physicalStartIndex(Ey, Direction::Y);
            iy1 = this->layout_->physicalEndIndex(Ey, Direction::Y);
            iz0 = this->layout_->physicalStartIndex(Ey, Direction::Z);
            iz1 = this->layout_->physicalEndIndex(Ey, Direction::Z);

            for (auto ix = ix0; ix <= ix1; ++ix)
            {
                for (auto iy = iy0; iy <= iy1; ++iy)
                {
                    for (auto iz = iz0; iz <= iz1; ++iz)
                    {
                        Ey(ix, iy, iz)
                            = ideal3D_(Ve, B, {ix, iy, iz}, ComponentTag<Component::Y>{})
                              + pressure_(n, Pe, {ix, iy, iz}, ComponentTag<Component::Y>{})
                              + resistive_(J, {ix, iy, iz}, ComponentTag<Component::Y>{})
                              + hyperresistive_(J, {ix, iy, iz}, ComponentTag<Component::Y>{});
                    }
                }
            }

            ix0 = this->layout_->physicalStartIndex(Ez, Direction::X);
            ix1 = this->layout_->physicalEndIndex(Ez, Direction::X);
            iy0 = this->layout_->physicalStartIndex(Ez, Direction::Y);
            iy1 = this->layout_->physicalEndIndex(Ez, Direction::Y);
            iz0 = this->layout_->physicalStartIndex(Ez, Direction::Z);
            iz1 = this->layout_->physicalEndIndex(Ez, Direction::Z);

            for (auto ix = ix0; ix <= ix1; ++ix)
            {
                for (auto iy = iy0; iy <= iy1; ++iy)
                {
                    for (auto iz = iz0; iz <= iz1; ++iz)
                    {
                        Ez(ix, iy, iz)
                            = ideal3D_(Ve, B, {ix, iy, iz}, ComponentTag<Component::Z>{})
                              + pressure_(n, Pe, {ix, iy, iz}, ComponentTag<Component::Z>{})
                              + resistive_(J, {ix, iy, iz}, ComponentTag<Component::Z>{})
                              + hyperresistive_(J, {ix, iy, iz}, ComponentTag<Component::Z>{});
                    }
                }
            }
        }

    public:
        template<typename VecField>
        void operator()(typename VecField::field_type const& n, VecField const& Ve,
                        typename VecField::field_type const& Pe, VecField const& B,
                        VecField const& J, VecField& Enew)
        {
            if (!this->hasLayout())
            {
                throw std::runtime_error(
                    "Error - Ohm - GridLayout not set, cannot proceed to calculate ohm()");
            }

            compute_(n, Ve, Pe, B, J, Enew);
        }
    };
} // namespace core
} // namespace PHARE



#endif
