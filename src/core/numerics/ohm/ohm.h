#ifndef PHARE_CORE_NUMERICS_OHM_OHM_H
#define PHARE_CORE_NUMERICS_OHM_OHM_H

#include <cstddef>
#include <iostream>

#include "data/grid/gridlayout.h"
#include "data/grid/gridlayoutdefs.h"
#include "data/vecfield/vecfield_component.h"
#include "utilities/index/index.h"

namespace PHARE
{
namespace core
{
    /** @brief Actual implementation of the ohm law
     *
     * This implementation is used by the Ohm object to calculate the electric field
     * It is templated by the layout and the dimension so specialized template actually do the job
     *
     */
    template<typename GridLayout, std::size_t dim>
    class OhmImpl
    {
    };


    /** @brief Base class for OhmImpl specialized classes
     * to factorize the code that deals with the GridLayout.
     *
     * The purpose of this class is to store the GridLayout
     * pointer of Ohm, and propose the interface to check
     * whether it's been set or not.
     *
     */
    template<typename GridLayout>
    class OhmImplInternals
    {
    protected:
        GridLayout* layout_{nullptr};

    public:
        /**
         * @brief hasLayoutSet returns true if the Layout has been set
         */
        bool hasLayoutSet() const { return (layout_ == nullptr) ? false : true; }


        /**
         * @brief setLayout is used to give Ohm a pointer to a gridlayout
         */
        void setLayout(GridLayout* layout)
        {
            if (layout_ != nullptr)
            {
                throw std::runtime_error(
                    "Error - Ohm - cannot set layout_ because it is already set");
            }
            else
            {
                layout_ = layout;
            }
        }
    };


    /** @brief 1D specialization of the ohm's law implementation
     */
    template<typename GridLayout>
    class OhmImpl<GridLayout, 1> : public OhmImplInternals<GridLayout>
    {
        static_assert(GridLayout::dimension == 1, "Error: Passed non-1D GridLayout to 1D OhmImpl");

    public:
        template<typename Field, typename VecField>
        void operator()(Field const& n, VecField const& Ve, Field const& Pe, VecField const& B,
                        VecField const& J, VecField& Enew)
        {
            auto& Ex = Enew.getComponent(Component::X);
            auto& Ey = Enew.getComponent(Component::Y);
            auto& Ez = Enew.getComponent(Component::Z);

            auto ix0 = this->layout_->physicalStartIndex(Ex, Direction::X);
            auto ix1 = this->layout_->physicalEndIndex(Ex, Direction::X);

            for (auto ix = ix0; ix <= ix1; ++ix)
            {
                Ex(ix) = ideal_(Ve, B, {ix}, ComponentTag<Component::X>{})
                         + pressure_(n, Pe, {ix}, ComponentTag<Component::X>{})
                         + resistive_(J, {ix}, ComponentTag<Component::X>{})
                         + hyperresistive_(J, {ix}, ComponentTag<Component::X>{});
            }

            ix0 = this->layout_->physicalStartIndex(Ey, Direction::X);
            ix1 = this->layout_->physicalEndIndex(Ey, Direction::X);

            for (auto ix = ix0; ix <= ix1; ++ix)
            {
                Ey(ix) = ideal_(Ve, B, {ix}, ComponentTag<Component::Y>{})
                         + pressure_(n, Pe, {ix}, ComponentTag<Component::Y>{})
                         + resistive_(J, {ix}, ComponentTag<Component::Y>{})
                         + hyperresistive_(J, {ix}, ComponentTag<Component::Y>{});
            }

            ix0 = this->layout_->physicalStartIndex(Ez, Direction::X);
            ix1 = this->layout_->physicalEndIndex(Ez, Direction::X);

            for (auto ix = ix0; ix <= ix1; ++ix)
            {
                Ez(ix) = ideal_(Ve, B, {ix}, ComponentTag<Component::Z>{})
                         + pressure_(n, Pe, {ix}, ComponentTag<Component::Z>{})
                         + resistive_(J, {ix}, ComponentTag<Component::Z>{})
                         + hyperresistive_(J, {ix}, ComponentTag<Component::Z>{});
            }
        }


    private:
        template<typename VecField, typename ComponentTag>
        auto ideal_(VecField const& Ve, VecField const& B, MeshIndex<1> index, ComponentTag)
        {
            if constexpr (ComponentTag::component == Component::X)
            {
                auto const& Vy = Ve.getComponent(Component::Y);
                auto const& Vz = Ve.getComponent(Component::Z);

                auto const& By = B.getComponent(Component::Y);
                auto const& Bz = B.getComponent(Component::Z);

                auto vyOnEx = 0.0;
                for (auto const& wp : GridLayout::momentsToEx())
                {
                    vyOnEx += wp.coef * Vy(index[0] + wp.indexes[0]);
                }

                auto vzOnEx = 0.0;
                for (auto const& wp : GridLayout::momentsToEx())
                {
                    vzOnEx += wp.coef * Vz(index[0] + wp.indexes[0]);
                }

                auto byOnEx = 0.0;
                for (auto const& wp : GridLayout::ByToEx())
                {
                    byOnEx += wp.coef * By(index[0] + wp.indexes[0]);
                }

                auto bzOnEx = 0.0;
                for (auto const& wp : GridLayout::BzToEx())
                {
                    bzOnEx += wp.coef * Bz(index[0] + wp.indexes[0]);
                }

                return -vyOnEx * bzOnEx + vzOnEx * byOnEx;
            }

            if constexpr (ComponentTag::component == Component::Y)
            {
                auto const& Vx = Ve.getComponent(Component::X);
                auto const& Vz = Ve.getComponent(Component::Z);

                auto const& Bx = B.getComponent(Component::X);
                auto const& Bz = B.getComponent(Component::Z);

                auto vxOnEy = 0.0;
                for (auto const& wp : GridLayout::momentsToEy())
                {
                    vxOnEy += wp.coef * Vx(index[0] + wp.indexes[0]);
                }

                auto vzOnEy = 0.0;
                for (auto const& wp : GridLayout::momentsToEy())
                {
                    vzOnEy += wp.coef * Vz(index[0] + wp.indexes[0]);
                }

                auto bxOnEy = 0.0;
                for (auto const& wp : GridLayout::BxToEy())
                {
                    bxOnEy += wp.coef * Bx(index[0] + wp.indexes[0]);
                }

                auto bzOnEy = 0.0;
                for (auto const& wp : GridLayout::BzToEy())
                {
                    bzOnEy += wp.coef * Bz(index[0] + wp.indexes[0]);
                }

                return -vzOnEy * bxOnEy + vxOnEy * bzOnEy;
            }

            if constexpr (ComponentTag::component == Component::Z)
            {
                auto const& Vx = Ve.getComponent(Component::X);
                auto const& Vy = Ve.getComponent(Component::Y);

                auto const& Bx = B.getComponent(Component::X);
                auto const& By = B.getComponent(Component::Y);

                auto vxOnEz = 0.0;
                for (auto const& wp : GridLayout::momentsToEz())
                {
                    vxOnEz += wp.coef * Vx(index[0] + wp.indexes[0]);
                }

                auto vyOnEz = 0.0;
                for (auto const& wp : GridLayout::momentsToEz())
                {
                    vyOnEz += wp.coef * Vy(index[0] + wp.indexes[0]);
                }

                auto bxOnEz = 0.0;
                for (auto const& wp : GridLayout::BxToEz())
                {
                    bxOnEz += wp.coef * Bx(index[0] + wp.indexes[0]);
                }

                auto byOnEz = 0.0;
                for (auto const& wp : GridLayout::ByToEz())
                {
                    byOnEz += wp.coef * By(index[0] + wp.indexes[0]);
                }

                return -vxOnEz * byOnEz + vyOnEz * bxOnEz;
            }
        }

        template<typename Field, typename ComponentTag>
        auto pressure_(Field const& n, Field const& Pe, MeshIndex<1> index, ComponentTag)
        {
            if constexpr (ComponentTag::component == Component::X)
            {
                auto nOnEx = 0.0;
                for (auto const& wp : GridLayout::momentsToEx())
                {
                    nOnEx += wp.coef * n(index[0] + wp.indexes[0]);
                }

                auto gradPOnEx = this->layout_->deriv(
                    Pe, index, DirectionTag<Direction::X>{}); // TODO : issue 3391

                return gradPOnEx / nOnEx;
            }

            if constexpr (ComponentTag::component == Component::Y)
            {
                return 0.0;
            }

            if constexpr (ComponentTag::component == Component::Z)
            {
                return 0.0;
            }
        }

        template<typename VecField, typename ComponentTag>
        auto resistive_(VecField const& J, MeshIndex<1> index, ComponentTag)
        {
            auto const eta = 1.0; // TODO : eta should comme from input file

            if constexpr (ComponentTag::component == Component::X)
            {
                auto const& Jx = J.getComponent(Component::X);

                auto jxOnEx = 0.0;
                for (auto const& wp : GridLayout::JxToEx())
                {
                    jxOnEx += wp.coef * Jx(index[0] + wp.indexes[0]);
                }

                return eta * jxOnEx;
            }

            if constexpr (ComponentTag::component == Component::Y)
            {
                auto const& Jy = J.getComponent(Component::Y);

                auto jyOnEy = 0.0;
                for (auto const& wp : GridLayout::JyToEy())
                {
                    jyOnEy += wp.coef * Jy(index[0] + wp.indexes[0]);
                }

                return eta * jyOnEy;
            }

            if constexpr (ComponentTag::component == Component::Z)
            {
                auto const& Jz = J.getComponent(Component::Z);

                auto jzOnEz = 0.0;
                for (auto const& wp : GridLayout::JzToEz())
                {
                    jzOnEz += wp.coef * Jz(index[0] + wp.indexes[0]);
                }

                return eta * jzOnEz;
            }
        }

        template<typename VecField, typename ComponentTag>
        auto hyperresistive_(VecField const& J, MeshIndex<1> index, ComponentTag)
        {
            auto const nu = 0.001; // TODO : nu should comme from input file

            if constexpr (ComponentTag::component == Component::X)
            {
                auto const& Jx = J.getComponent(Component::X);

                auto lapJx = this->layout_->laplacian(Jx, index); // TODO : issue 3391

                return -nu * lapJx;
            }

            if constexpr (ComponentTag::component == Component::Y)
            {
                auto const& Jy = J.getComponent(Component::Y);

                auto lapJy = this->layout_->laplacian(Jy, index); // TODO : issue 3391

                return -nu * lapJy;
            }

            if constexpr (ComponentTag::component == Component::Z)
            {
                auto const& Jz = J.getComponent(Component::Z);

                auto lapJz = this->layout_->laplacian(Jz, index); // TODO : issue 3391

                return -nu * lapJz;
            }
        }
    };



    template<typename GridLayout>
    class OhmImpl<GridLayout, 2> : public OhmImplInternals<GridLayout>
    {
        static_assert(GridLayout::dimension == 2, "Error: Passed non-2D GridLayout to 2D OhmImpl");

    public:
        template<typename Field, typename VecField>
        void operator()(Field const& n, VecField const& Ve, Field const& Pe, VecField const& B,
                        VecField const& J, VecField& Enew)
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
                    Ex(ix, iy) = ideal_(Ve, B, {ix, iy}, ComponentTag<Component::X>{})
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
                    Ey(ix, iy) = ideal_(Ve, B, {ix, iy}, ComponentTag<Component::Y>{})
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
                    Ez(ix, iy) = ideal_(Ve, B, {ix, iy}, ComponentTag<Component::Z>{})
                                 + pressure_(n, Pe, {ix, iy}, ComponentTag<Component::Z>{})
                                 + resistive_(J, {ix, iy}, ComponentTag<Component::Z>{})
                                 + hyperresistive_(J, {ix, iy}, ComponentTag<Component::Z>{});
                }
            }
        }

    private:
        template<typename VecField, typename ComponentTag>
        auto ideal_(VecField const& Ve, VecField const& B, MeshIndex<2> index, ComponentTag)
        {
            if constexpr (ComponentTag::component == Component::X)
            {
                auto const& Vy = Ve.getComponent(Component::Y);
                auto const& Vz = Ve.getComponent(Component::Z);

                auto const& By = B.getComponent(Component::Y);
                auto const& Bz = B.getComponent(Component::Z);

                auto vyOnEx = 0.0;
                for (auto const& wp : GridLayout::momentsToEx())
                {
                    vyOnEx += wp.coef * Vy(index[0] + wp.indexes[0], index[1] + wp.indexes[1]);
                }

                // vyOnEx = linearCombination(Vy, index, GridLayout::momentsToEx());

                auto vzOnEx = 0.0;
                for (auto const& wp : GridLayout::momentsToEx())
                {
                    vzOnEx += wp.coef * Vz(index[0] + wp.indexes[0], index[1] + wp.indexes[1]);
                }

                auto byOnEx = 0.0;
                for (auto const& wp : GridLayout::ByToEx())
                {
                    byOnEx += wp.coef * By(index[0] + wp.indexes[0], index[1] + wp.indexes[1]);
                }

                auto bzOnEx = 0.0;
                for (auto const& wp : GridLayout::BzToEx())
                {
                    bzOnEx += wp.coef * Bz(index[0] + wp.indexes[0], index[1] + wp.indexes[1]);
                }

                return -vyOnEx * bzOnEx + vzOnEx * byOnEx;
            }

            if constexpr (ComponentTag::component == Component::Y)
            {
                auto const& Vx = Ve.getComponent(Component::X);
                auto const& Vz = Ve.getComponent(Component::Z);

                auto const& Bx = B.getComponent(Component::X);
                auto const& Bz = B.getComponent(Component::Z);

                auto vxOnEy = 0.0;
                for (auto const& wp : GridLayout::momentsToEy())
                {
                    vxOnEy += wp.coef * Vx(index[0] + wp.indexes[0], index[1] + wp.indexes[1]);
                }

                auto vzOnEy = 0.0;
                for (auto const& wp : GridLayout::momentsToEy())
                {
                    vzOnEy += wp.coef * Vz(index[0] + wp.indexes[0], index[1] + wp.indexes[1]);
                }

                auto bxOnEy = 0.0;
                for (auto const& wp : GridLayout::BxToEy())
                {
                    bxOnEy += wp.coef * Bx(index[0] + wp.indexes[0], index[1] + wp.indexes[1]);
                }

                auto bzOnEy = 0.0;
                for (auto const& wp : GridLayout::BzToEy())
                {
                    bzOnEy += wp.coef * Bz(index[0] + wp.indexes[0], index[1] + wp.indexes[1]);
                }

                return -vzOnEy * bxOnEy + vxOnEy * bzOnEy;
            }

            if constexpr (ComponentTag::component == Component::Z)
            {
                auto const& Vx = Ve.getComponent(Component::X);
                auto const& Vy = Ve.getComponent(Component::Y);

                auto const& Bx = B.getComponent(Component::X);
                auto const& By = B.getComponent(Component::Y);

                auto vxOnEz = 0.0;
                for (auto const& wp : GridLayout::momentsToEz())
                {
                    vxOnEz += wp.coef * Vx(index[0] + wp.indexes[0], index[1] + wp.indexes[1]);
                }

                auto vyOnEz = 0.0;
                for (auto const& wp : GridLayout::momentsToEz())
                {
                    vyOnEz += wp.coef * Vy(index[0] + wp.indexes[0], index[1] + wp.indexes[1]);
                }

                auto bxOnEz = 0.0;
                for (auto const& wp : GridLayout::BxToEz())
                {
                    bxOnEz += wp.coef * Bx(index[0] + wp.indexes[0], index[1] + wp.indexes[1]);
                }

                auto byOnEz = 0.0;
                for (auto const& wp : GridLayout::ByToEz())
                {
                    byOnEz += wp.coef * By(index[0] + wp.indexes[0], index[1] + wp.indexes[1]);
                }

                return -vxOnEz * byOnEz + vyOnEz * bxOnEz;
            }
        }

        template<typename Field, typename ComponentTag>
        auto pressure_(Field const& n, Field const& Pe, MeshIndex<2> index, ComponentTag)
        {
            if constexpr (ComponentTag::component == Component::X)
            {
                auto nOnEx = 0.0;
                for (auto const& wp : GridLayout::momentsToEx())
                {
                    nOnEx += wp.coef * n(index[0] + wp.indexes[0], index[1] + wp.indexes[1]);
                }

                auto gradPOnEx = this->layout_->deriv(
                    Pe, index, DirectionTag<Direction::X>{}); // TODO : issue 3391

                return gradPOnEx / nOnEx;
            }

            if constexpr (ComponentTag::component == Component::Y)
            {
                auto nOnEy = 0.0;
                for (auto const& wp : GridLayout::momentsToEy())
                {
                    nOnEy += wp.coef * n(index[0] + wp.indexes[0], index[1] + wp.indexes[1]);
                }

                auto gradPOnEy = this->layout_->deriv(
                    Pe, index, DirectionTag<Direction::Y>{}); // TODO : issue 3391

                return gradPOnEy / nOnEy;
            }

            if constexpr (ComponentTag::component == Component::Z)
            {
                return 0.0;
            }
        }



        template<typename VecField, typename ComponentTag>
        auto resistive_(VecField const& J, MeshIndex<2> index, ComponentTag)
        {
            auto const eta = 1.0; // TODO : eta should comme from input file

            if constexpr (ComponentTag::component == Component::X)
            {
                auto const& Jx = J.getComponent(Component::X);

                auto jxOnEx = 0.0;
                for (auto const& wp : GridLayout::JxToEx())
                {
                    jxOnEx += wp.coef * Jx(index[0] + wp.indexes[0], index[1] + wp.indexes[1]);
                }

                return eta * jxOnEx;
            }

            if constexpr (ComponentTag::component == Component::Y)
            {
                auto const& Jy = J.getComponent(Component::Y);

                auto jyOnEy = 0.0;
                for (auto const& wp : GridLayout::JyToEy())
                {
                    jyOnEy += wp.coef * Jy(index[0] + wp.indexes[0], index[1] + wp.indexes[1]);
                }

                return eta * jyOnEy;
            }

            if constexpr (ComponentTag::component == Component::Z)
            {
                auto const& Jz = J.getComponent(Component::Z);

                auto jzOnEz = 0.0;
                for (auto const& wp : GridLayout::JzToEz())
                {
                    jzOnEz += wp.coef * Jz(index[0] + wp.indexes[0], index[1] + wp.indexes[1]);
                }

                return eta * jzOnEz;
            }
        }


        template<typename VecField, typename ComponentTag>
        auto hyperresistive_(VecField const& J, MeshIndex<2> index, ComponentTag)
        {
            auto const nu = 0.001; // TODO : nu should comme from input file

            if constexpr (ComponentTag::component == Component::X)
            {
                auto const& Jx = J.getComponent(Component::X);

                auto lapJx = this->layout_->laplacian(Jx, index); // TODO : issue 3391

                return -nu * lapJx;
            }

            if constexpr (ComponentTag::component == Component::Y)
            {
                auto const& Jy = J.getComponent(Component::Y);

                auto lapJy = this->layout_->laplacian(Jy, index); // TODO : issue 3391

                return -nu * lapJy;
            }

            if constexpr (ComponentTag::component == Component::Z)
            {
                auto const& Jz = J.getComponent(Component::Z);

                auto lapJz = this->layout_->laplacian(Jz, index); // TODO : issue 3391

                return -nu * lapJz;
            }
        }
    };



    template<typename GridLayout>
    class OhmImpl<GridLayout, 3> : public OhmImplInternals<GridLayout>
    {
        static_assert(GridLayout::dimension == 3, "Error: Passed non-3D GridLayout to 3D OhmImpl");

    public:
        template<typename Field, typename VecField>
        void operator()(Field const& n, VecField const& Ve, Field const& Pe, VecField const& B,
                        VecField const& J, VecField& Enew)
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
                            = ideal_(Ve, B, {ix, iy, iz}, ComponentTag<Component::X>{})
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
                            = ideal_(Ve, B, {ix, iy, iz}, ComponentTag<Component::Y>{})
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
                            = ideal_(Ve, B, {ix, iy, iz}, ComponentTag<Component::Z>{})
                              + pressure_(n, Pe, {ix, iy, iz}, ComponentTag<Component::Z>{})
                              + resistive_(J, {ix, iy, iz}, ComponentTag<Component::Z>{})
                              + hyperresistive_(J, {ix, iy, iz}, ComponentTag<Component::Z>{});
                    }
                }
            }
        }

    private:
        template<typename VecField, typename ComponentTag>
        auto ideal_(VecField const& Ve, VecField const& B, MeshIndex<3> index, ComponentTag)
        {
            if constexpr (ComponentTag::component == Component::X)
            {
                auto const& Vy = Ve.getComponent(Component::Y);
                auto const& Vz = Ve.getComponent(Component::Z);

                auto const& By = B.getComponent(Component::Y);
                auto const& Bz = B.getComponent(Component::Z);

                auto vyOnEx = 0.0;
                for (auto const& wp : GridLayout::momentsToEx())
                {
                    vyOnEx += wp.coef
                              * Vy(index[0] + wp.indexes[0], index[1] + wp.indexes[1],
                                   index[2] + wp.indexes[2]);
                }

                auto vzOnEx = 0.0;
                for (auto const& wp : GridLayout::momentsToEx())
                {
                    vzOnEx += wp.coef
                              * Vz(index[0] + wp.indexes[0], index[1] + wp.indexes[1],
                                   index[1] + wp.indexes[2]);
                }

                auto byOnEx = 0.0;
                for (auto const& wp : GridLayout::ByToEx())
                {
                    byOnEx += wp.coef
                              * By(index[0] + wp.indexes[0], index[1] + wp.indexes[1],
                                   index[2] + wp.indexes[2]);
                }

                auto bzOnEx = 0.0;
                for (auto const& wp : GridLayout::BzToEx())
                {
                    bzOnEx += wp.coef
                              * Bz(index[0] + wp.indexes[0], index[1] + wp.indexes[1],
                                   index[2] + wp.indexes[2]);
                }

                return -vyOnEx * bzOnEx + vzOnEx * byOnEx;
            }

            if constexpr (ComponentTag::component == Component::Y)
            {
                auto const& Vx = Ve.getComponent(Component::X);
                auto const& Vz = Ve.getComponent(Component::Z);

                auto const& Bx = B.getComponent(Component::X);
                auto const& Bz = B.getComponent(Component::Z);

                auto vxOnEy = 0.0;
                for (auto const& wp : GridLayout::momentsToEy())
                {
                    vxOnEy += wp.coef
                              * Vx(index[0] + wp.indexes[0], index[1] + wp.indexes[1],
                                   index[2] + wp.indexes[2]);
                }

                auto vzOnEy = 0.0;
                for (auto const& wp : GridLayout::momentsToEy())
                {
                    vzOnEy += wp.coef
                              * Vz(index[0] + wp.indexes[0], index[1] + wp.indexes[1],
                                   index[2] + wp.indexes[2]);
                }

                auto bxOnEy = 0.0;
                for (auto const& wp : GridLayout::BxToEy())
                {
                    bxOnEy += wp.coef
                              * Bx(index[0] + wp.indexes[0], index[1] + wp.indexes[1],
                                   index[2] + wp.indexes[2]);
                }

                auto bzOnEy = 0.0;
                for (auto const& wp : GridLayout::BzToEy())
                {
                    bzOnEy += wp.coef
                              * Bz(index[0] + wp.indexes[0], index[1] + wp.indexes[1],
                                   index[2] + wp.indexes[2]);
                }

                return -vzOnEy * bxOnEy + vxOnEy * bzOnEy;
            }

            if constexpr (ComponentTag::component == Component::Z)
            {
                auto const& Vx = Ve.getComponent(Component::X);
                auto const& Vy = Ve.getComponent(Component::Y);

                auto const& Bx = B.getComponent(Component::X);
                auto const& By = B.getComponent(Component::Y);

                auto vxOnEz = 0.0;
                for (auto const& wp : GridLayout::momentsToEz())
                {
                    vxOnEz += wp.coef
                              * Vx(index[0] + wp.indexes[0], index[1] + wp.indexes[1],
                                   index[2] + wp.indexes[2]);
                }

                auto vyOnEz = 0.0;
                for (auto const& wp : GridLayout::momentsToEz())
                {
                    vyOnEz += wp.coef
                              * Vy(index[0] + wp.indexes[0], index[1] + wp.indexes[1],
                                   index[2] + wp.indexes[2]);
                }

                auto bxOnEz = 0.0;
                for (auto const& wp : GridLayout::BxToEz())
                {
                    bxOnEz += wp.coef
                              * Bx(index[0] + wp.indexes[0], index[1] + wp.indexes[1],
                                   index[2] + wp.indexes[2]);
                }

                auto byOnEz = 0.0;
                for (auto const& wp : GridLayout::ByToEz())
                {
                    byOnEz += wp.coef
                              * By(index[0] + wp.indexes[0], index[1] + wp.indexes[1],
                                   index[2] + wp.indexes[2]);
                }

                return -vxOnEz * byOnEz + vyOnEz * bxOnEz;
            }
        }

        template<typename Field, typename ComponentTag>
        auto pressure_(Field const& n, Field const& Pe, MeshIndex<3> index, ComponentTag)
        {
            if constexpr (ComponentTag::component == Component::X)
            {
                auto nOnEx = 0.0;
                for (auto const& wp : GridLayout::momentsToEx())
                {
                    nOnEx += wp.coef
                             * n(index[0] + wp.indexes[0], index[1] + wp.indexes[1],
                                 index[2] + wp.indexes[2]);
                }

                auto gradPOnEx = this->layout_->deriv(
                    Pe, index, DirectionTag<Direction::X>{}); // TODO : issue 3391

                return gradPOnEx / nOnEx;
            }

            if constexpr (ComponentTag::component == Component::Y)
            {
                auto nOnEy = 0.0;
                for (auto const& wp : GridLayout::momentsToEy())
                {
                    nOnEy += wp.coef
                             * n(index[0] + wp.indexes[0], index[1] + wp.indexes[1],
                                 index[2] + wp.indexes[2]);
                }

                auto gradPOnEy = this->layout_->deriv(
                    Pe, index, DirectionTag<Direction::Y>{}); // TODO : issue 3391

                return gradPOnEy / nOnEy;
            }

            if constexpr (ComponentTag::component == Component::Z)
            {
                auto nOnEz = 0.0;
                for (auto const& wp : GridLayout::momentsToEz())
                {
                    nOnEz += wp.coef
                             * n(index[0] + wp.indexes[0], index[1] + wp.indexes[1],
                                 index[2] + wp.indexes[2]);
                }

                auto gradPOnEz = this->layout_->deriv(
                    Pe, index, DirectionTag<Direction::Z>{}); // TODO : issue 3391

                return gradPOnEz / nOnEz;
            }
        }


        template<typename VecField, typename ComponentTag>
        auto resistive_(VecField const& J, MeshIndex<3> index, ComponentTag)
        {
            auto const eta = 1.0; // TODO : eta should comme from input file

            if constexpr (ComponentTag::component == Component::X)
            {
                auto const& Jx = J.getComponent(Component::X);

                auto jxOnEx = 0.0;
                for (auto const& wp : GridLayout::JxToEx())
                {
                    jxOnEx += wp.coef
                              * Jx(index[0] + wp.indexes[0], index[1] + wp.indexes[1],
                                   index[2] + wp.indexes[2]);
                }

                return eta * jxOnEx;
            }

            if constexpr (ComponentTag::component == Component::Y)
            {
                auto const& Jy = J.getComponent(Component::Y);

                auto jyOnEy = 0.0;
                for (auto const& wp : GridLayout::JyToEy())
                {
                    jyOnEy += wp.coef
                              * Jy(index[0] + wp.indexes[0], index[1] + wp.indexes[1],
                                   index[2] + wp.indexes[2]);
                }

                return eta * jyOnEy;
            }

            if constexpr (ComponentTag::component == Component::Z)
            {
                auto const& Jz = J.getComponent(Component::Z);

                auto jzOnEz = 0.0;
                for (auto const& wp : GridLayout::JzToEz())
                {
                    jzOnEz += wp.coef
                              * Jz(index[0] + wp.indexes[0], index[2] + wp.indexes[1],
                                   index[2] + wp.indexes[2]);
                }

                return eta * jzOnEz;
            }
        }


        template<typename VecField, typename ComponentTag>
        auto hyperresistive_(VecField const& J, MeshIndex<3> index, ComponentTag)
        {
            auto const nu = 0.001; // TODO : nu should comme from input file

            if constexpr (ComponentTag::component == Component::X)
            {
                auto const& Jx = J.getComponent(Component::X);

                auto lapJx = this->layout_->laplacian(Jx, index); // TODO : issue 3391

                return -nu * lapJx;
            }

            if constexpr (ComponentTag::component == Component::Y)
            {
                auto const& Jy = J.getComponent(Component::Y);

                auto lapJy = this->layout_->laplacian(Jy, index); // TODO : issue 3391

                return -nu * lapJy;
            }

            if constexpr (ComponentTag::component == Component::Z)
            {
                auto const& Jz = J.getComponent(Component::Z);

                auto lapJz = this->layout_->laplacian(Jz, index); // TODO : issue 3391

                return -nu * lapJz;
            }
        }
    };




    template<typename GridLayout>
    class Ohm
    {
    private:
        OhmImpl<GridLayout, GridLayout::dimension> impl_;

    public:
        template<typename Field, typename VecField>
        void operator()(Field const& n, VecField const& Ve, Field const& Pe, VecField const& B,
                        VecField const& J, VecField& Enew)
        {
            if (!impl_.hasLayoutSet())
            {
                throw std::runtime_error(
                    "Error - Ohm - GridLayout not set, cannot proceed to calculate ohm()");
            }

            impl_(n, Ve, Pe, B, J, Enew);
        }


        void setLayout(GridLayout* layout) { impl_.setLayout(layout); }
    };
} // namespace core
} // namespace PHARE



#endif
