#ifndef PHARE_CORE_NUMERICS_AMPERE_AMPERE_H
#define PHARE_CORE_NUMERICS_AMPERE_AMPERE_H

#include <cstddef>
#include <iostream>

#include "core/data/grid/gridlayoutdefs.h"
#include "core/data/grid/gridlayout_utils.h"
#include "core/data/vecfield/vecfield_component.h"
#include "core/utilities/index/index.h"

namespace PHARE
{
namespace core
{
    template<typename GridLayout>
    class Ampere : public LayoutHolder<GridLayout>
    {
    private:
        template<typename VecField, std::size_t dim, std::enable_if_t<dim == 1, int> = 0>
        void operator()(VecField const& B, VecField& J)
        {
            // auto &Jx = J.getComponent(Component::X); // =  0
            auto& Jy = J.getComponent(Component::Y); // = -dxBz
            auto& Jz = J.getComponent(Component::Z); // =  dxBy

            // auto const &Bx = B.getComponent(Component::X);
            auto const& By = B.getComponent(Component::Y);
            auto const& Bz = B.getComponent(Component::Z);

            // TODO Direction should not be in gridlayoutdef but in utilities somehow
            // TODO 1st arg of physicalStartIndex( could be QtyCentering::primal ?
            auto start = this->layout_->physicalStartIndex(Jy, Direction::X);
            auto end   = this->layout_->physicalEndIndex(Jy, Direction::X);

            for (auto ix = start; ix <= end; ++ix)
            {
                Jy(ix) = -this->layout_->deriv(Bz, {ix}, DirectionTag<Direction::X>{});
            }

            start = this->layout_->physicalStartIndex(Jz, Direction::X);
            end   = this->layout_->physicalEndIndex(Jz, Direction::X);

            for (auto ix = start; ix <= end; ++ix)
            {
                Jz(ix) = this->layout_->deriv(By, {ix}, DirectionTag<Direction::X>{});
            }
        }




        template<typename VecField, std::size_t dim, std::enable_if_t<dim == 2, int> = 0>
        void operator()(VecField const& B, VecField& J)
        {
            auto& Jx = J.getComponent(Component::X); // =  dyBz
            auto& Jy = J.getComponent(Component::Y); // = -dxBz
            auto& Jz = J.getComponent(Component::Z); // =  dxBy - dyBx

            auto const& Bx = B.getComponent(Component::X);
            auto const& By = B.getComponent(Component::Y);
            auto const& Bz = B.getComponent(Component::Z);

            auto psi_X = this->layout_->physicalStartIndex(Jx, Direction::X);
            auto pei_X = this->layout_->physicalEndIndex(Jx, Direction::X);
            auto psi_Y = this->layout_->physicalStartIndex(Jx, Direction::Y);
            auto pei_Y = this->layout_->physicalEndIndex(Jx, Direction::Y);

            for (auto ix = psi_X; ix <= pei_X; ++ix)
            {
                for (auto iy = psi_Y; iy <= pei_Y; ++iy)
                {
                    Jx(ix, iy) = this->layout_->deriv(Bz, {ix, iy}, DirectionTag<Direction::Y>{});
                }
            }

            psi_X = this->layout_->physicalStartIndex(Jy, Direction::X);
            pei_X = this->layout_->physicalEndIndex(Jy, Direction::X);
            psi_Y = this->layout_->physicalStartIndex(Jy, Direction::Y);
            pei_Y = this->layout_->physicalEndIndex(Jy, Direction::Y);

            for (auto ix = psi_X; ix <= pei_X; ++ix)
            {
                for (auto iy = psi_Y; iy <= pei_Y; ++iy)
                {
                    Jy(ix, iy) = -this->layout_->deriv(Bz, {ix, iy}, DirectionTag<Direction::X>{});
                }
            }

            psi_X = this->layout_->physicalStartIndex(Jz, Direction::X);
            pei_X = this->layout_->physicalEndIndex(Jz, Direction::X);
            psi_Y = this->layout_->physicalStartIndex(Jz, Direction::Y);
            pei_Y = this->layout_->physicalEndIndex(Jz, Direction::Y);

            for (uint32 ix = psi_X; ix <= pei_X; ++ix)
            {
                for (uint32 iy = psi_Y; iy <= pei_Y; ++iy)
                {
                    Jz(ix, iy) = this->layout_->deriv(By, {ix, iy}, DirectionTag<Direction::X>{})
                                 - this->layout_->deriv(Bx, {ix, iy}, DirectionTag<Direction::Y>{});
                }
            }
        }




        template<typename VecField, std::size_t dim, std::enable_if_t<dim == 3, int> = 0>
        void operator()(VecField const& B, VecField& J)
        {
            auto& Jx = J.getComponent(Component::X); // =  dyBz - dzBx
            auto& Jy = J.getComponent(Component::Y); // =  dzBx - dxBz
            auto& Jz = J.getComponent(Component::Z); // =  dxBy - dyBx

            auto const& Bx = B.getComponent(Component::X);
            auto const& By = B.getComponent(Component::Y);
            auto const& Bz = B.getComponent(Component::Z);

            auto psi_X = this->layout_->physicalStartIndex(Jx, Direction::X);
            auto pei_X = this->layout_->physicalEndIndex(Jx, Direction::X);
            auto psi_Y = this->layout_->physicalStartIndex(Jx, Direction::Y);
            auto pei_Y = this->layout_->physicalEndIndex(Jx, Direction::Y);
            auto psi_Z = this->layout_->physicalStartIndex(Jx, Direction::Z);
            auto pei_Z = this->layout_->physicalEndIndex(Jx, Direction::Z);

            for (auto ix = psi_X; ix <= pei_X; ++ix)
            {
                for (auto iy = psi_Y; iy <= pei_Y; ++iy)
                {
                    for (auto iz = psi_Z; iz <= pei_Z; ++iz)
                    {
                        Jx(ix, iy, iz)
                            = this->layout_->deriv(Bz, {ix, iy, iz}, DirectionTag<Direction::Y>{})
                              - this->layout_->deriv(By, {ix, iy, iz},
                                                     DirectionTag<Direction::Z>{});
                    }
                }
            }

            psi_X = this->layout_->physicalStartIndex(Jy, Direction::X);
            pei_X = this->layout_->physicalEndIndex(Jy, Direction::X);
            psi_Y = this->layout_->physicalStartIndex(Jy, Direction::Y);
            pei_Y = this->layout_->physicalEndIndex(Jy, Direction::Y);
            psi_Z = this->layout_->physicalStartIndex(Jy, Direction::Z);
            pei_Z = this->layout_->physicalEndIndex(Jy, Direction::Z);

            for (auto ix = psi_X; ix <= pei_X; ++ix)
            {
                for (auto iy = psi_Y; iy <= pei_Y; ++iy)
                {
                    for (auto iz = psi_Z; iz <= pei_Z; ++iz)
                    {
                        Jy(ix, iy, iz)
                            = this->layout_->deriv(Bx, {ix, iy, iz}, DirectionTag<Direction::Z>{})
                              - this->layout_->deriv(Bz, {ix, iy, iz},
                                                     DirectionTag<Direction::X>{});
                    }
                }
            }

            psi_X = this->layout_->physicalStartIndex(Jz, Direction::X);
            pei_X = this->layout_->physicalEndIndex(Jz, Direction::X);
            psi_Y = this->layout_->physicalStartIndex(Jz, Direction::Y);
            pei_Y = this->layout_->physicalEndIndex(Jz, Direction::Y);
            psi_Z = this->layout_->physicalStartIndex(Jz, Direction::Z);
            pei_Z = this->layout_->physicalEndIndex(Jz, Direction::Z);

            for (auto ix = psi_X; ix <= pei_X; ++ix)
            {
                for (auto iy = psi_Y; iy <= pei_Y; ++iy)
                {
                    for (auto iz = psi_Z; iz <= pei_Z; ++iz)
                    {
                        Jz(ix, iy, iz)
                            = this->layout_->deriv(By, {ix, iy, iz}, DirectionTag<Direction::X>{})
                              - this->layout_->deriv(Bx, {ix, iy, iz},
                                                     DirectionTag<Direction::Y>{});
                    }
                }
            }
        }




    public:
        template<typename VecField>
        void operator()(VecField const& B, VecField& J)
        {
            if (!this->hasLayout())
            {
                throw std::runtime_error(
                    "Error - Ampere - GridLayout not set, cannot proceed to calculate ampere()");
            }
            std::cout << "I'm solving ampere " << GridLayout::dimension << "\n";

            this->operator()<VecField, GridLayout::dimension>(B, J);
        }
    };
} // namespace core
} // namespace PHARE



#endif
