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
    class Faraday : public LayoutHolder<GridLayout>
    {
    private:
        template<typename VecField, std::size_t dim, std::enable_if_t<dim == 1, int> = 0>
        void operator()(VecField const& B, VecField const& E, VecField& Bnew, double dt)
        {
            // dBxdt =  0
            // dBydt =  dxEz
            // dBzdt = -dxEy

            auto const& By = B.getComponent(Component::Y);
            auto const& Bz = B.getComponent(Component::Z);

            auto const& Ey = E.getComponent(Component::Y);
            auto const& Ez = E.getComponent(Component::Z);

            auto& Bynew = Bnew.getComponent(Component::Y);
            auto& Bznew = Bnew.getComponent(Component::Z);

            // Direction should not be in gridlayoutdef but in utilities somehow
            auto start = this->layout_->physicalStartIndex(Bynew, Direction::X);
            auto end   = this->layout_->physicalEndIndex(Bynew, Direction::X);

            for (auto ix = start; ix <= end; ++ix)
            {
                Bynew(ix)
                    = By(ix) + dt * this->layout_->deriv(Ez, {ix}, DirectionTag<Direction::X>{});
            }

            start = this->layout_->physicalStartIndex(Bznew, Direction::X);
            end   = this->layout_->physicalEndIndex(Bznew, Direction::X);

            for (auto ix = start; ix <= end; ++ix)
            {
                Bznew(ix)
                    = Bz(ix) - dt * this->layout_->deriv(Ey, {ix}, DirectionTag<Direction::X>{});
            }
        }



        template<typename VecField, std::size_t dim, std::enable_if_t<dim == 2, int> = 0>
        void operator()(VecField const& B, VecField const& E, VecField& Bnew, double dt)
        {
            // dBxdt =  -dyEz
            // dBydt =  dxEz
            // dBzdt = -dxEy + dyEx

            auto const& Bx = B.getComponent(Component::X);
            auto const& By = B.getComponent(Component::Y);
            auto const& Bz = B.getComponent(Component::Z);

            auto const& Ex = E.getComponent(Component::X);
            auto const& Ey = E.getComponent(Component::Y);
            auto const& Ez = E.getComponent(Component::Z);

            auto& Bxnew = Bnew.getComponent(Component::X);
            auto& Bynew = Bnew.getComponent(Component::Y);
            auto& Bznew = Bnew.getComponent(Component::Z);

            auto psi_X = this->layout_->physicalStartIndex(Bxnew, Direction::X);
            auto pei_X = this->layout_->physicalEndIndex(Bxnew, Direction::X);
            auto psi_Y = this->layout_->physicalStartIndex(Bxnew, Direction::Y);
            auto pei_Y = this->layout_->physicalEndIndex(Bxnew, Direction::Y);

            for (auto ix = psi_X; ix <= pei_X; ++ix)
            {
                for (auto iy = psi_Y; iy <= pei_Y; ++iy)
                {
                    Bxnew(ix, iy)
                        = Bx(ix, iy)
                          - dt * this->layout_->deriv(Ez, {ix, iy}, DirectionTag<Direction::Y>{});
                }
            }

            psi_X = this->layout_->physicalStartIndex(Bynew, Direction::X);
            pei_X = this->layout_->physicalEndIndex(Bynew, Direction::X);
            psi_Y = this->layout_->physicalStartIndex(Bynew, Direction::Y);
            pei_Y = this->layout_->physicalEndIndex(Bynew, Direction::Y);

            for (auto ix = psi_X; ix <= pei_X; ++ix)
            {
                for (auto iy = psi_Y; iy <= pei_Y; ++iy)
                {
                    Bynew(ix, iy)
                        = By(ix, iy)
                          + dt * this->layout_->deriv(Ez, {ix, iy}, DirectionTag<Direction::X>{});
                }
            }

            psi_X = this->layout_->physicalStartIndex(Bznew, Direction::X);
            pei_X = this->layout_->physicalEndIndex(Bznew, Direction::X);
            psi_Y = this->layout_->physicalStartIndex(Bznew, Direction::Y);
            pei_Y = this->layout_->physicalEndIndex(Bznew, Direction::Y);

            for (auto ix = psi_X; ix <= pei_X; ++ix)
            {
                for (auto iy = psi_Y; iy <= pei_Y; ++iy)
                {
                    Bznew(ix, iy)
                        = Bz(ix, iy)
                          - dt * this->layout_->deriv(Ey, {ix, iy}, DirectionTag<Direction::X>{})
                          + dt * this->layout_->deriv(Ex, {ix, iy}, DirectionTag<Direction::Y>{});
                }
            }
        }



        template<typename VecField, std::size_t dim, std::enable_if_t<dim == 3, int> = 0>
        void operator()(VecField const& B, VecField const& E, VecField& Bnew, double dt)
        {
            // dBxdt = -dyEz + dzEy
            // dBydt = -dzEx + dxEz
            // dBzdt = -dxEy + dyEx

            auto const& Bx = B.getComponent(Component::X);
            auto const& By = B.getComponent(Component::Y);
            auto const& Bz = B.getComponent(Component::Z);

            auto const& Ex = E.getComponent(Component::X);
            auto const& Ey = E.getComponent(Component::Y);
            auto const& Ez = E.getComponent(Component::Z);

            auto& Bxnew = Bnew.getComponent(Component::X);
            auto& Bynew = Bnew.getComponent(Component::Y);
            auto& Bznew = Bnew.getComponent(Component::Z);

            auto psi_X = this->layout_->physicalStartIndex(Bxnew, Direction::X);
            auto pei_X = this->layout_->physicalEndIndex(Bxnew, Direction::X);
            auto psi_Y = this->layout_->physicalStartIndex(Bxnew, Direction::Y);
            auto pei_Y = this->layout_->physicalEndIndex(Bxnew, Direction::Y);
            auto psi_Z = this->layout_->physicalStartIndex(Bxnew, Direction::Z);
            auto pei_Z = this->layout_->physicalEndIndex(Bxnew, Direction::Z);

            for (auto ix = psi_X; ix <= pei_X; ++ix)
            {
                for (auto iy = psi_Y; iy <= pei_Y; ++iy)
                {
                    for (auto iz = psi_Z; iz <= pei_Z; ++iz)
                    {
                        Bxnew(ix, iy, iz)
                            = Bx(ix, iy, iz)
                              - dt
                                    * this->layout_->deriv(Ez, {ix, iy, iz},
                                                           DirectionTag<Direction::Y>{})
                              + dt
                                    * this->layout_->deriv(Ey, {ix, iy, iz},
                                                           DirectionTag<Direction::Z>{});
                    }
                }
            }

            psi_X = this->layout_->physicalStartIndex(Bynew, Direction::X);
            pei_X = this->layout_->physicalEndIndex(Bynew, Direction::X);
            psi_Y = this->layout_->physicalStartIndex(Bynew, Direction::Y);
            pei_Y = this->layout_->physicalEndIndex(Bynew, Direction::Y);
            psi_Z = this->layout_->physicalStartIndex(Bynew, Direction::Z);
            pei_Z = this->layout_->physicalEndIndex(Bynew, Direction::Z);

            for (auto ix = psi_X; ix <= pei_X; ++ix)
            {
                for (auto iy = psi_Y; iy <= pei_Y; ++iy)
                {
                    for (auto iz = psi_Z; iz <= pei_Z; ++iz)
                    {
                        Bynew(ix, iy, iz)
                            = By(ix, iy, iz)
                              - dt
                                    * this->layout_->deriv(Ex, {ix, iy, iz},
                                                           DirectionTag<Direction::Z>{})
                              + dt
                                    * this->layout_->deriv(Ez, {ix, iy, iz},
                                                           DirectionTag<Direction::X>{});
                    }
                }
            }

            psi_X = this->layout_->physicalStartIndex(Bznew, Direction::X);
            pei_X = this->layout_->physicalEndIndex(Bznew, Direction::X);
            psi_Y = this->layout_->physicalStartIndex(Bznew, Direction::Y);
            pei_Y = this->layout_->physicalEndIndex(Bznew, Direction::Y);
            psi_Z = this->layout_->physicalStartIndex(Bznew, Direction::Z);
            pei_Z = this->layout_->physicalEndIndex(Bznew, Direction::Z);

            for (auto ix = psi_X; ix <= pei_X; ++ix)
            {
                for (auto iy = psi_Y; iy <= pei_Y; ++iy)
                {
                    for (auto iz = psi_Z; iz <= pei_Z; ++iz)
                    {
                        Bznew(ix, iy, iz)
                            = Bz(ix, iy, iz)
                              - dt
                                    * this->layout_->deriv(Ey, {ix, iy, iz},
                                                           DirectionTag<Direction::X>{})
                              + dt
                                    * this->layout_->deriv(Ex, {ix, iy, iz},
                                                           DirectionTag<Direction::Y>{});
                    }
                }
            }
        }


    public:
        template<typename VecField>
        void operator()(VecField const& B, VecField const& E, VecField& Bnew, double dt)
        {
            if (!this->hasLayout())
            {
                throw std::runtime_error(
                    "Error - Faraday - GridLayout not set, cannot proceed to calculate faraday()");
            }

            this->operator()<VecField, GridLayout::dimension>(B, E, Bnew, dt);
        }
    };
} // namespace core
} // namespace PHARE



#endif
