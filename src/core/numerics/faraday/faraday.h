#ifndef PHARE_CORE_NUMERICS_AMPERE_AMPERE_H
#define PHARE_CORE_NUMERICS_AMPERE_AMPERE_H

#include <cstddef>
#include <iostream>

#include "data/grid/gridlayoutdefs.h"
#include "data/vecfield/vecfield_component.h"
#include "utilities/index/index.h"

namespace PHARE
{
namespace core
{
    /** @brief Actual implementation of the Faraday equation
     *
     * This implementation is used by the Faraday object to calculate the magentic field
     * It is templated by the layout and the dimension so specialized template actually do the job
     *
     */
    template<typename GridLayout, std::size_t dim>
    class FaradayImpl
    {
    };


    /** @brief Base class for FaradayImpl specialized classes
     * to factorize the code that deals with the GridLayout.
     *
     * The purpose of this class is to store the GridLayout
     * pointer of Faraday, and propose the interface to check
     * whether it's been set or not.
     *
     */
    template<typename GridLayout>
    class FaradayImplInternals
    {
    protected:
        GridLayout* layout_{nullptr};
        const double dt_ = 1.0; // TODO is it where time step should be ?

    public:
        /**
         * @brief hasLayoutSet returns true if the Layout has been set
         */
        bool hasLayoutSet() const { return (layout_ == nullptr) ? false : true; }


        /**
         * @brief setLayout is used to give Ampere a pointer to a gridlayout
         */
        void setLayout(GridLayout* layout)
        {
            if (layout_ != nullptr)
            {
                throw std::runtime_error(
                    "Error - ampere - cannot set layout_ because it is already set");
            }
            else
            {
                layout_ = layout;
            }
        }
    };


    /** @brief 1D specialization of the ampere equation solver implementation
     */
    template<typename GridLayout>
    class FaradayImpl<GridLayout, 1> : public FaradayImplInternals<GridLayout>
    {
        static_assert(GridLayout::dimension == 1,
                      "Error: Passed non-1D GridLayout to 1D FaradayImpl");

    public:
        template<typename VecField>
        void operator()(VecField const& B, VecField const& E, VecField& Bnew)
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
                    = By(ix)
                      + this->dt_ * this->layout_->deriv(Ez, {ix}, DirectionTag<Direction::X>{});
            }

            start = this->layout_->physicalStartIndex(Bznew, Direction::X);
            end   = this->layout_->physicalEndIndex(Bznew, Direction::X);

            for (auto ix = start; ix <= end; ++ix)
            {
                Bznew(ix)
                    = Bz(ix)
                      - this->dt_ * this->layout_->deriv(Ey, {ix}, DirectionTag<Direction::X>{});
            }
        }
    };



    template<typename GridLayout>
    class FaradayImpl<GridLayout, 2> : public FaradayImplInternals<GridLayout>
    {
        static_assert(GridLayout::dimension == 2,
                      "Error: Passed non-2D GridLayout to 2D FaradayImpl");

    public:
        template<typename VecField>
        void operator()(VecField const& B, VecField const& E, VecField& Bnew)
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
                          - this->dt_
                                * this->layout_->deriv(Ez, {ix, iy}, DirectionTag<Direction::Y>{});
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
                          + this->dt_
                                * this->layout_->deriv(Ez, {ix, iy}, DirectionTag<Direction::X>{});
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
                          - this->dt_
                                * this->layout_->deriv(Ey, {ix, iy}, DirectionTag<Direction::X>{})
                          + this->dt_
                                * this->layout_->deriv(Ex, {ix, iy}, DirectionTag<Direction::Y>{});
                }
            }
        }
    };



    template<typename GridLayout>
    class FaradayImpl<GridLayout, 3> : public FaradayImplInternals<GridLayout>
    {
        static_assert(GridLayout::dimension == 3,
                      "Error: Passed non-3D GridLayout to 3D FaradayImpl");

    public:
        template<typename VecField>
        void operator()(VecField const& B, VecField const& E, VecField& Bnew)
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
                              - this->dt_
                                    * this->layout_->deriv(Ez, {ix, iy, iz},
                                                           DirectionTag<Direction::Y>{})
                              + this->dt_
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
                              - this->dt_
                                    * this->layout_->deriv(Ex, {ix, iy, iz},
                                                           DirectionTag<Direction::Z>{})
                              + this->dt_
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
                              - this->dt_
                                    * this->layout_->deriv(Ey, {ix, iy, iz},
                                                           DirectionTag<Direction::X>{})
                              + this->dt_
                                    * this->layout_->deriv(Ex, {ix, iy, iz},
                                                           DirectionTag<Direction::Y>{});
                    }
                }
            }
        }
    };




    template<typename GridLayout>
    class Faraday
    {
    private:
        FaradayImpl<GridLayout, GridLayout::dimension> impl_;

    public:
        template<typename VecField>
        void operator()(VecField const& B, VecField const& E, VecField& Bnew)
        {
            if (!impl_.hasLayoutSet())
            {
                throw std::runtime_error(
                    "Error - Faraday - GridLayout not set, cannot proceed to calculate faraday()");
            }

            impl_(B, E, Bnew);
        }


        void setLayout(GridLayout* layout) { impl_.setLayout(layout); }
    };
} // namespace core
} // namespace PHARE



#endif
