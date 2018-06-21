#ifndef PHARE_CORE_NUMERICS_AMPERE_AMPERE_H
#define PHARE_CORE_NUMERICS_AMPERE_AMPERE_H

#include <cstddef>
#include <iostream>

#include "data/grid/gridlayoutdefs.h"
#include "data/vecfield/vecfield_component.h"
#include "utilities/index/index.h"

namespace PHARE
{
/** @brief Actual implementation of the Ampere equation
 *
 * This implementation is used by the Ampere object to calculate the current density J.
 * It is templated by the layout and the dimension so specialized template actually do the job
 *
 */
template<typename GridLayout, std::size_t dim>
class AmpereImpl
{
};


/** @brief Base class for AmpereImpl specialized classes
 * to factorize the code that deals with the GridLayout.
 *
 * The purpose of this class is to store the GridLayout
 * pointer of Ampere, and propose the interface to check
 * whether it's been set or not.
 *
 */
template<typename GridLayout>
class AmpereImplInternals
{
protected:
    GridLayout *layout_{nullptr};

public:
    /**
     * @brief hasLayoutSet returns true if the Layout has been set
     */
    bool hasLayoutSet() const { return (layout_ == nullptr) ? false : true; }


    /**
     * @brief setLayout is used to give Ampere a pointer to a gridlayout
     */
    void setLayout(GridLayout *layout)
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
class AmpereImpl<GridLayout, 1> : public AmpereImplInternals<GridLayout>
{
    static_assert(GridLayout::dimension == 1, "Error: Passed non-1D GridLayout to 1D AmpereImpl");

public:
    template<typename VecField>
    void operator()(VecField const &B, VecField &J)
    {
        //auto &Jx = J.getComponent(Component::X); // =  0
        auto &Jy = J.getComponent(Component::Y); // = -dxBz
        auto &Jz = J.getComponent(Component::Z); // =  dxBy

        //auto const &Bx = B.getComponent(Component::X);
        auto const &By = B.getComponent(Component::Y);
        auto const &Bz = B.getComponent(Component::Z);

        // TODO Direction should not be in gridlayoutdef but in utilities somehow
        // TODO 1st arg of physicalStartIndex( could be QtyCentering::primal ?
        auto start = this->layout_->physicalStartIndex(Jy, Direction::X);
        auto end   = this->layout_->physicalEndIndex(Jy, Direction::X);

        for (auto ix = start; ix <= end; ++ix)
        {
            Jy(ix) = -this->layout_->deriv(Bz, make_index(ix), DirectionTag<Direction::X>{});
        }

        start = this->layout_->physicalStartIndex(Jz, Direction::X);
        end   = this->layout_->physicalEndIndex(Jz, Direction::X);

        for (auto ix = start; ix <= end; ++ix)
        {
            Jz(ix) = this->layout_->deriv(By, make_index(ix), DirectionTag<Direction::X>{});
        }
    }
};



template<typename GridLayout>
class AmpereImpl<GridLayout, 2> : public AmpereImplInternals<GridLayout>
{
    static_assert(GridLayout::dimension == 2, "Error: Passed non-2D GridLayout to 2D AmpereImpl");

public:
    template<typename VecField>
    void operator()(VecField const &B, VecField &J)
    {
        auto &Jx = J.getComponent(Component::X); // =  dyBz
        auto &Jy = J.getComponent(Component::Y); // = -dxBz
        auto &Jz = J.getComponent(Component::Z); // =  dxBy - dyBx

        auto const &Bx = B.getComponent(Component::X);
        auto const &By = B.getComponent(Component::Y);
        auto const &Bz = B.getComponent(Component::Z);

        auto psi_X = this->layout_->physicalStartIndex(Jx, Direction::X);
        auto pei_X = this->layout_->physicalEndIndex(Jx, Direction::X);
        auto psi_Y = this->layout_->physicalStartIndex(Jx, Direction::Y);
        auto pei_Y = this->layout_->physicalEndIndex(Jx, Direction::Y);

        for (auto ix = psi_X; ix <= pei_X; ++ix)
        {
            for (auto iy = psi_Y; iy <= pei_Y; ++iy)
            {
                Jx(ix, iy) = this->layout_->deriv(Bz, make_index(ix, iy), DirectionTag<Direction::Y>{});
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
                Jy(ix, iy) = -this->layout_->deriv(Bz, make_index(ix, iy), DirectionTag<Direction::X>{});
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
                Jz(ix, iy) = this->layout_->deriv(By, make_index(ix, iy), DirectionTag<Direction::X>{})
                            -this->layout_->deriv(Bx, make_index(ix, iy), DirectionTag<Direction::Y>{});
            }
        }
    }
};



template<typename GridLayout>
class AmpereImpl<GridLayout, 3> : public AmpereImplInternals<GridLayout>
{
    static_assert(GridLayout::dimension == 3, "Error: Passed non-3D GridLayout to 3D AmpereImpl");

public:
    template<typename VecField>
    void operator()(VecField const &B, VecField &J)
    {
        auto &Jx = J.getComponent(Component::X); // =  dyBz - dzBx
        auto &Jy = J.getComponent(Component::Y); // =  dzBx - dxBz
        auto &Jz = J.getComponent(Component::Z); // =  dxBy - dyBx

        auto const &Bx = B.getComponent(Component::X);
        auto const &By = B.getComponent(Component::Y);
        auto const &Bz = B.getComponent(Component::Z);

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
                    Jx(ix, iy, iz) = this->layout_->deriv(Bz, make_index(ix, iy, iz), DirectionTag<Direction::Y>{})
                                    -this->layout_->deriv(By, make_index(ix, iy, iz), DirectionTag<Direction::Z>{});
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
                    Jy(ix, iy, iz) = this->layout_->deriv(Bx, make_index(ix, iy, iz), DirectionTag<Direction::Z>{})
                                    -this->layout_->deriv(Bz, make_index(ix, iy, iz), DirectionTag<Direction::X>{});
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
                    Jz(ix, iy, iz) = this->layout_->deriv(By, make_index(ix, iy, iz), DirectionTag<Direction::X>{})
                                    -this->layout_->deriv(Bx, make_index(ix, iy, iz), DirectionTag<Direction::Y>{});
                }
            }
        }
    }
};




template<typename GridLayout>
class Ampere
{
private:
    AmpereImpl<GridLayout, GridLayout::dimension> impl_;

public:
    template<typename VecField>
    void operator()(VecField const &B, VecField &J)
    {
        if (!impl_.hasLayoutSet())
        {
            throw std::runtime_error(
                "Error - Ampere - GridLayout not set, cannot proceed to calculate ampere()");
        }
        std::cout << "I'm solving ampere " << GridLayout::dimension << "\n";

        impl_(B, J);
    }


    void setLayout(GridLayout *layout) { impl_.setLayout(layout); }
};
} // namespace PHARE



#endif
