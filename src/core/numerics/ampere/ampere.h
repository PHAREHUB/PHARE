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
        // Jx = dyBz - dzBy = 0
        // Jy = dzBx - dxBz = -dxBz
        // Jz = dxBy - dyBx = dxBy
        auto &Jx = J.getComponent(Component::X);
        auto &Jy = J.getComponent(Component::Y);
        auto &Jz = J.getComponent(Component::Z);

        auto const &Bx = B.getComponent(Component::X);
        auto const &By = B.getComponent(Component::Y);
        auto const &Bz = B.getComponent(Component::Z);

        // TODO Direction should not be in gridlayoutdef but in utilities somehow
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
        this->layout_.deriv();
        std::cout << "solving ampere2D\n";
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
        this->layout_.deriv();
        std::cout << "solving ampere3D\n";
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
