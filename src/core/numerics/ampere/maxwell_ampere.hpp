#ifndef PHARE_CORE_NUMERICS_MAXWELL_AMPERE_HPP
#define PHARE_CORE_NUMERICS_MAXWELL_AMPERE_HPP

#include <cstddef>

#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/data/grid/gridlayout_utils.hpp"
#include "core/data/vecfield/vecfield_component.hpp"


namespace PHARE::core
{
template<typename GridLayout>
class MaxwellAmpere : public LayoutHolder<GridLayout>
{
    constexpr static auto dimension = GridLayout::dimension;
    using LayoutHolder<GridLayout>::layout_;

public:
    template<typename VecField>
    void operator()(VecField const& B, VecField const& E, VecField const& J, VecField& Enew, double dt)
    {
        if (!this->hasLayout())
            throw std::runtime_error(
                "Error - MaxwellAmpere - GridLayout not set, cannot proceed to calculate MaxwellAmpere()");

        if (!(B.isUsable() && E.isUsable() && J.isUsable() && Bnew.isUsable()))
            throw std::runtime_error("Error - MaxwellAmpere - not all VecField parameters are usable");

        this->dt_ = dt;

        // can't use structured bindings because
        //   "reference to local binding declared in enclosing function"
        auto const& Ex = E(Component::X);
        auto const& Ey = E(Component::Y);
        auto const& Ez = E(Component::Z);

        auto const& Jx = J(Component::X);
        auto const& Jy = J(Component::Y);
        auto const& Jz = J(Component::Z);

        auto& Exnew = Enew(Component::X);
        auto& Eynew = Enew(Component::Y);
        auto& Eznew = Enew(Component::Z);

        layout_->evalOnBox(Exnew,
                                [&](auto&... args) mutable { ExEq_(Ex, B, Jx, Exnew, args...); });
        layout_->evalOnBox(Eynew,
                                [&](auto&... args) mutable { EyEq_(Ey, B, Jx, Eynew, args...); });
        layout_->evalOnBox(Eznew,
                                [&](auto&... args) mutable { EzEq_(Ez, B, Jx, Eznew, args...); });
    }


private:
    double dt_;
    double c_ = 299792458.0; // m/s
    double m_p = 1.6726219e-27; // kg
    double mu_0 = 1.25663706e-6; // N/A^2
    double B_0 = 1.; // placeholder reference magnetic field, SI units
    double n_0 = 100.; // placeholder reference density, SI units
    double Va = B_0 / std::sqrt(mu_0 * m_p * n_0); // normalized Alfven speed, non-relativistic (CHECK)
    double c_norm = c_ / Va; // normalized velocity
    double c2 = c_norm * c_norm;
    double inv_c2 = 1.0 / c2;

    template<typename VecField, typename Field, typename... Indexes>
    void ExEq_(Field const& Ex, VecField const& E, VecField const& J, Field& Exnew, Indexes const&... ijk) const
    {
        auto const& [_, By, Bz] = B();

        if constexpr (dimension == 1)
            Exnew(ijk...) = Ex(ijk...) - dt_ * inv_c2 * Jx(ijk...) ; 

        if constexpr (dimension == 2)
            Exnew(ijk...) = Ex(ijk...) + dt_ * inv_c2 * (layout_->template deriv<Direction::Y>(Bz, {ijk...}) 
                            - Jx(ijk...) );

        if constexpr (dimension == 3)
            Exnew(ijk...) = Ex(ijk...) + dt_ * inv_c2 * (layout_->template deriv<Direction::Y>(Bz, {ijk...})
                            - layout_->template deriv<Direction::Z>(By, {ijk...}) - Jx(ijk...));
    }

    template<typename VecField, typename Field, typename... Indexes>
    void EyEq_(Field const& Ey, VecField const& E, VecField const& J, Field& Eynew, Indexes const&... ijk) const
    {
        auto const& [Bx, _, Bz] = B();

        if constexpr (dimension == 1 || dimension == 2)
            Eynew(ijk...) = Ey(ijk...) - dt_ * inv_c2 *( layout_->template deriv<Direction::X>(Bz, {ijk...}) 
                            + Jy(ijk...));

        if constexpr (dimension == 3)
            Eynew(ijk...) = Ey(ijk...) + dt_ * inv_c2 * (layout_->template deriv<Direction::Z>(Bx, {ijk...})
                            - layout_->template deriv<Direction::X>(Bz, {ijk...}) -  Jy(ijk...));
    }

    template<typename VecField, typename Field, typename... Indexes>
    void EzEq_(Field const& Ez, VecField const& E, VecField const& J, Field& Eznew, Indexes const&... ijk) const
    {
        auto const& [Bx, By, _] = B();

        if constexpr (dimension == 1)
            Eznew(ijk...) = Ez(ijk...) + dt_ * inv_c2 * (layout_->template deriv<Direction::X>(By, {ijk...}) 
                            - Jz(ijk...));

        if constexpr (dimension == 2 || dimension == 3)
            Eznew(ijk...) = Ez(ijk...) - dt_ * inv_c2 * (layout_->template deriv<Direction::X>(By, {ijk...})
                            - layout_->template deriv<Direction::Y>(Bx, {ijk...})- Jz(ijk...));
    }
};

} // namespace PHARE::core

#endif
