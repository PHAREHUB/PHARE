#ifndef PHARE_FARADAY_HPP
#define PHARE_FARADAY_HPP

#include <cstddef>

#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/data/vecfield/vecfield_component.hpp"
#include "core/numerics/differential_operators.hpp"


namespace PHARE::core
{



template<typename GridLayout>
class Faraday
{
    constexpr static auto dimension = GridLayout::dimension;

public:
    Faraday(GridLayout const& layout)
        : layout_{layout}
    {
    }


    template<typename VecField>
    void operator()(VecField const& B, VecField const& E, VecField& Bnew, double dt)
    {
        if (!(B.isUsable() && E.isUsable() && Bnew.isUsable()))
            throw std::runtime_error("Error - Faraday - not all VecField parameters are usable");

        this->dt_ = dt;

        auto const& Bx = B(Component::X);
        auto const& By = B(Component::Y);
        auto const& Bz = B(Component::Z);

        auto& Bxnew = Bnew(Component::X);
        auto& Bynew = Bnew(Component::Y);
        auto& Bznew = Bnew(Component::Z);

        auto const curlX = curl<Component::X>(layout_, E, Bxnew);
        auto const curlY = curl<Component::Y>(layout_, E, Bynew);
        auto const curlZ = curl<Component::Z>(layout_, E, Bznew);

        layout_.evalOnBoxIndependent(
            Bxnew, [&](auto const&... ijk) { Bxnew(ijk...) = Bx(ijk...) - dt_ * curlX({ijk...}); });
        layout_.evalOnBoxIndependent(
            Bynew, [&](auto const&... ijk) { Bynew(ijk...) = By(ijk...) - dt_ * curlY({ijk...}); });
        layout_.evalOnBoxIndependent(
            Bznew, [&](auto const&... ijk) { Bznew(ijk...) = Bz(ijk...) - dt_ * curlZ({ijk...}); });
    }


private:
    double dt_;
    GridLayout layout_;
};

} // namespace PHARE::core


#endif
