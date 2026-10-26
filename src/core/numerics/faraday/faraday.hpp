#ifndef PHARE_FARADAY_HPP
#define PHARE_FARADAY_HPP

#include <cstddef>

#include "core/def.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/data/vecfield/vecfield_component.hpp"
#include "core/numerics/differential_operators.hpp"


namespace PHARE::core
{



template<typename GridLayout>
class Faraday
{
public:
    Faraday(GridLayout const& layout)
        : layout_{layout}
    {
    }


    template<typename VecField>
    void operator()(VecField const& B, VecField const& E, VecField& Bnew, double dt)
    {
        if (!isUsable(B, E, Bnew))
            throw std::runtime_error("Error - Faraday - not all VecField parameters are usable");

        this->dt_ = dt;

        // can't use structured bindings because
        //   "reference to local binding declared in enclosing function"
        auto const& Bx = B(Component::X);
        auto const& By = B(Component::Y);
        auto const& Bz = B(Component::Z);

        auto& Bxnew = Bnew(Component::X);
        auto& Bynew = Bnew(Component::Y);
        auto& Bznew = Bnew(Component::Z);

        auto const curlE = curl(layout_, E);

        layout_.evalOnBox(Bxnew, [&](auto const&... ijk) {
            Bxnew(ijk...) = Bx(ijk...) - dt_ * curlE.x({ijk...});
        });
        layout_.evalOnBox(Bynew, [&](auto const&... ijk) {
            Bynew(ijk...) = By(ijk...) - dt_ * curlE.y({ijk...});
        });
        layout_.evalOnBox(Bznew, [&](auto const&... ijk) {
            Bznew(ijk...) = Bz(ijk...) - dt_ * curlE.z({ijk...});
        });
    }


private:
    double dt_;
    GridLayout layout_;
};

} // namespace PHARE::core


#endif
