#ifndef PHARE_FARADAY_HPP
#define PHARE_FARADAY_HPP

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
    void operator()(VecField const& B, VecField const& E, VecField& Bnew, double dt) _PHARE_ALL_FN_
    {
        dt_ = dt;
        // can't use structured bindings because
        //   "reference to local binding declared in enclosing function"
        auto const& Bx = B(Component::X);
        auto const& By = B(Component::Y);
        auto const& Bz = B(Component::Z);

        auto& Bxnew = Bnew(Component::X);
        auto& Bynew = Bnew(Component::Y);
        auto& Bznew = Bnew(Component::Z);

        auto const curlE      = curl(layout_, E);
        basic::Field const bx = Bx, by = By, bz = Bz;
        basic::Field bxnew = Bxnew, bynew = Bynew, bznew = Bznew;

        layout_.evalOnBox(Bxnew, [=] _PHARE_ALL_FN_(auto const&... ijk) mutable {
            bxnew(ijk...) = bx(ijk...) - dt * curlE.x({ijk...});
        });
        layout_.evalOnBox(Bynew, [=] _PHARE_ALL_FN_(auto const&... ijk) mutable {
            bynew(ijk...) = by(ijk...) - dt * curlE.y({ijk...});
        });
        layout_.evalOnBox(Bznew, [=] _PHARE_ALL_FN_(auto const&... ijk) mutable {
            bznew(ijk...) = bz(ijk...) - dt * curlE.z({ijk...});
        });
    }

private:
    GridLayout layout_;
    double dt_;
};

} // namespace PHARE::core

#endif
