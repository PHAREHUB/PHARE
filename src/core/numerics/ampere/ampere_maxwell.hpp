#ifndef PHARE_CORE_NUMERICS_AMPERE_AMPERE_MAXWELL_HPP
#define PHARE_CORE_NUMERICS_AMPERE_AMPERE_MAXWELL_HPP

#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/data/vecfield/vecfield_component.hpp"
#include "core/numerics/differential_operators.hpp"


namespace PHARE::core
{
/**
 * @brief AmpereMaxwell advances the electric field with the full Ampere equation, in PHARE
 * normalized units: dE/dt = (c/V_A)^2 (curl B - J).
 * Only the physical nodes of E are updated, ghost nodes are left untouched and must be
 * filled by communications afterwards.
 */
template<typename GridLayout>
class AmpereMaxwell
{
public:
    AmpereMaxwell(GridLayout const& layout, double const c2)
        : layout_{layout}
        , c2_{c2}
    {
    }


    template<typename VecField>
    void operator()(VecField const& B, VecField const& J, VecField& E, double const dt)
    {
        // can't use structured bindings because
        //   "reference to local binding declared in enclosing function"
        auto const& Jx = J(Component::X);
        auto const& Jy = J(Component::Y);
        auto const& Jz = J(Component::Z);

        auto& Ex = E(Component::X);
        auto& Ey = E(Component::Y);
        auto& Ez = E(Component::Z);

        auto const curlB = curl(layout_, B);
        auto const coef  = dt * c2_;

        layout_.evalOnBox(
            Ex, [&](auto const&... ijk) { Ex(ijk...) += coef * (curlB.x({ijk...}) - Jx(ijk...)); });
        layout_.evalOnBox(
            Ey, [&](auto const&... ijk) { Ey(ijk...) += coef * (curlB.y({ijk...}) - Jy(ijk...)); });
        layout_.evalOnBox(
            Ez, [&](auto const&... ijk) { Ez(ijk...) += coef * (curlB.z({ijk...}) - Jz(ijk...)); });
    }


private:
    GridLayout layout_;
    double c2_;
};

} // namespace PHARE::core
#endif
