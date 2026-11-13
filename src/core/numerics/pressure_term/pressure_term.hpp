#ifndef PHARE_CORE_NUMERICS_PRESSURE_TERM_HPP
#define PHARE_CORE_NUMERICS_PRESSURE_TERM_HPP

#include "core/data/field/field.hpp"
#include "core/data/vecfield/vecfield_component.hpp"
#include "core/numerics/differential_operators.hpp"

#include <cmath>
#include <stdexcept>
#include <vector>


namespace PHARE::core
{
// The electron pressure term of Ohm's law, -grad(Pe)/n, on the E edges.
//
// The electron pressure closure computes it, not Ohm, because the discrete property it must keep
// depends on the pressure law. When Pe is a function of n only, -grad(Pe)/n is the gradient of a
// node quantity (the enthalpy). Written as a discrete gradient, its discrete curl vanishes exactly;
// any other form lets density noise create magnetic field (#1527). Ohm only adds the result to E,
// so it is computed on the box Ohm evaluates E on: the domain box of each E component.


// Isothermal closure, Pe = Te n with a uniform Te: the pressure term is -Te grad(ln n), the
// discrete gradient of the node quantity Te ln(n).
// n is primal, so grad(ln n) lands on the E edges; it must be > 0 on the domain nodes.
template<typename GridLayout>
class IsothermalPressureTerm
{
public:
    explicit IsothermalPressureTerm(double const Te)
        : Te_{Te}
    {
    }

    template<typename Density, typename VecField>
    void operator()(GridLayout const& layout, Density const& n, VecField& pressureTerm)
    {
        using NodeField = Field<GridLayout::dimension, typename Density::physical_quantity_type>;

        // E's domain box reads n only on its domain nodes: ln(n) is computed there only, and its
        // ghost nodes are left as the previous patch wrote them
        lnNData_.resize(n.size());
        auto lnN = NodeField{"lnN", n.physicalQuantity(), lnNData_.data(), n.shape()};
        layout.evalOnBox(n, [&](auto const&... ijk) { lnN(ijk...) = std::log(n(ijk...)); });

        auto const gradLnN = grad(layout, lnN);
        auto const fill    = [&](auto& component, auto const& partial) {
            if (layout.centering(component) != partial.centering())
                throw std::runtime_error(
                    "isothermal pressure term: grad(ln n) does not land on "
                       "the pressure term (n must be primal, the term E-centred)");
            layout.evalOnBox(component, [&](auto const&... ijk) {
                component(ijk...) = -Te_ * partial({ijk...});
            });
        };
        fill(pressureTerm(Component::X), gradLnN.x);
        fill(pressureTerm(Component::Y), gradLnN.y);
        fill(pressureTerm(Component::Z), gradLnN.z);
    }

private:
    double Te_;
    std::vector<double> lnNData_; // reused from one patch to the next
};

} // namespace PHARE::core

#endif
