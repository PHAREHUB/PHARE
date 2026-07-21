#ifndef PHARE_CORE_BOUNDARY_BOUNDARY_INFLOW_COMPOSE_HPP
#define PHARE_CORE_BOUNDARY_BOUNDARY_INFLOW_COMPOSE_HPP

#include "initializer/data_provider.hpp"

#include <array>
#include <cstddef>
#include <memory>
#include <utility>
#include <vector>

namespace PHARE::core::inflow_compose
{
// All helpers preserve the batch-evaluation invariant: each input function is invoked
// exactly once over the whole coordinate batch, then combined element-wise. They mirror
// the shape of the earlier BoundaryFactory constant-B motional-E composition.

/** @brief A constant lifted into a space-time function: returns c at every node. Sized to
 * the first spatial coordinate span (the node count). Time is ignored. */
template<std::size_t dim>
initializer::SpaceTimeFunction<dim> constFunction(double const c)
{
    return [c](auto const&... args) -> std::shared_ptr<Span<double>> {
        auto const& first = std::get<0>(std::forward_as_tuple(args...));
        std::vector<double> out(first.size(), c);
        return std::make_shared<VectorSpan<double>>(std::move(out));
    };
}

/** @brief Element-wise product f*g of two space-time functions. */
template<std::size_t dim>
initializer::SpaceTimeFunction<dim> mulFunction(initializer::SpaceTimeFunction<dim> f,
                                                initializer::SpaceTimeFunction<dim> g)
{
    return [f = std::move(f), g = std::move(g)](
               auto const&... args) -> std::shared_ptr<Span<double>> {
        auto sf = f(args...);
        auto sg = g(args...);
        std::vector<double> out(sf->size());
        for (std::size_t k = 0; k < out.size(); ++k)
            out[k] = (*sf)[k] * (*sg)[k];
        return std::make_shared<VectorSpan<double>>(std::move(out));
    };
}

/** @brief a*f1*g1 + b*f2*g2, element-wise. Building block for -v x B components. */
template<std::size_t dim>
initializer::SpaceTimeFunction<dim>
prodComb2(double const a, initializer::SpaceTimeFunction<dim> f1,
          initializer::SpaceTimeFunction<dim> g1, double const b,
          initializer::SpaceTimeFunction<dim> f2, initializer::SpaceTimeFunction<dim> g2)
{
    return [a, b, f1 = std::move(f1), g1 = std::move(g1), f2 = std::move(f2),
            g2 = std::move(g2)](auto const&... args) -> std::shared_ptr<Span<double>> {
        auto s1 = f1(args...);
        auto h1 = g1(args...);
        auto s2 = f2(args...);
        auto h2 = g2(args...);
        std::vector<double> out(s1->size());
        for (std::size_t k = 0; k < out.size(); ++k)
            out[k] = a * (*s1)[k] * (*h1)[k] + b * (*s2)[k] * (*h2)[k];
        return std::make_shared<VectorSpan<double>>(std::move(out));
    };
}

/** @brief The three components of the ideal motional field E = -V x B from two
 * function-vectors. Generalizes the constant-V motional field to a time-varying V.
 *   E_x = -(V_y B_z - V_z B_y) = -V_y B_z + V_z B_y
 *   E_y = -(V_z B_x - V_x B_z) = -V_z B_x + V_x B_z
 *   E_z = -(V_x B_y - V_y B_x) = -V_x B_y + V_y B_x */
template<std::size_t dim>
std::array<initializer::SpaceTimeFunction<dim>, 3>
negCrossFunction(std::array<initializer::SpaceTimeFunction<dim>, 3> const& V,
                 std::array<initializer::SpaceTimeFunction<dim>, 3> const& B)
{
    return {prodComb2<dim>(-1.0, V[1], B[2], 1.0, V[2], B[1]),
            prodComb2<dim>(-1.0, V[2], B[0], 1.0, V[0], B[2]),
            prodComb2<dim>(-1.0, V[0], B[1], 1.0, V[1], B[0])};
}

} // namespace PHARE::core::inflow_compose

#endif // PHARE_CORE_BOUNDARY_BOUNDARY_INFLOW_COMPOSE_HPP
