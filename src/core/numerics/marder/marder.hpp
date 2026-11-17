#ifndef PHARE_CORE_NUMERICS_MARDER_MARDER_HPP
#define PHARE_CORE_NUMERICS_MARDER_MARDER_HPP

#include "core/data/field/field_traits.hpp"
#include "core/data/field/field_scratch.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/utilities/box/box.hpp"
#include "core/utilities/box/box_spans.hpp"

#include <array>
#include <cassert>
#include <cstddef>
#include <vector>


namespace PHARE::core
{
/**
 * @brief Marder applies one step of Marder (Langdon) divergence cleaning to the electric field:
 *
 *     F  = div E - c2 rho                on the nodes of rho (primal, as div E on Yee)
 *     E += d dt grad F                   on the nodes of each E component
 *
 * c2 = (c/V_A)^2 is the coefficient of Gauss's law in PHARE units, as in AmpereMaxwell.
 * F then diffuses: F <- F + d dt Laplacian(F), with the compact discrete Laplacian, so a mode
 * of wave vector k is multiplied by 1 - d dt sum_i 4/dl_i^2 sin^2(k_i dl_i / 2). With
 * s = d dt sum_i 2/dl_i^2, s <= 1 is stable, s = 1/2 removes the grid scale mode in one step.
 * The correction is a discrete gradient: its discrete curl is zero, B is not affected.
 *
 * Only the physical nodes of E are written. F is computed on the physical nodes of rho, which
 * reads one ghost layer of E (the dual node before the first primal one); rho is only read on
 * its physical nodes. The E ghosts must be filled again afterwards. F lives in the field scratch
 * buffer, which grows once and is then reused.
 *
 * F is set to 0 on the skipped nodes, given as boxes of AMR indexes of the primal nodes (node i is
 * the lower corner of cell i): their residual is neither removed nor spread, a Dirichlet condition
 * F = 0. The fully kinetic solver skips the nodes on the edge of its level and of the region a
 * finer level synchronized (fk_electric.hpp, regionEdgeNodes).
 */
template<typename GridLayout>
class Marder
{
    static constexpr auto dimension = GridLayout::dimension;
    using Nodes_t                   = std::array<std::size_t, 3>;

    // physical nodes and strides of a field, indexes padded to 3 directions at the front, the
    // last one being contiguous
    struct Nodes
    {
        Nodes_t lo{0, 0, 0}, hi{0, 0, 0}, strides{0, 0, 1};
    };

public:
    Marder(GridLayout const& layout, double const c2, double const d)
        : layout_{layout}
        , c2_{c2}
        , d_{d}
    {
    }


    template<typename VecField, IsField Field>
    void operator()(VecField& E, Field const& rho, double const dt,
                    std::vector<Box<int, dimension>> const& skipped = {}) const
    {
        constexpr auto pad = 3 - dimension;
        assert(GridLayout::centering(rho.physicalQuantity()) == allPrimal_());

        auto const rhoNodes = nodes_(rho);
        ScratchWrite scratch;
        auto* const F = scratch.get().get_no_copy(rho.size()).data();

        // F = -c2 rho, then + dE_d/dx_d direction by direction: dual node i is at i + 1/2
        forEach_(rhoNodes, [&](auto const n) { F[n] = -c2_ * rho.data()[n]; });
        for (std::size_t d = 0; d < dimension; ++d)
        {
            auto const& Ed     = E[d];
            auto const ENodes  = nodes_(Ed);
            auto const stride  = ENodes.strides[pad + d];
            auto const inverse = layout_.inverseMeshSize(static_cast<Direction>(d));
            forEach_(rhoNodes, ENodes, [&](auto const n, auto const m) {
                F[n] += (Ed.data()[m] - Ed.data()[m - stride]) * inverse;
            });
        }
        // no residual on the skipped nodes
        for (auto const& box : skipped)
            for_each_point(box, [&](auto const& node) { F[index_(rhoNodes, node)] = 0.; });

        // E_d += d dt dF/dx_d, from the primal nodes i and i + 1 around the dual node i. The
        // components along directions the grid does not have are not changed.
        auto const coef = d_ * dt;
        for (std::size_t d = 0; d < dimension; ++d)
        {
            auto& Ed          = E[d];
            auto const stride = rhoNodes.strides[pad + d];
            auto const factor = coef * layout_.inverseMeshSize(static_cast<Direction>(d));
            forEach_(nodes_(Ed), rhoNodes, [&](auto const m, auto const n) {
                Ed.data()[m] += factor * (F[n + stride] - F[n]);
            });
        }
    }


private:
    static constexpr auto allPrimal_()
    {
        std::array<QtyCentering, dimension> centering;
        centering.fill(QtyCentering::primal);
        return centering;
    }

    template<typename Field>
    Nodes nodes_(Field const& field) const
    {
        constexpr auto pad = 3 - dimension;
        Nodes_t shape{1, 1, 1};
        Nodes nodes;
        for (std::size_t d = 0; d < dimension; ++d)
        {
            auto const direction = static_cast<Direction>(d);
            shape[pad + d]       = field.shape()[d];
            nodes.lo[pad + d]    = layout_.physicalStartIndex(field, direction);
            nodes.hi[pad + d]    = layout_.physicalEndIndex(field, direction);
            assert(nodes.lo[pad + d] >= 1 and nodes.hi[pad + d] + 1 < shape[pad + d]);
        }
        nodes.strides = {shape[1] * shape[2], shape[2], 1};
        return nodes;
    }

    // index in the field of the primal node of AMR index `node`, a physical one
    template<typename Index>
    std::size_t index_(Nodes const& nodes, Index const& node) const
    {
        constexpr auto pad = 3 - dimension;
        std::size_t n      = 0;
        for (std::size_t d = 0; d < dimension; ++d)
        {
            auto const offset = node[d] - layout_.AMRBox().lower[d];
            assert(offset >= 0);
            auto const local = nodes.lo[pad + d] + static_cast<std::size_t>(offset);
            assert(local <= nodes.hi[pad + d]);
            n += local * nodes.strides[pad + d];
        }
        return n;
    }

    // fn(n) on the nodes lo ... hi (inclusive) of a field
    template<typename Fn>
    static void forEach_(Nodes const& nodes, Fn&& fn)
    {
        forEach_(nodes, nodes, [&](auto const n, auto) { fn(n); });
    }

    // fn(n, m) with n the index in the first field of the nodes lo ... hi of that field, m the
    // index of the node with the same local indexes in the second field
    template<typename Fn>
    static void forEach_(Nodes const& first, Nodes const& second, Fn&& fn)
    {
        auto const& [lo, hi, s] = first;
        auto const& t           = second.strides;
        for (auto i = lo[0]; i <= hi[0]; ++i)
            for (auto j = lo[1]; j <= hi[1]; ++j)
            {
                auto const row  = i * s[0] + j * s[1];
                auto const row2 = i * t[0] + j * t[1];
                for (auto k = lo[2]; k <= hi[2]; ++k)
                    fn(row + k, row2 + k);
            }
    }


    GridLayout const& layout_;
    double const c2_;
    double const d_;
};


/**
 * @brief marderCoefficient returns the coefficient d of Marder from the dimensionless cleaning
 * strength s = d dt sum_i 2/dl_i^2 on the mesh of the layout (s = 1/2 removes the grid scale
 * mode in one step, s <= 1 is stable).
 */
template<typename GridLayout>
double marderCoefficient(GridLayout const& layout, double const s, double const dt)
{
    double sum = 0.;
    for (auto const dl : layout.meshSize())
        sum += 2. / (dl * dl);
    return s / (dt * sum);
}

} // namespace PHARE::core

#endif
