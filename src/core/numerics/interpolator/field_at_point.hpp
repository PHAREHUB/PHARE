#ifndef PHARE_CORE_NUMERICS_INTERPOLATOR_FIELD_AT_POINT_HPP
#define PHARE_CORE_NUMERICS_INTERPOLATOR_FIELD_AT_POINT_HPP

#include <array>
#include <cmath>
#include <cstddef>
#include <tuple>

#include "core/numerics/interpolator/interpolator.hpp"
#include "core/utilities/point/point.hpp"

namespace PHARE::core
{

/**
 * @brief Interpolate a scalar field at an arbitrary physical-space point.
 *
 * This class performs the same polynomial interpolation as @ref Interpolator
 * but accepts a physical-space @ref Point instead of a particle's (iCell, delta)
 * pair.  It works for any mix of primal/dual centerings (node-, cell-, face-, or
 * edge-centred fields).
 *
 * Physical coordinates follow the same convention as @ref GridLayout::fieldNodeCoordinates:
 * the coordinate of AMR index @c i with spacing @c dx is simply @c i * dx (no
 * domain-origin offset is applied).
 *
 * @tparam dim         Spatial dimension (1, 2, or 3).
 * @tparam interpOrder Interpolation order (1, 2, or 3).
 */
template<std::size_t dim, std::size_t interpOrder>
class FieldAtPoint : private Interpolator<dim, interpOrder>
{
    using Base = Interpolator<dim, interpOrder>;

public:
    /**
     * @brief Return `true` iff @p point sits inside the (ghost-included) domain with
     * enough slack on every direction for this interpolator's stencil to fit.
     *
     * Currently the stencil uses **two consecutive grid values per direction**
     * (so 4 values in 2D, 8 in 3D). Per direction the slack required at each end of
     * the ghost box depends on the centering of the field being interpolated:
     *   - primal: stencil straddles integer node positions — no slack (margin = 0);
     *   - dual:   stencil straddles half-integer cell centres — 1 cell slack (margin = 1).
     */
    template<typename GridLayout>
    static bool pointIsInterpolable(GridLayout const& layout, Point<double, dim> const& point,
                                    std::array<QtyCentering, dim> const& centerings)
    {
        auto const& dx     = layout.meshSize();
        auto const& amrBox = layout.AMRBox();
        auto const nGhosts = static_cast<int>(layout.nbrGhosts());
        for (auto d = 0u; d < dim; ++d)
        {
            int const margin = (centerings[d] == QtyCentering::dual) ? 1 : 0;
            int const iCell  = static_cast<int>(std::floor(point[d] / dx[d]));
            if (iCell < amrBox.lower[d] - nGhosts + margin
                || iCell > amrBox.upper[d] + nGhosts - margin)
                return false;
        }
        return true;
    }

    /**
     * @brief Evaluate @p field at @p physPoint using order-interpOrder interpolation.
     *
     * The field centering is resolved at compile time from @p quantity.
     *
     * @tparam quantity   Physical quantity enum value (determines centering per dimension).
     * @param  layout     Grid layout of the patch that owns @p field.
     * @param  field      Field to interpolate from.
     * @param  physPoint  Target point in physical coordinates (AMR-index * meshSize).
     * @return            Interpolated value.
     */
    template<auto quantity, typename GridLayout, typename Field>
    double operator()(GridLayout const& layout, Field const& field,
                      Point<double, dim> const& physPoint)
    {
        auto const& dx = layout.meshSize();

        // Convert physical point to cell index and fractional delta.
        // Cell i spans [i*dx, (i+1)*dx], so iCell = floor(x/dx), delta = x/dx - iCell.
        Point<int, dim> iCell;
        std::array<double, dim> delta;
        for (auto d = 0u; d < dim; ++d)
        {
            double const normalizedPos = physPoint[d] / dx[d];
            iCell[d]                   = static_cast<int>(std::floor(normalizedPos));
            delta[d]                   = normalizedPos - iCell[d];
        }

        // Compute start indices and weights for both centerings.
        // MeshToParticle will select the appropriate one per field component dimension.
        this->template indexAndWeights_<QtyCentering, QtyCentering::dual>(layout, iCell, delta);
        this->template indexAndWeights_<QtyCentering, QtyCentering::primal>(layout, iCell, delta);

        auto const indexWeights = std::forward_as_tuple(this->dual_startIndex_, this->dual_weights_,
                                                        this->primal_startIndex_,
                                                        this->primal_weights_);

        return this->meshToParticle_.template operator()<GridLayout, quantity>(field, indexWeights);
    }

    /**
     * @brief Evaluate @p field at @p physPoint with the field centering resolved at runtime.
     *
     * Use this overload when the centering of the field is not known at compile time
     * (e.g. when iterating over heterogeneous field components inside a BC applier).
     * The centering is obtained from @p layout and the field's physical quantity.
     *
     * @param  layout    Grid layout of the patch that owns @p field.
     * @param  field     Field to interpolate from.
     * @param  physPoint Target point in physical coordinates (AMR-index * meshSize).
     * @return           Interpolated value.
     */
    template<typename GridLayout, typename Field>
    double operator()(GridLayout const& layout, Field const& field,
                      Point<double, dim> const& physPoint)
    {
        auto const centering = GridLayout::centering(field);

        auto const& dx = layout.meshSize();

        Point<int, dim> iCell;
        std::array<double, dim> delta;
        for (auto d = 0u; d < dim; ++d)
        {
            double const normalizedPos = physPoint[d] / dx[d];
            iCell[d]                   = static_cast<int>(std::floor(normalizedPos));
            delta[d]                   = normalizedPos - iCell[d];
        }

        this->template indexAndWeights_<QtyCentering, QtyCentering::dual>(layout, iCell, delta);
        this->template indexAndWeights_<QtyCentering, QtyCentering::primal>(layout, iCell, delta);

        // Select the correct start index and weights per dimension at runtime.
        auto startFor = [&](std::size_t d) -> std::uint32_t {
            return centering[d] == QtyCentering::primal ? this->primal_startIndex_[d]
                                                        : this->dual_startIndex_[d];
        };
        auto weightsFor = [&](std::size_t d) -> auto const& {
            return centering[d] == QtyCentering::primal ? this->primal_weights_[d]
                                                        : this->dual_weights_[d];
        };

        if constexpr (dim == 1)
        {
            auto const s0         = startFor(0);
            auto const& w0        = weightsFor(0);
            double fieldAtPoint   = 0.;
            for (auto i0 = 0u; i0 < w0.size(); ++i0)
                fieldAtPoint += field(s0 + i0) * w0[i0];
            return fieldAtPoint;
        }
        else if constexpr (dim == 2)
        {
            auto const s0 = startFor(0);
            auto const s1 = startFor(1);
            auto const& w0        = weightsFor(0);
            auto const& w1        = weightsFor(1);
            double fieldAtPoint   = 0.;
            for (auto i0 = 0u; i0 < w0.size(); ++i0)
                for (auto i1 = 0u; i1 < w1.size(); ++i1)
                    fieldAtPoint += field(s0 + i0, s1 + i1) * w0[i0] * w1[i1];
            return fieldAtPoint;
        }
        else
        {
            auto const s0 = startFor(0);
            auto const s1 = startFor(1);
            auto const s2 = startFor(2);
            auto const& w0        = weightsFor(0);
            auto const& w1        = weightsFor(1);
            auto const& w2        = weightsFor(2);
            double fieldAtPoint   = 0.;
            for (auto i0 = 0u; i0 < w0.size(); ++i0)
                for (auto i1 = 0u; i1 < w1.size(); ++i1)
                    for (auto i2 = 0u; i2 < w2.size(); ++i2)
                        fieldAtPoint += field(s0 + i0, s1 + i1, s2 + i2) * w0[i0] * w1[i1] * w2[i2];
            return fieldAtPoint;
        }
    }
};

} // namespace PHARE::core

#endif // PHARE_CORE_NUMERICS_INTERPOLATOR_FIELD_AT_POINT_HPP
