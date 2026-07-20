#ifndef PHARE_CORE_INNER_BOUNDARY_FIELD_ADAPTIVE_DIRICHLET_OR_NEUMANN_INNER_BOUNDARY_CONDITION_HPP
#define PHARE_CORE_INNER_BOUNDARY_FIELD_ADAPTIVE_DIRICHLET_OR_NEUMANN_INNER_BOUNDARY_CONDITION_HPP

#include "core/inner_boundary/field_inner_boundary_condition.hpp"
#include "core/mhd/mhd_quantities.hpp"

#include <array>

namespace PHARE::core
{
/**
 * @brief Per-ghost-element characteristic switch between Dirichlet and Neumann.
 *
 * For each ghost element the BC samples a target vector quantity (the @p criterion, e.g.
 * momentum rhoV) at the boundary *surface* point and projects it on the outward boundary
 * normal. The sign of that normal flux selects the condition:
 *
 *   - criterion·n > 0  (flux entering the fluid domain from the body) → Dirichlet:
 *                       ghost = 2*value - mirror  (prescribe a reservoir value)
 *   - criterion·n ≤ 0  (flux into the body, leaving the domain)       → Neumann:
 *                       ghost = mirror            (zero-gradient extrapolation)
 *
 * The criterion is interpolated at the boundary point — the projection of the ghost onto the
 * surface, midway between the ghost node and its mirror — because the in/out decision is a
 * physical condition on the surface, not at the reflected fluid point. Since
 * mirrorPoint = 2*project(ghost) - ghost, that surface point is 0.5*(ghostCoord + mirrorPoint).
 *
 * Only the flux *sign* uses the boundary point; the field value itself is still
 * extrapolated/reflected from the mirror point, exactly like the standalone Dirichlet/Neumann
 * BCs this composes.
 *
 * @note The criterion vector is resolved from @p criterion via PhysicalStateT::getVector(...),
 *       so the state type must expose that accessor (MHDState does).
 */
template<typename ScalarOrTensorFieldT, typename GridLayoutT, typename PhysicalStateT>
class FieldAdaptiveDirichletOrNeumannInnerBoundaryCondition
    : public FieldInnerBoundaryCondition<ScalarOrTensorFieldT, GridLayoutT, PhysicalStateT>
{
public:
    using Super = FieldInnerBoundaryCondition<ScalarOrTensorFieldT, GridLayoutT, PhysicalStateT>;

    using value_type                    = typename Super::field_type::value_type;
    using field_type                    = Super::field_type;
    using inner_boundary_mesh_data_type = Super::inner_boundary_mesh_data_type;
    using ghost_elem_data_type          = Super::ghost_elem_data_type;
    using inner_boundary_type           = Super::inner_boundary_type;
    using interpolator_type             = Super::interpolator_type;
    using state_type                    = Super::state_type;
    using context_type                  = Super::context_type;

    static constexpr size_t dimension = Super::dimension;
    static constexpr size_t N         = Super::N;
    static constexpr bool is_scalar   = Super::is_scalar;

    FieldAdaptiveDirichletOrNeumannInnerBoundaryCondition(MHDQuantity::Vector criterion,
                                                          value_type dirichletValue,
                                                          int priority = 0)
        : criterion_{criterion}
        , values_{dirichletValue}
        , priority_{priority}
    {
    }

    FieldAdaptiveDirichletOrNeumannInnerBoundaryCondition(MHDQuantity::Vector criterion,
                                                          std::array<value_type, N> dirichletValues,
                                                          int priority = 0)
        : criterion_{criterion}
        , values_{dirichletValues}
        , priority_{priority}
    {
    }

    FieldInnerBoundaryConditionType getType() const override
    {
        return FieldInnerBoundaryConditionType::AdaptiveDirichletOrNeumann;
    }

    // applyToMoments ordering: a value > 0 lets this BC run after lower-priority ones
    // (e.g. density after momentum). See InnerBoundaryManager::applyToMoments.
    int priority() const override { return priority_; }

    void apply(ScalarOrTensorFieldT& scalarOrTensorField, GridLayoutT const& layout,
               inner_boundary_mesh_data_type const& boundaryMeshData,
               context_type const& ctx) override
    {
        auto fields = [&]() {
            if constexpr (is_scalar)
                return std::make_tuple(scalarOrTensorField);
            else
                return scalarOrTensorField.components();
        }();

        // target vector quantity used as the in/out criterion (e.g. rhoV)
        auto criterionComps = ctx.statenew.getVector(criterion_).components();

        for_N<N>([&](auto i) {
            auto& field            = std::get<i>(fields);
            auto const centering   = GridLayoutT::centering(field);
            auto const& ghostElems = boundaryMeshData.getGhostDataFromCentering(centering);

            for (ghost_elem_data_type const& ghostElem : ghostElems)
            {
                // WARNING: when the mirror is not interpolable, the ghost cell is left
                // untouched. This may be the cause of issues — TBD. If so, a lower-order
                // interpolation could be applied instead. See GhostElemData::mirrorIsInterpolable.
                if (!ghostElem.mirrorIsInterpolable)
                    continue;

                // boundary (surface) point = projection of the ghost onto the boundary.
                // mirrorPoint = 2*project(ghost) - ghost  ⇒  project = 0.5*(ghost + mirror).
                auto const amrIdxU = layout.localToAMR(ghostElem.index);
                Point<int, dimension> amrIdx;
                for_N<dimension>([&](auto kc) { amrIdx[kc()] = static_cast<int>(amrIdxU[kc()]); });
                auto const ghostCoord = layout.fieldNodeCoordinates(field, amrIdx);

                Point<double, dimension> boundaryPoint;
                for_N<dimension>([&](auto kc) {
                    constexpr auto k = kc();
                    boundaryPoint[k] = 0.5 * (ghostCoord[k] + ghostElem.mirrorPoint[k]);
                });

                // criterion must be interpolable at the boundary point for every component we dot
                bool criterionInterpolable = true;
                for_N<dimension>([&](auto kc) {
                    constexpr auto k      = kc();
                    auto const& critComp  = std::get<k>(criterionComps);
                    auto const cCentering = GridLayoutT::centering(critComp);
                    if (!interpolator_type::pointIsInterpolable(layout, boundaryPoint, cCentering))
                        criterionInterpolable = false;
                });
                if (!criterionInterpolable)
                    continue;

                // normal flux of the criterion at the surface: sum_k criterion_k * n_k
                double fluxNormal = 0.0;
                for_N<dimension>([&](auto kc) {
                    constexpr auto k = kc();
                    fluxNormal
                        += this->interpolator_(layout, std::get<k>(criterionComps), boundaryPoint)
                           * ghostElem.normal[k];
                });

                double const mirrorValue
                    = this->interpolator_(layout, field, ghostElem.mirrorPoint);

                if (fluxNormal > 0.0)
                    // Dirichlet: linear reflection so the surface (half-point between ghost and
                    // mirror) equals the prescribed value values_[i].
                    field(ghostElem.index) = 2.0 * values_[i] - mirrorValue;
                else
                    field(ghostElem.index) = mirrorValue; // Neumann (zero-gradient)
            }
        });
    }

private:
    MHDQuantity::Vector criterion_;
    std::array<value_type, N> values_{0};
    int priority_{0};
};

} // namespace PHARE::core
#endif // PHARE_CORE_INNER_BOUNDARY_FIELD_ADAPTIVE_DIRICHLET_OR_NEUMANN_INNER_BOUNDARY_CONDITION_HPP
