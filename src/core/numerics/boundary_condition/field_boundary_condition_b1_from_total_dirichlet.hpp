#ifndef PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_FIELD_BOUNDARY_CONDITION_B1_FROM_TOTAL_DIRICHLET_HPP
#define PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_FIELD_BOUNDARY_CONDITION_B1_FROM_TOTAL_DIRICHLET_HPP

#include "core/boundary/boundary_defs.hpp"
#include "core/numerics/boundary_condition/divergence_free_transverse_common.hpp"
#include "core/numerics/boundary_condition/external_field_on_layout.hpp"
#include "core/numerics/boundary_condition/field_boundary_condition.hpp"
#include "core/utilities/box/box_spans.hpp"

#include <array>
#include <cstddef>
#include <memory>

namespace PHARE::core
{
/**
 * @brief Boundary condition for the perturbation B1 = B - B0 that imposes the total field B on
 * the tangential components and sets the normal component so that the numerical divergence of B1
 * is zero.
 *
 * The external field B0 is recomputed from the model's updater on the layout and at the time of
 * the fill, so the tangential ghosts hold @f$B_{in} - B_0@f$ exactly where B0 is stored.
 *
 * @warning Only valid for vector fields with the same centering as the magnetic field.
 *
 * @tparam VecFieldT Type of the vector field.
 * @tparam GridLayoutT Grid layout configuration.
 */
template<typename VecFieldT, typename GridLayoutT>
class FieldBoundaryConditionB1FromTotalDirichlet
    : public IFieldBoundaryCondition<VecFieldT, GridLayoutT>
{
public:
    using Super                = IFieldBoundaryCondition<VecFieldT, GridLayoutT>;
    using tensor_quantity_type = Super::tensor_quantity_type;
    using field_type           = Super::field_type;
    using value_type           = field_type::value_type;
    using external_field_type  = ExternalFieldOnLayout<VecFieldT, GridLayoutT>;
    using updater_type         = external_field_type::updater_type;

    static constexpr std::size_t dimension = Super::dimension;
    static constexpr std::size_t N         = Super::N;
    static_assert(N == 3,
                  "B1-from-total Dirichlet boundary condition only applies to vector fields.");

    FieldBoundaryConditionB1FromTotalDirichlet(std::array<value_type, N> const& Btotal,
                                               std::shared_ptr<updater_type> updater)
        : Btotal_{Btotal}
        , externalField_{std::move(updater)}
    {
    }

    FieldBoundaryConditionType getType() const override
    {
        return FieldBoundaryConditionType::B1FromTotalDirichlet;
    }

    void apply(VecFieldT& vecField, BoundaryLocation const boundaryLocation,
               Box<std::uint32_t, dimension> const& localGhostBox, GridLayoutT const& gridLayout,
               double const time) override
    {
        std::size_t const iNormal = static_cast<std::size_t>(getDirection(boundaryLocation));

        assert(gridLayout.centering(vecField) == gridLayout.centering(tensor_quantity_type::B));

        auto const& B0 = externalField_.B0(gridLayout, time);
        auto fields    = vecField.components();

        for_N<N>([&](auto iTransverse) {
            if (static_cast<std::size_t>(iTransverse) != iNormal)
            {
                field_type& B1c       = std::get<iTransverse>(fields);
                field_type const& B0c = B0[iTransverse];
                for_each_point(
                    gridLayout.toFieldBox(localGhostBox, B1c.physicalQuantity()),
                    [&](auto const& index) { B1c(index) = Btotal_[iTransverse] - B0c(index); });
            }
        });

        applyDivergenceFreeNormalComponent<dimension>(fields, iNormal, getSide(boundaryLocation),
                                                      gridLayout, localGhostBox);
    }

private:
    std::array<value_type, N> Btotal_;
    external_field_type externalField_;
};

} // namespace PHARE::core
#endif // PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_FIELD_BOUNDARY_CONDITION_B1_FROM_TOTAL_DIRICHLET_HPP
