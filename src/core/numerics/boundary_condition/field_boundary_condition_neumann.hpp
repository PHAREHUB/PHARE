#ifndef PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_FIELD_BOUNDARY_CONDITION_NEUMANN_HPP
#define PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_FIELD_BOUNDARY_CONDITION_NEUMANN_HPP

#include "core/boundary/boundary_defs.hpp"
#include "core/data/grid/gridlayout.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/numerics/boundary_condition/field_boundary_condition.hpp"

#include <cstddef>
#include <tuple>

namespace PHARE::core
{
/**
 * @brief Neumann boundary condition implementation for fields and tensor fields.
 *
 * Implements a zero-gradient boundary condition by mirroring values from the physical domain
 * into the ghost regions.
 *
 * @tparam ScalarOrTensorFieldT Type of the field or tensor field.
 * @tparam GridLayoutT Grid layout configuration.
 *
 */
template<typename ScalarOrTensorFieldT, typename GridLayoutT>
class FieldBoundaryConditionNeumann
    : public IFieldBoundaryCondition<ScalarOrTensorFieldT, GridLayoutT>
{
public:
    using Super                = IFieldBoundaryCondition<ScalarOrTensorFieldT, GridLayoutT>;
    using tensor_quantity_type = Super::tensor_quantity_type;
    using field_type           = Super::field_type;

    static constexpr std::size_t dimension = Super::dimension;
    static constexpr std::size_t N         = Super::N;
    static constexpr bool is_scalar        = Super::is_scalar;

    FieldBoundaryConditionNeumann() = default;

    FieldBoundaryConditionNeumann(FieldBoundaryConditionNeumann const&)            = default;
    FieldBoundaryConditionNeumann& operator=(FieldBoundaryConditionNeumann const&) = default;
    FieldBoundaryConditionNeumann(FieldBoundaryConditionNeumann&&)                 = default;
    FieldBoundaryConditionNeumann& operator=(FieldBoundaryConditionNeumann&&)      = default;

    virtual ~FieldBoundaryConditionNeumann() = default;

    FieldBoundaryConditionType getType() const override
    {
        return FieldBoundaryConditionType::Neumann;
    }

    void apply(ScalarOrTensorFieldT& scalarOrTensorField, BoundaryLocation const boundaryLocation,
               Box<std::uint32_t, dimension> const& localGhostBox, GridLayoutT const& gridLayout,
               [[maybe_unused]] double const time) override
    {
        using Index               = Point<std::uint32_t, dimension>;
        Direction const direction = getDirection(boundaryLocation);
        Side const side           = getSide(boundaryLocation);

        auto fields = Super::asComponentTuple(scalarOrTensorField);

        for_N<N>([&](auto i) {
            field_type& field            = std::get<i>(fields);
            QtyCentering const centering = GridLayoutT::centering(
                field.physicalQuantity())[static_cast<std::size_t>(direction)];
            auto fieldBox = gridLayout.toFieldBox(localGhostBox, field.physicalQuantity());
            for (Index const& index : fieldBox)
            {
                Index mirrorIndex = gridLayout.boundaryMirrored(direction, side, centering, index);
                field(index)      = field(mirrorIndex);
            }
        });
    }
}; // class FieldBoundaryConditionNeumann

} // namespace PHARE::core
#endif // PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_FIELD_BOUNDARY_CONDITION_NEUMANN_HPP
