#ifndef PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_FIELD_NEUMANN_BOUNDARY_CONDITION_HPP
#define PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_FIELD_NEUMANN_BOUNDARY_CONDITION_HPP

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
class FieldNeumannBoundaryCondition
    : public IFieldBoundaryCondition<ScalarOrTensorFieldT, GridLayoutT>
{
public:
    using Super                = IFieldBoundaryCondition<ScalarOrTensorFieldT, GridLayoutT>;
    using tensor_quantity_type = Super::tensor_quantity_type;
    using field_type           = Super::field_type;

    static constexpr size_t dimension = Super::dimension;
    static constexpr size_t N         = Super::N;
    static constexpr bool is_scalar   = Super::is_scalar;

    FieldNeumannBoundaryCondition() = default;

    FieldNeumannBoundaryCondition(FieldNeumannBoundaryCondition const&)            = default;
    FieldNeumannBoundaryCondition& operator=(FieldNeumannBoundaryCondition const&) = default;
    FieldNeumannBoundaryCondition(FieldNeumannBoundaryCondition&&)                 = default;
    FieldNeumannBoundaryCondition& operator=(FieldNeumannBoundaryCondition&&)      = default;

    virtual ~FieldNeumannBoundaryCondition() = default;

    FieldBoundaryConditionType getType() const override
    {
        return FieldBoundaryConditionType::Neumann;
    }

    void apply(ScalarOrTensorFieldT& scalarOrTensorField,
               BoundaryLocation const boundaryLocation,
               Box<std::uint32_t, dimension> const& localGhostBox, GridLayoutT const& gridLayout,
               [[maybe_unused]] Super::boundary_condition_context_type const& ctx) override
    {
        using Index             = Point<std::uint32_t, dimension>;
        Direction const direction = getDirection(boundaryLocation);
        Side const side           = getSide(boundaryLocation);

        auto fields = [&]() {
            if constexpr (is_scalar)
                return std::make_tuple(scalarOrTensorField);
            else
                return scalarOrTensorField.components();
        }();

        for_N<N>([&](auto i) {
            field_type& field = std::get<i>(fields);
            QtyCentering const centering
                = GridLayoutT::centering(field.physicalQuantity())[static_cast<size_t>(direction)];
            auto fieldBox = gridLayout.toFieldBox(localGhostBox, field.physicalQuantity());
            for (Index const& index : fieldBox)
            {
                Index mirrorIndex = gridLayout.boundaryMirrored(direction, side, centering, index);
                field(index)      = field(mirrorIndex);
            }
        });
    }
}; // class FieldNeumannBoundaryCondition

} // namespace PHARE::core
#endif // PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_FIELD_NEUMANN_BOUNDARY_CONDITION_HPP
