#ifndef PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_FIELD_BOUNDARY_CONDITION_DIRICHLET_HPP
#define PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_FIELD_BOUNDARY_CONDITION_DIRICHLET_HPP

#include "core/boundary/boundary_defs.hpp"
#include "core/data/grid/gridlayout.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/numerics/boundary_condition/field_boundary_condition.hpp"

#include <array>
#include <cstddef>

namespace PHARE::core
{
enum class DirichletExtrapolation { Constant = 0, Linear = 1 };

/**
 * @brief Dirichlet boundary condition for scalar and vector fields.
 *
 * Impose a value on the boundary by extrapolating the (tensor) field in the ghost cells, either
 * linearly through the boundary value or as a constant equal to it. The imposed value is a
 * constant, per component for a tensor field.
 *
 * @tparam ScalarOrTensorFieldT Type of the field or tensor field.
 * @tparam GridLayoutT Grid layout configuration.
 *
 */
template<typename ScalarOrTensorFieldT, typename GridLayoutT>
class FieldBoundaryConditionDirichlet
    : public IFieldBoundaryCondition<ScalarOrTensorFieldT, GridLayoutT>
{
public:
    using Super                = IFieldBoundaryCondition<ScalarOrTensorFieldT, GridLayoutT>;
    using tensor_quantity_type = Super::tensor_quantity_type;
    using field_type           = Super::field_type;
    using value_type           = field_type::value_type;

    static constexpr std::size_t dimension = Super::dimension;
    static constexpr std::size_t N         = Super::N;
    static constexpr bool is_scalar        = Super::is_scalar;

    FieldBoundaryConditionDirichlet() = default;

    FieldBoundaryConditionDirichlet(value_type value, DirichletExtrapolation extrapolation)
        requires(is_scalar)
        : value_{value}
        , extrapolation_{extrapolation} {};

    FieldBoundaryConditionDirichlet(std::array<value_type, N> value,
                                    DirichletExtrapolation extrapolation)
        : value_{value}
        , extrapolation_{extrapolation} {};

    FieldBoundaryConditionDirichlet(FieldBoundaryConditionDirichlet const&)            = default;
    FieldBoundaryConditionDirichlet& operator=(FieldBoundaryConditionDirichlet const&) = default;
    FieldBoundaryConditionDirichlet(FieldBoundaryConditionDirichlet&&)                 = default;
    FieldBoundaryConditionDirichlet& operator=(FieldBoundaryConditionDirichlet&&)      = default;

    virtual ~FieldBoundaryConditionDirichlet() = default;

    FieldBoundaryConditionType getType() const override
    {
        return FieldBoundaryConditionType::Dirichlet;
    }

    void apply(ScalarOrTensorFieldT& scalarOrTensorField, BoundaryLocation const boundaryLocation,
               Box<std::uint32_t, dimension> const& localGhostBox, GridLayoutT const& gridLayout,
               [[maybe_unused]] double const time) override
    {
        Direction const direction = getDirection(boundaryLocation);
        Side const side           = getSide(boundaryLocation);

        auto fields = Super::asComponentTuple(scalarOrTensorField);

        std::size_t const iDir = static_cast<std::size_t>(direction);

        for_N<N>([&](auto i) {
            field_type& field            = std::get<i>(fields);
            QtyCentering const centering = GridLayoutT::centering(field.physicalQuantity())[iDir];
            auto fieldBox = gridLayout.toFieldBox(localGhostBox, field.physicalQuantity());

            auto extrapolate = [&](auto const& index, value_type const v) {
                auto const mirrorIndex
                    = gridLayout.boundaryMirrored(direction, side, centering, index);
                field(index) = (extrapolation_ == DirichletExtrapolation::Constant
                                || mirrorIndex[iDir] == index[iDir])
                                   ? v
                                   : 2.0 * v - field(mirrorIndex);
            };

            for (auto const& index : fieldBox)
                extrapolate(index, value_[i]);
        });
    }

private:
    std::array<value_type, N> value_{0};
    DirichletExtrapolation extrapolation_{DirichletExtrapolation::Linear};

}; // class FieldBoundaryConditionDirichlet

} // namespace PHARE::core
#endif // PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_FIELD_BOUNDARY_CONDITION_DIRICHLET_HPP
