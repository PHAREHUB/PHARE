#ifndef PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_FIELD_DIRICHLET_BOUNDARY_CONDITION_HPP
#define PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_FIELD_DIRICHLET_BOUNDARY_CONDITION_HPP

#include "core/boundary/boundary_defs.hpp"
#include "core/data/grid/gridlayout.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/numerics/boundary_condition/field_boundary_condition.hpp"
#include "initializer/data_provider.hpp"

#include <array>
#include <cstddef>
#include <memory>
#include <optional>
#include <tuple>
#include <vector>

namespace PHARE::core
{
/**
 * @brief Dirichlet boundary condition for scalar and vector fields.
 *
 * Impose a value on the boundary by linearly extrapolating the (tensor) field in the ghost
 * cells. The imposed value can be either a compile-time constant (per component) or a user
 * function of space and time, evaluated at each ghost-box field node at the current time
 * (@c ctx.time). See @c initializer::SpaceTimeFunction.
 *
 * @tparam ScalarOrTensorFieldT Type of the field or tensor field.
 * @tparam GridLayoutT Grid layout configuration.
 *
 */
template<typename ScalarOrTensorFieldT, typename GridLayoutT>
class FieldDirichletBoundaryCondition
    : public IFieldBoundaryCondition<ScalarOrTensorFieldT, GridLayoutT>
{
public:
    using Super                    = IFieldBoundaryCondition<ScalarOrTensorFieldT, GridLayoutT>;
    using tensor_quantity_type     = Super::tensor_quantity_type;
    using field_type               = Super::field_type;
    using value_type               = field_type::value_type;

    static constexpr size_t dimension = Super::dimension;
    static constexpr size_t N         = Super::N;
    static constexpr bool is_scalar   = Super::is_scalar;

    FieldDirichletBoundaryCondition() = default;

    FieldDirichletBoundaryCondition(value_type value)
        requires(is_scalar)
        : value_{value} {};

    FieldDirichletBoundaryCondition(std::array<value_type, N> value)
        : value_{value} {};

    // Function-based Dirichlet value: f(x[,y[,z]], t). Scalar field overload.
    FieldDirichletBoundaryCondition(initializer::SpaceTimeFunction<dimension> fn)
        requires(is_scalar)
        : hasFn_{true}
    {
        fn_[0] = std::move(fn);
    }

    // Function-based Dirichlet value, one function per component (vector field).
    FieldDirichletBoundaryCondition(std::array<initializer::SpaceTimeFunction<dimension>, N> fns)
        : hasFn_{true}
    {
        for (std::size_t i = 0; i < N; ++i)
            fn_[i] = std::move(fns[i]);
    }

    FieldDirichletBoundaryCondition(FieldDirichletBoundaryCondition const&)            = default;
    FieldDirichletBoundaryCondition& operator=(FieldDirichletBoundaryCondition const&) = default;
    FieldDirichletBoundaryCondition(FieldDirichletBoundaryCondition&&)                 = default;
    FieldDirichletBoundaryCondition& operator=(FieldDirichletBoundaryCondition&&)      = default;

    virtual ~FieldDirichletBoundaryCondition() = default;

    FieldBoundaryConditionType getType() const override
    {
        return FieldBoundaryConditionType::Dirichlet;
    }

    void apply(ScalarOrTensorFieldT& scalarOrTensorField,
               BoundaryLocation const boundaryLocation,
               Box<std::uint32_t, dimension> const& localGhostBox, GridLayoutT const& gridLayout,
               Super::boundary_condition_context_type const& ctx) override
    {
        Direction const direction = getDirection(boundaryLocation);
        Side const side           = getSide(boundaryLocation);

        if (static_cast<size_t>(direction) >= dimension)
            return;

        auto fields = [&]() {
            if constexpr (is_scalar)
                return std::make_tuple(scalarOrTensorField);
            else
                return scalarOrTensorField.components();
        }();

        size_t const iDir = static_cast<size_t>(direction);

        for_N<N>([&](auto i) {
            field_type& field = std::get<i>(fields);
            QtyCentering const centering
                = GridLayoutT::centering(field.physicalQuantity())[iDir];
            auto fieldBox = gridLayout.toFieldBox(localGhostBox, field.physicalQuantity());

            auto extrapolate = [&](_index_type const& index, value_type const v) {
                _index_type mirrorIndex
                    = gridLayout.boundaryMirrored(direction, side, centering, index);
                field(index) = (mirrorIndex[iDir] == index[iDir]) ? v
                                                                  : 2.0 * v - field(mirrorIndex);
            };

            if (hasFn_)
            {
                // Evaluate the user function at every ghost-box field node, at ctx.time.
                // Local ghost node 0 maps to the (centering-aware) lower AMR ghost index,
                // so the global physical coordinate of a local node is taken from
                // fieldNodeCoordinates(field, ghostAMRlower + localIndex).
                auto const ghostAMRlower = gridLayout.AMRGhostBoxFor(field).lower;
                tuple_fixed_type<std::vector<double>, dimension> coords{};
                std::vector<_index_type> nodes;

                for (_index_type const& index : fieldBox)
                {
                    Point<int, dimension> amrIndex;
                    for (std::size_t d = 0; d < dimension; ++d)
                        amrIndex[d] = ghostAMRlower[d] + static_cast<int>(index[d]);
                    auto const xyz = gridLayout.fieldNodeCoordinates(field, amrIndex);
                    for_N<dimension>([&](auto d) { std::get<d>(coords).push_back(xyz[d]); });
                    nodes.push_back(index);
                }

                std::shared_ptr<Span<double>> valuesPtr = std::apply(
                    [&](auto const&... xs) { return (*fn_[i])(xs..., ctx.time); }, coords);
                Span<double> const& values = *valuesPtr;

                for (std::size_t k = 0; k < nodes.size(); ++k)
                    extrapolate(nodes[k], values[k]);
            }
            else
            {
                for (_index_type const& index : fieldBox)
                    extrapolate(index, value_[i]);
            }
        });
    }

private:
    using _index_type = Point<std::uint32_t, dimension>;

    std::array<value_type, N> value_{0};
    std::array<std::optional<initializer::SpaceTimeFunction<dimension>>, N> fn_{};
    bool hasFn_ = false;

}; // class FieldDirichletBoundaryCondition

} // namespace PHARE::core
#endif // PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_FIELD_DIRICHLET_BOUNDARY_CONDITION_HPP
