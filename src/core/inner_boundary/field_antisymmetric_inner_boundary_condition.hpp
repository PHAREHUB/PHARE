#ifndef PHARE_CORE_INNER_BOUNDARY_FIELD_ANTISYMMETRIC_INNER_BOUNDARY_CONDITION_HPP
#define PHARE_CORE_INNER_BOUNDARY_FIELD_ANTISYMMETRIC_INNER_BOUNDARY_CONDITION_HPP

#include "core/inner_boundary/field_inner_boundary_condition.hpp"
#include "core/inner_boundary/field_dirichlet_inner_boundary_condition.hpp"

namespace PHARE::core
{
template<typename ScalarOrTensorFieldT, typename GridLayoutT, typename PhysicalStateT>
class FieldAntisymmetricInnerBoundaryCondition
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


    FieldAntisymmetricInnerBoundaryCondition() = default;


    FieldInnerBoundaryConditionType getType() const override
    {
        return FieldInnerBoundaryConditionType::Antisymmetric;
    }

    void apply(ScalarOrTensorFieldT& scalarOrTensorField, GridLayoutT const& layout,
               inner_boundary_mesh_data_type const& boundaryMeshData,
               context_type const& ctx) override
    {
        // if scalar, fallback to dirichlet with null value
        if constexpr (is_scalar)
        {
            scalarDirichletCondition_.apply(scalarOrTensorField, layout, boundaryMeshData, ctx);
        }
        else
        {
        // handling of the vector case
        auto fields = scalarOrTensorField.components();

        for_N<N>([&](auto ic) {
            constexpr auto i       = ic();
            auto& currentField     = std::get<i>(fields);
            auto const centering   = GridLayoutT::centering(currentField);
            auto const& ghostElems = boundaryMeshData.getGhostDataFromCentering(centering);

            for (ghost_elem_data_type const& ghostElem : ghostElems)
            {
                // WARNING: when the mirror is not interpolable, the ghost cell is left
                // untouched. This may be the cause of issues — TBD. If so, a lower-order
                // interpolation could be applied instead. See GhostElemData::mirrorIsInterpolable.
                if (!ghostElem.mirrorIsInterpolable)
                    continue;

                // ghostElem.mirrorIsInterpolable was vetted only for component i's centering.
                // Antisymmetric BC interpolates ALL sibling components at the same mirror, so
                // every component j must also be interpolable for its own centering.
                bool allInterpolable = true;
                for_N<N>([&](auto jc) {
                    constexpr auto j = jc();
                    if constexpr (j != i)
                    {
                        auto const centering_j = GridLayoutT::centering(std::get<j>(fields));
                        if (!interpolator_type::pointIsInterpolable(layout, ghostElem.mirrorPoint,
                                                                    centering_j))
                            allInterpolable = false;
                    }
                });
                if (!allInterpolable)
                    continue;

                // get interpolated value of all components at the mirror point
                Point<value_type, N> interpolatedComponents;
                for_N<N>([&](auto jc) {
                    constexpr auto j = jc();
                    interpolatedComponents[j]
                        = this->interpolator_(layout, std::get<j>(fields), ghostElem.mirrorPoint);
                });

                // Extend the dim-dimensional boundary normal to N dimensions,
                // padding with zeros for components beyond the spatial dimension
                // (e.g. z-component is 0 in 2D since the boundary has no z-normal).
                Point<value_type, N> normal_N{};
                for_N<dimension>([&](auto kc) {
                    constexpr auto k = kc();
                    normal_N[k]      = static_cast<value_type>(ghostElem.normal[k]);
                });

                // compute normal and tangential parts of the interpolated vector with respect to
                // the boundary. For antisymmetric BC: keep normal (Neumann), reverse tangential
                // so that linear interpolation at boundary gives zero tangential.
                auto const v_n = normal_N * dot_product(interpolatedComponents, normal_N);
                auto const v_t = interpolatedComponents - v_n;

                // ghost cell value: normal component preserved (Neumann), tangential reversed
                auto const v                  = v_n - v_t;
                currentField(ghostElem.index) = v[i];
            }
        });
        } // end else (vector case)
    }

private:
    FieldDirichletInnerBoundaryCondition<field_type, GridLayoutT, PhysicalStateT>
        scalarDirichletCondition_;
};

} // namespace PHARE::core
#endif // PHARE_CORE_INNER_BOUNDARY_FIELD_ANTISYMMETRIC_INNER_BOUNDARY_CONDITION_HPP
