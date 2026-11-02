#ifndef PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_FIELD_BOUNDARY_CONDITION_ETOT1_FROM_INFLOW_STATE_HPP
#define PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_FIELD_BOUNDARY_CONDITION_ETOT1_FROM_INFLOW_STATE_HPP

#include "core/boundary/boundary_defs.hpp"
#include "core/data/vecfield/vecfield.hpp"
#include "core/models/quantities/mhd_quantities.hpp"
#include "core/numerics/boundary_condition/external_field_on_layout.hpp"
#include "core/numerics/boundary_condition/field_boundary_condition.hpp"
#include "core/numerics/primite_conservative_converter/mhd_conversion.hpp"
#include "core/utilities/box/box_spans.hpp"

#include <array>
#include <cstddef>
#include <memory>

namespace PHARE::core
{
/**
 * @brief Boundary condition for the reduced energy Etot1 that imposes the inflow state, the
 * magnetic part being built from the perturbation @f$B_1 = B_{in} - B_0@f$.
 *
 * The external field B0 is recomputed from the model's updater on the layout and at the time of
 * the fill, then projected to cell centers the way the to-conservative converter projects B1.
 *
 * @tparam FieldT Type of the scalar field.
 * @tparam GridLayoutT Grid layout configuration.
 */
template<typename FieldT, typename GridLayoutT>
class FieldBoundaryConditionEtot1FromInflowState
    : public IFieldBoundaryCondition<FieldT, GridLayoutT>
{
public:
    using Super                  = IFieldBoundaryCondition<FieldT, GridLayoutT>;
    using physical_quantity_type = Super::physical_quantity_type;
    using vecfield_type          = VecField<FieldT, physical_quantity_type>;
    using external_field_type    = ExternalFieldOnLayout<vecfield_type, GridLayoutT>;
    using updater_type           = external_field_type::updater_type;

    static constexpr std::size_t dimension = Super::dimension;

    FieldBoundaryConditionEtot1FromInflowState(double const gamma, double const rho,
                                               std::array<double, 3> const& v,
                                               std::array<double, 3> const& Btotal, double const P,
                                               std::shared_ptr<updater_type> updater)
        : gamma_{gamma}
        , rho_{rho}
        , v_{v}
        , Btotal_{Btotal}
        , P_{P}
        , externalField_{std::move(updater)}
    {
    }

    FieldBoundaryConditionType getType() const override
    {
        return FieldBoundaryConditionType::Etot1FromInflowState;
    }

    void apply(FieldT& Etot1, BoundaryLocation const /*boundaryLocation*/,
               Box<std::uint32_t, dimension> const& localGhostBox, GridLayoutT const& gridLayout,
               double const time) override
    {
        using Scalar = MHDQuantity::Scalar;

        auto const& B0 = externalField_.B0(gridLayout, time);

        for_each_point(
            gridLayout.toFieldBox(localGhostBox, Etot1.physicalQuantity()), [&](auto const& index) {
                auto const idx = index.template as<std::size_t>();
                auto const b1x
                    = Btotal_[0]
                      - GridLayoutT::template project<Scalar::Bx, Scalar::rho>(B0[0], idx);
                auto const b1y
                    = Btotal_[1]
                      - GridLayoutT::template project<Scalar::By, Scalar::rho>(B0[1], idx);
                auto const b1z
                    = Btotal_[2]
                      - GridLayoutT::template project<Scalar::Bz, Scalar::rho>(B0[2], idx);
                Etot1(index) = eosPToEtot1(gamma_, rho_, v_[0], v_[1], v_[2], b1x, b1y, b1z, P_);
            });
    }

private:
    double gamma_;
    double rho_;
    std::array<double, 3> v_;
    std::array<double, 3> Btotal_;
    double P_;
    external_field_type externalField_;
};

} // namespace PHARE::core
#endif // PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_FIELD_BOUNDARY_CONDITION_ETOT1_FROM_INFLOW_STATE_HPP
