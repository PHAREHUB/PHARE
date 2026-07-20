#ifndef PHARE_CORE_INNER_BOUNDARY_INNER_BOUNDARY_CONDITION_FACTORY_HPP
#define PHARE_CORE_INNER_BOUNDARY_INNER_BOUNDARY_CONDITION_FACTORY_HPP

#include "core/data/field/field_traits.hpp"
#include "core/data/grid/gridlayout_traits.hpp"
#include "core/data/vecfield/vecfield.hpp"
#include "core/inner_boundary/field_none_inner_boundary_condition.hpp"
#include "core/inner_boundary/field_antisymmetric_inner_boundary_condition.hpp"
#include "core/inner_boundary/field_neumann_inner_boundary_condition.hpp"
#include "core/inner_boundary/field_symmetric_inner_boundary_condition.hpp"
#include "core/inner_boundary/field_dirichlet_inner_boundary_condition.hpp"
#include "core/inner_boundary/field_adaptive_dirichlet_or_neumann_inner_boundary_condition.hpp"
#include "core/inner_boundary/field_ionospheric_convection_momentum_inner_boundary_condition.hpp"
#include "core/inner_boundary/field_total_energy_from_pressure_inner_boundary_condition.hpp"
#include "core/inner_boundary/field_inner_boundary_condition.hpp"
#include "core/inner_boundary/inner_boundary_geometry.hpp"
#include "core/numerics/thermo/thermo.hpp"

#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace PHARE::core
{

/**
 * @brief Taxonomy of supported inner-boundary physics types.
 *
 * Each value maps to a set of per-quantity boundary conditions with physics-based
 * rules, analogous to BoundaryType for outer boundaries.
 *
 * - **Reflective** — specular reflection.
 *     scalars: Neumann, B: Symmetric, rhoV: Symmetric, E: Antisymmetric,
 *     others: Neumann.
 * - **IonosphericConvection** — basic ionospheric inner boundary.
 *     E: Dirichlet(0), rhoV: Neumann, rho: adaptive Dirichlet/Neumann (criterion rhoV,
 *     prescribed density), Etot: TotalEnergyFromPressure with an adaptive Dirichlet/Neumann
 *     pressure (criterion rhoV, prescribed pressure). Moments are enforced momentum → density
 *     → energy via priorities.
 */
enum class InnerBoundaryConditionType { Reflective, IonosphericConvection };

inline InnerBoundaryConditionType getInnerBoundaryConditionTypeFromString(std::string const& name)
{
    static std::unordered_map<std::string, InnerBoundaryConditionType> const typeMap{
        {"reflective", InnerBoundaryConditionType::Reflective},
        {"ionospheric-convection", InnerBoundaryConditionType::IonosphericConvection},
    };
    auto it = typeMap.find(name);
    if (it == typeMap.end())
        throw std::runtime_error("Unknown inner boundary condition type: " + name);
    return it->second;
}


/**
 * @brief Factory that creates per-quantity inner boundary condition objects.
 *
 * Mirrors BoundaryFactory for outer boundaries. For a given InnerBoundaryConditionType
 * it registers the appropriate FieldInnerBoundaryCondition for every scalar and
 * vector quantity supplied by the caller.
 *
 * @tparam PhysicalQuantityT  Quantity traits (MHDQuantity, HybridQuantity, …).
 * @tparam FieldT             Scalar field type.
 * @tparam GridLayoutT        Grid layout type.
 * @tparam PhysicalStateT     Physical state type forwarded to BC appliers.
 */
template<typename PhysicalQuantityT, IsField FieldT, IsGridLayout GridLayoutT,
         typename PhysicalStateT>
class InnerBoundaryConditionFactory
{
public:
    using vecfield_type  = VecField<FieldT, PhysicalQuantityT>;
    using scalar_qty     = typename PhysicalQuantityT::Scalar;
    using vector_qty     = typename PhysicalQuantityT::Vector;
    using scalar_bc_type = FieldInnerBoundaryCondition<FieldT, GridLayoutT, PhysicalStateT>;
    using vector_bc_type = FieldInnerBoundaryCondition<vecfield_type, GridLayoutT, PhysicalStateT>;
    using scalar_bc_map_type = std::unordered_map<scalar_qty, std::unique_ptr<scalar_bc_type>>;
    using vector_bc_map_type = std::unordered_map<vector_qty, std::unique_ptr<vector_bc_type>>;
    using geometry_type      = InnerBoundaryGeometry<GridLayoutT::dimension>;

    InnerBoundaryConditionFactory() = delete;

    /**
     * @brief Populate @p scalarBCs and @p vectorBCs following the rules for @p type.
     *
     * @param type             Inner boundary condition type.
     * @param scalarQuantities Scalar quantities that need a BC.
     * @param vectorQuantities Vector quantities that need a BC.
     * @param scalarBCs        Output map: scalar quantity → BC (written in place).
     * @param vectorBCs        Output map: vector quantity → BC (written in place).
     */
    static void create(InnerBoundaryConditionType type,
                       std::vector<scalar_qty> const& scalarQuantities,
                       std::vector<vector_qty> const& vectorQuantities,
                       scalar_bc_map_type& scalarBCs, vector_bc_map_type& vectorBCs,
                       std::shared_ptr<Thermo> thermo, double prescribedDensity = 0.,
                       double prescribedPressure                      = 0.,
                       [[maybe_unused]] geometry_type const* geometry = nullptr)
    {
        switch (type)
        {
            case InnerBoundaryConditionType::Reflective:
                register_reflective_(scalarQuantities, vectorQuantities, scalarBCs, vectorBCs,
                                     std::move(thermo));
                break;
            case InnerBoundaryConditionType::IonosphericConvection:
                register_ionospheric_convection_(scalarQuantities, vectorQuantities, scalarBCs,
                                                 vectorBCs, std::move(thermo), prescribedDensity,
                                                 prescribedPressure);
                break;
            default: throw std::runtime_error("InnerBoundaryConditionFactory: unknown type");
        }
    }

private:
    template<typename ScalarOrTensorFieldT>
    using None = FieldNoneInnerBoundaryCondition<ScalarOrTensorFieldT, GridLayoutT, PhysicalStateT>;

    template<typename ScalarOrTensorFieldT>
    using Neumann
        = FieldNeumannInnerBoundaryCondition<ScalarOrTensorFieldT, GridLayoutT, PhysicalStateT>;

    template<typename ScalarOrTensorFieldT>
    using Symmetric
        = FieldSymmetricInnerBoundaryCondition<ScalarOrTensorFieldT, GridLayoutT, PhysicalStateT>;

    template<typename ScalarOrTensorFieldT>
    using Antisymmetric = FieldAntisymmetricInnerBoundaryCondition<ScalarOrTensorFieldT,
                                                                   GridLayoutT, PhysicalStateT>;

    template<typename ScalarOrTensorFieldT>
    using TotalEnergyFromPressure
        = FieldTotalEnergyFromPressureInnerBoundaryCondition<ScalarOrTensorFieldT, GridLayoutT,
                                                             PhysicalStateT>;

    template<typename ScalarOrTensorFieldT>
    using IonosphericConvectionMomentum
        = FieldIonosphericConvectionMomentumInnerBoundaryCondition<ScalarOrTensorFieldT,
                                                                   GridLayoutT, PhysicalStateT>;

    template<typename ScalarOrTensorFieldT>
    using Dirichlet
        = FieldDirichletInnerBoundaryCondition<ScalarOrTensorFieldT, GridLayoutT, PhysicalStateT>;

    template<typename ScalarOrTensorFieldT>
    using AdaptiveDirichletOrNeumann
        = FieldAdaptiveDirichletOrNeumannInnerBoundaryCondition<ScalarOrTensorFieldT, GridLayoutT,
                                                                PhysicalStateT>;

    /**
     * @brief Reflective body BC rules:
     *   Etot   → TotalEnergyFromPressure (zero-gradient on pressure, energy reconstructed)
     *   scalars → Neumann
     *   B       → None           (B handled via E + constrained transport, never written here)
     *   rhoV    → Symmetric      (no normal flow, free-slip)
     *   E       → Antisymmetric  (tangential E = 0)
     *   others  → Neumann
     */
    static void register_reflective_(std::vector<scalar_qty> const& scalars,
                                     std::vector<vector_qty> const& vectors,
                                     scalar_bc_map_type& scalarBCs, vector_bc_map_type& vectorBCs,
                                     std::shared_ptr<Thermo> thermo)
    {
        for (auto const qty : scalars)
        {
            if (qty == PhysicalQuantityT::Scalar::Etot)
                // Etot ghosts are derived from a Neumann pressure condition rather than mirrored
                // directly, to avoid spurious heating where strong gradients cross the boundary.
                scalarBCs[qty] = std::make_unique<TotalEnergyFromPressure<FieldT>>(
                    std::make_unique<Neumann<FieldT>>(), thermo);
            else
                scalarBCs[qty] = std::make_unique<Neumann<FieldT>>();
        }

        for (auto const qty : vectors)
        {
            switch (qty)
            {
                case PhysicalQuantityT::Vector::B:
                    vectorBCs[qty] = std::make_unique<None<vecfield_type>>();
                    break;
                case PhysicalQuantityT::Vector::rhoV:
                    vectorBCs[qty] = std::make_unique<Symmetric<vecfield_type>>();
                    break;
                case PhysicalQuantityT::Vector::E:
                    vectorBCs[qty] = std::make_unique<Antisymmetric<vecfield_type>>();
                    break;
                default: vectorBCs[qty] = std::make_unique<Neumann<vecfield_type>>(); break;
            }
        }
    }

    /**
     * @brief Basic ionospheric-convection body BC rules:
     *   E       → Dirichlet(0)   (no corotation: V→0 at the surface ⇒ E = -V×B = 0)
     *   rhoV    → Dirichlet       (momentum, enforced first)
     *   rho     → Dirichlet (constant extrapolation, prescribed density; enforced second)
     *   Etot   → TotalEnergyFromPressure wrapping a Dirichlet (constant extrapolation) pressure
     *             (prescribed p_in, held at the surface; enforced last)
     *   B       → None           (B handled via E + constrained transport, never written here)
     *   others  → Neumann
     *
     * Moment priorities (momentum 0 < density < energy) make applyToMoments run them in the
     * order momentum → density → energy.
     */
    static void register_ionospheric_convection_(
        std::vector<scalar_qty> const& scalars, std::vector<vector_qty> const& vectors,
        scalar_bc_map_type& scalarBCs, vector_bc_map_type& vectorBCs,
        std::shared_ptr<Thermo> thermo, double prescribedDensity, double prescribedPressure)
    {
        for (auto const qty : scalars)
        {
            if (qty == PhysicalQuantityT::Scalar::Etot)
            {
                std::unique_ptr<scalar_bc_type> bc_pressure = std::make_unique<Dirichlet<FieldT>>(
                    prescribedPressure, Dirichlet<FieldT>::ExtrapolationType::Constant);
                scalarBCs[qty] = std::make_unique<TotalEnergyFromPressure<FieldT>>(
                    std::move(bc_pressure), thermo);
            }
            else if (qty == PhysicalQuantityT::Scalar::rho)
            {
                std::unique_ptr<scalar_bc_type> density_bc = std::make_unique<Dirichlet<FieldT>>(
                    prescribedDensity, Dirichlet<FieldT>::ExtrapolationType::Constant);
                scalarBCs[qty] = std::move(density_bc);
            }
            else
                scalarBCs[qty] = std::make_unique<Neumann<FieldT>>();
        }

        for (auto const qty : vectors)
        {
            switch (qty)
            {
                case PhysicalQuantityT::Vector::B:
                    vectorBCs[qty] = std::make_unique<None<vecfield_type>>();
                    break;
                case PhysicalQuantityT::Vector::rhoV:
                    vectorBCs[qty] = std::make_unique<Dirichlet<vecfield_type>>(
                        0.0, Dirichlet<vecfield_type>::ExtrapolationType::Constant);
                    break;
                case PhysicalQuantityT::Vector::E:
                    vectorBCs[qty] = std::make_unique<Dirichlet<vecfield_type>>(
                        0.0, Dirichlet<vecfield_type>::ExtrapolationType::Constant);
                    break;
                default: vectorBCs[qty] = std::make_unique<Neumann<vecfield_type>>(); break;
            }
        }
    }
};

} // namespace PHARE::core

#endif // PHARE_CORE_INNER_BOUNDARY_INNER_BOUNDARY_CONDITION_FACTORY_HPP
