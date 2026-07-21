#ifndef PHARE_CORE_NUMERICS_FIELD_BOUNDARY_CONDITION_FACTORY
#define PHARE_CORE_NUMERICS_FIELD_BOUNDARY_CONDITION_FACTORY

#include "core/data/tensorfield/tensorfield_traits.hpp"
#include "core/data/vecfield/vecfield_traits.hpp"
#include "core/numerics/boundary_condition/field_antisymmetric_boundary_condition.hpp"
#include "core/numerics/boundary_condition/field_boundary_condition.hpp"
#include "core/numerics/boundary_condition/field_dirichlet_boundary_condition.hpp"
#include "core/numerics/boundary_condition/field_divergence_free_transverse_dirichlet_boundary_condition.hpp"
#include "core/numerics/boundary_condition/field_divergence_free_transverse_neumann_boundary_condition.hpp"
#include "core/numerics/boundary_condition/field_neumann_boundary_condition.hpp"
#include "core/numerics/boundary_condition/field_none_boundary_condition.hpp"
#include "core/numerics/boundary_condition/field_symmetric_boundary_condition.hpp"
#include "core/numerics/boundary_condition/field_total_energy_from_pressure_boundary_condition.hpp"

#include <memory>
#include <stdexcept>

namespace PHARE::core
{
/**
 * @brief Factory for creating field boundary condition objects.
 *
 */
class FieldBoundaryConditionFactory
{
public:
    /**
     * @brief Main function to create field boundary conditions
     *
     * It passes the required arguments for the constructor of the actual field boundary condition
     * by perfect forwarding.
     *
     * @tparam type The value of the enum @c FieldBoundaryConditionType.
     * @tparam ScalarOrTensorFieldT Field or TensorField.
     * @tparam GridLayoutT The type of grid layout.
     * @tparam Args The types of the arguments for the constructor of the boundary condition
     * corresponding to @p type.
     *
     * @param args Arguments passed by perfect forwarding to the boundary condition constructor
     * corresponding to @p type.
     *
     * @return A unique pointer to the created field boundary condition.
     */
    template<FieldBoundaryConditionType type, IsScalarOrTensorField ScalarOrTensorFieldT,
             IsGridLayout GridLayoutT, typename... Args>
    static std::unique_ptr<IFieldBoundaryCondition<ScalarOrTensorFieldT, GridLayoutT>>
    create(Args&&... args)
    {
        if constexpr (type == FieldBoundaryConditionType::None)
        {
            return std::make_unique<FieldNoneBoundaryCondition<ScalarOrTensorFieldT, GridLayoutT>>(
                std::forward<Args>(args)...);
        }
        else if constexpr (type == FieldBoundaryConditionType::Neumann)
        {
            return std::make_unique<
                FieldNeumannBoundaryCondition<ScalarOrTensorFieldT, GridLayoutT>>(
                std::forward<Args>(args)...);
        }
        else if constexpr (type == FieldBoundaryConditionType::Dirichlet)
        {
            return std::make_unique<
                FieldDirichletBoundaryCondition<ScalarOrTensorFieldT, GridLayoutT>>(
                std::forward<Args>(args)...);
        }
        else if constexpr (type == FieldBoundaryConditionType::Symmetric)
        {
            return std::make_unique<
                FieldSymmetricBoundaryCondition<ScalarOrTensorFieldT, GridLayoutT>>(
                std::forward<Args>(args)...);
        }
        else if constexpr (type == FieldBoundaryConditionType::AntiSymmetric)
        {
            return std::make_unique<
                FieldAntiSymmetricBoundaryCondition<ScalarOrTensorFieldT, GridLayoutT>>(
                std::forward<Args>(args)...);
        }
        else if constexpr (type == FieldBoundaryConditionType::DivergenceFreeTransverseNeumann)
        {
            if constexpr (IsVecField<ScalarOrTensorFieldT>)
            {
                return std::make_unique<FieldDivergenceFreeTransverseNeumannBoundaryCondition<
                    ScalarOrTensorFieldT, GridLayoutT>>(std::forward<Args>(args)...);
            }
            else
            {
                throw std::runtime_error("Divergence-free transverse Neumann condition only "
                                         "applies to vector fields.");
            }
        }
        else if constexpr (type == FieldBoundaryConditionType::DivergenceFreeTransverseDirichlet)
        {
            if constexpr (IsVecField<ScalarOrTensorFieldT>)
            {
                return std::make_unique<FieldDivergenceFreeTransverseDirichletBoundaryCondition<
                    ScalarOrTensorFieldT, GridLayoutT>>(std::forward<Args>(args)...);
            }
            else
            {
                throw std::runtime_error("Divergence-free transverse Dirichlet condition only "
                                         "applies to vector fields.");
            }
        }
        else if constexpr (type == FieldBoundaryConditionType::TotalEnergyFromPressure)
        {
            if constexpr (IsField<ScalarOrTensorFieldT>)
            {
                return std::make_unique<FieldTotalEnergyFromPressureBoundaryCondition<
                    ScalarOrTensorFieldT, GridLayoutT>>(std::forward<Args>(args)...);
            }
            else
            {
                throw std::runtime_error(
                    "TotalEnergyFromPressure condition only applies to scalar fields.");
            }
        }
        else
        {
            // static_assert(false, "Unhandled FieldBoundaryConditionType");
            throw std::runtime_error("Unhandled FieldBoundaryConditionType");
        };
    }
};

} // namespace PHARE::core

#endif // PHARE_CORE_NUMERICS_FIELD_BOUNDARY_CONDITION_FACTORY
