#ifndef PHARE_CORE_BOUNDARY_BOUNDARY_FACTORY
#define PHARE_CORE_BOUNDARY_BOUNDARY_FACTORY

#include "core/boundary/boundary.hpp"
#include "core/boundary/boundary_defs.hpp"
#include "core/data/field/field_traits.hpp"
#include "core/models/quantities/mhd_quantities.hpp"
#include "core/numerics/primite_conservative_converter/to_conservative_converter.hpp"

#include "initializer/data_provider.hpp"
#include "initializer/dict_utils.hpp"

#include <array>
#include <concepts>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace PHARE::core
{

template<typename T>
concept IsMHDQuantity = std::same_as<T, MHDQuantity>;

/**
 * @brief Contains all the recipes to create a boundary object according to the desired
 * type of physical boundary (reflective, open, ...). It extracts all the necessary data from
 * the input data dict associated to the boundary (value of physical quantities on the boundary for
 * an Inflow condition for instance), and create the right boundary conditions associated to each
 * physical quantity that requires one.
 *
 * @tparam GridLayoutT The type for the grid layout.
 */
template<typename GridLayoutT, IsField FieldT>
class BoundaryFactory
{
public:
    using physical_quantity_type    = decltype(GridLayoutT::options.field_options)::Quantity;
    using field_type                = FieldT;
    using boundary_type             = Boundary<GridLayoutT, FieldT>;
    using boundary_ptr_type         = std::unique_ptr<boundary_type>;
    using scalar_quantity_list_type = std::vector<typename physical_quantity_type::Scalar>;
    using vector_quantity_list_type = std::vector<typename physical_quantity_type::Vector>;

    static constexpr std::size_t dimension = GridLayoutT::dimension;

    BoundaryFactory() = delete;

    /**
     * @brief Create a boundary with the type indicated in the input dict, and register to it all
     * corresponding field boundary conditions.
     *
     * @param location The location of the boundary.
     * @param dict Input dictionnary related to the boundary.
     * @param scalars Scalar quantities for which it is necessary to register a field boundary
     *                condition.
     * @param vectors Vector quantities for which it is necessary to register a field boundary
     *                condition.
     * @param gamma Heat capacity ratio, used to compute the inflow total energy.
     *
     * @return A unique pointer to the created @c Boundary object.
     */
    static boundary_ptr_type create(BoundaryLocation location, initializer::PHAREDict dict,
                                    scalar_quantity_list_type const& scalars,
                                    vector_quantity_list_type const& vectors,
                                    double const gamma = 0.0)
    {
        BoundaryType const type = cppdict::get_value(dict, "type", BoundaryType::None);
        _model_menu_type const quantities{scalars, vectors};

        // initialize the boundary
        boundary_ptr_type boundary = std::make_unique<boundary_type>(type, location);

        // register the right boundary condition per physical quantity following the boundary type
        switch (type)
        {
            case BoundaryType::None: register_none_conditions_(boundary, quantities); break;
            case BoundaryType::Reflective:
                register_reflective_conditions_(boundary, quantities);
                break;
            case BoundaryType::SuperMagnetofastInflow:
                register_inflow_conditions_(boundary, dict, quantities, gamma);
                break;
            case BoundaryType::Open: register_open_conditions_(boundary, quantities); break;
            default: throw std::runtime_error("Boundary type not implemented.");
        }
        return boundary;
    }

private:
    /** @brief Utility struct to group scalar and vector quantities together */
    struct _model_menu_type
    {
        scalar_quantity_list_type const& scalars;
        vector_quantity_list_type const& vectors;
    };


    /** @brief Register no-op (None) conditions so a "none" boundary leaves ghosts untouched
     * rather than falling through to another type or throwing "condition not found". */
    static void register_none_conditions_(boundary_ptr_type& boundary,
                                          _model_menu_type const& quantities)
    {
        for (auto const quantity : quantities.scalars)
            boundary->template registerFieldCondition<FieldBoundaryConditionType::None>(quantity);
        for (auto const quantity : quantities.vectors)
            boundary->template registerFieldCondition<FieldBoundaryConditionType::None>(quantity);
    }

    /** @brief Register boundary conditions to make a reflective boundary */
    static void register_reflective_conditions_(boundary_ptr_type& boundary,
                                                _model_menu_type const& quantities)
    {
        if constexpr (!IsMHDQuantity<physical_quantity_type>)
            throw std::runtime_error(
                "Reflective boundary type is only supported by the MHD model.");
        else
        {
            for (auto const quantity : quantities.scalars)
            {
                boundary->template registerFieldCondition<FieldBoundaryConditionType::Neumann>(
                    quantity);
            }
            for (auto const quantity : quantities.vectors)
            {
                switch (quantity)
                {
                    case (physical_quantity_type::Vector::B):
                        // Fill outside-domain B ghosts with a divergence-free transverse Neumann
                        // extrapolation of the interior field.
                        boundary->template registerFieldCondition<
                            FieldBoundaryConditionType::DivergenceFreeTransverseNeumann>(quantity);
                        break;
                    case (physical_quantity_type::Vector::E):
                        boundary->template registerFieldCondition<
                            FieldBoundaryConditionType::AntiSymmetric>(quantity);
                        break;
                    case (physical_quantity_type::Vector::rhoV):
                        boundary->template registerFieldCondition<
                            FieldBoundaryConditionType::Symmetric>(quantity);
                        break;
                    default:
                        boundary->template registerFieldCondition<FieldBoundaryConditionType::None>(
                            quantity);
                        break;
                }
            }
        }
    }

    /** @brief Register boundary conditions to make an open boundary */
    static void register_open_conditions_(boundary_ptr_type& boundary,
                                          _model_menu_type const& quantities)
    {
        if constexpr (!IsMHDQuantity<physical_quantity_type>)
            throw std::runtime_error("Open boundary type is only supported by the MHD model.");
        else
        {
            for (auto const quantity : quantities.scalars)
            {
                switch (quantity)
                {
                    case (physical_quantity_type::Scalar::rho):
                    case (physical_quantity_type::Scalar::Etot):
                        boundary
                            ->template registerFieldCondition<FieldBoundaryConditionType::Neumann>(
                                quantity);
                        break;
                    default:
                        boundary->template registerFieldCondition<FieldBoundaryConditionType::None>(
                            quantity);
                }
            }
            for (auto const quantity : quantities.vectors)
            {
                switch (quantity)
                {
                    case (physical_quantity_type::Vector::rhoV):
                        boundary
                            ->template registerFieldCondition<FieldBoundaryConditionType::Neumann>(
                                quantity);
                        break;
                    case (physical_quantity_type::Vector::B):
                        boundary->template registerFieldCondition<
                            FieldBoundaryConditionType::DivergenceFreeTransverseNeumann>(quantity);
                        break;
                    default:
                        boundary->template registerFieldCondition<FieldBoundaryConditionType::None>(
                            quantity);
                        break;
                }
            }
        }
    }

    /** @brief Register boundary conditions to make a super-magnetofast inflow boundary.
     *
     *  Density, momentum and total energy are imposed, the latter being computed from the
     *  prescribed inflow state (density, velocity, B, pressure). The magnetic field is imposed
     *  through the motional electric field E = -v x B on the boundary, and its ghosts are filled
     *  with the inflow B by a divergence-free transverse Dirichlet condition.
     */
    static void register_inflow_conditions_(boundary_ptr_type& boundary,
                                            initializer::PHAREDict const& data,
                                            _model_menu_type const& quantities, double const gamma)
    {
        if constexpr (!IsMHDQuantity<physical_quantity_type>)
            throw std::runtime_error(
                "SuperMagnetofastInflow boundary type is only supported by the MHD model.");
        else
        {
            if (!(gamma > 1.0))
                throw std::runtime_error(
                    "BoundaryFactory: a heat capacity ratio > 1 is required for "
                    "SuperMagnetofastInflow boundaries, got "
                    + std::to_string(gamma) + ".");

            if (!data.contains("B"))
                throw std::runtime_error(
                    "BoundaryFactory: SuperMagnetofastInflow requires the magnetic field 'B'.");

            auto const rho  = data["density"].template to<double>();
            auto const P    = data["pressure"].template to<double>();
            auto const v    = initializer::parseDimXYZType<double, 3>(data, "velocity");
            auto const B    = initializer::parseDimXYZType<double, 3>(data, "B");
            auto const Etot = eosPToEtot(gamma, rho, v[0], v[1], v[2], B[0], B[1], B[2], P);

            for (auto const quantity : quantities.scalars)
            {
                switch (quantity)
                {
                    case (physical_quantity_type::Scalar::rho):
                        boundary->template registerFieldCondition<
                            FieldBoundaryConditionType::Dirichlet>(
                            quantity, rho, DirichletExtrapolation::Constant);
                        break;
                    case (physical_quantity_type::Scalar::Etot):
                        boundary->template registerFieldCondition<
                            FieldBoundaryConditionType::Dirichlet>(
                            quantity, Etot, DirichletExtrapolation::Constant);
                        break;
                    default:
                        boundary->template registerFieldCondition<FieldBoundaryConditionType::None>(
                            quantity);
                        break;
                }
            }

            for (auto const quantity : quantities.vectors)
            {
                switch (quantity)
                {
                    case (physical_quantity_type::Vector::rhoV):
                        boundary->template registerFieldCondition<
                            FieldBoundaryConditionType::Dirichlet>(
                            quantity, vToRhoV(rho, v), DirichletExtrapolation::Constant);
                        break;
                    case (physical_quantity_type::Vector::B):
                        boundary->template registerFieldCondition<
                            FieldBoundaryConditionType::DivergenceFreeTransverseDirichlet>(
                            quantity, B, DirichletExtrapolation::Constant);
                        break;
                    case (physical_quantity_type::Vector::E):
                        boundary->template registerFieldCondition<
                            FieldBoundaryConditionType::Dirichlet>(
                            quantity,
                            std::array<double, 3>{v[2] * B[1] - v[1] * B[2],
                                                  v[0] * B[2] - v[2] * B[0],
                                                  v[1] * B[0] - v[0] * B[1]},
                            DirichletExtrapolation::Constant);
                        break;
                    default:
                        boundary->template registerFieldCondition<FieldBoundaryConditionType::None>(
                            quantity);
                        break;
                }
            }
        }
    }
};

} // namespace PHARE::core

#endif // PHARE_CORE_BOUNDARY_BOUNDARY_FACTORY
