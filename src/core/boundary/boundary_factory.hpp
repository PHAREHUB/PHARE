#ifndef PHARE_CORE_BOUNDARY_BOUNDARY_FACTORY
#define PHARE_CORE_BOUNDARY_BOUNDARY_FACTORY

#include "core/boundary/boundary.hpp"
#include "core/boundary/boundary_defs.hpp"
#include "core/boundary/boundary_inflow_compose.hpp"
#include "core/data/field/field_traits.hpp"
#include "core/data/grid/gridlayout_traits.hpp"
#include "core/numerics/boundary_condition/field_boundary_condition_factory.hpp"
#include "core/numerics/primite_conservative_converter/to_conservative_converter.hpp"
#include "core/numerics/thermo/thermo.hpp"

#include "initializer/data_provider.hpp"
#include "initializer/dict_utils.hpp"

#include <memory>
#include <stdexcept>
#include <vector>

namespace PHARE::core
{

/**
 * @brief Concept that detects whether a physical quantity type carries the conserved-variable set
 * required by super-magnetofast inflow boundary conditions (momentum vector @c rhoV and total
 * energy @c Etot). Satisfied by MHDQuantity, not by HybridQuantity.
 */
template<typename T>
concept HasInflowQuantities = requires {
    { T::Vector::rhoV };
    { T::Scalar::Etot };
};

/**
 * @brief Contains all the recipes to create a boundary object according to the desired
 * type of physical boundary (reflective, open, ...). It can extracts all the necessary data from
 * the input data dict associated to the boundary (value of physical quantities on the boundary for
 * an Inflow condition for instance), and create the right boundary conditions associated to each
 * physical quantity that requires one.
 *
 * @tparam PhysicalQuantityT The model category of physical quantities (MHDQuantity or
 * HybridQuantity).
 * @tparam FieldT The type for scalar fields.
 * @tparam GridLayoutT The type for the grid layout.
 */
template<typename PhysicalQuantityT, IsField FieldT, IsGridLayout GridLayoutT>
class BoundaryFactory
{
public:
    using boundary_type             = Boundary<PhysicalQuantityT, FieldT, GridLayoutT>;
    using boundary_ptr_type         = std::unique_ptr<boundary_type>;
    using scalar_quantity_list_type = std::vector<typename PhysicalQuantityT::Scalar>;
    using vector_quantity_list_type = std::vector<typename PhysicalQuantityT::Vector>;

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
     * @param thermo Optional thermodynamic model used by EOS-dependent boundary conditions
     *               (e.g. inflow). May be nullptr for boundary types that do not require an EOS.
     *
     * @return A unique pointer to the created @c Boundary object.
     */
    static boundary_ptr_type create(BoundaryLocation location, initializer::PHAREDict dict,
                                    scalar_quantity_list_type const& scalars,
                                    vector_quantity_list_type const& vectors,
                                    std::shared_ptr<Thermo> thermo = nullptr)
    {
        std::string typeName = dict["type"].to<std::string>();
        BoundaryType type    = getBoundaryTypeFromString(typeName);
        _model_menu_type const quantities{scalars, vectors};
        initializer::PHAREDict const data
            = (dict.contains("data")) ? dict["data"] : initializer::PHAREDict{};

        // initialize the boundary
        boundary_ptr_type boundary = std::make_unique<boundary_type>(type, location);

        // register the right boundary condition per physical quantity following the boundary type
        switch (type)
        {
            case BoundaryType::None: register_none_conditions_(boundary, quantities); break;
            case BoundaryType::Reflective:
                register_reflective_conditions_(boundary, data, quantities);
                break;
            case BoundaryType::SuperMagnetofastInflow:
                if constexpr (HasInflowQuantities<PhysicalQuantityT>)
                    register_super_magnetofast_inflow_conditions_(boundary, data, quantities,
                                                                  thermo);
                else
                    throw std::runtime_error(
                        "SuperMagnetofastInflow boundary type is not supported for this physical "
                        "model.");
                break;
            case BoundaryType::SuperMagnetofastOutflow:
            case BoundaryType::Open:
                if constexpr (HasInflowQuantities<PhysicalQuantityT>)
                    register_open_conditions_(boundary, data, quantities, thermo);
                else
                    throw std::runtime_error(
                        "Open boundary type is not supported for this physical "
                        "model.");
                break;
            case BoundaryType::FreePressureInflow:
                if constexpr (HasInflowQuantities<PhysicalQuantityT>)
                    register_free_pressure_inflow_conditions_(boundary, data, quantities, thermo);
                else
                    throw std::runtime_error(
                        "FreePressureInflow boundary type is not supported for this physical "
                        "model.");
                break;
            case BoundaryType::FixedPressureOutflow:
                if constexpr (HasInflowQuantities<PhysicalQuantityT>)
                    register_fixed_pressure_outflow_conditions_(boundary, data, quantities, thermo);
                else
                    throw std::runtime_error(
                        "FixedPressureOutflow boundary type is not supported for this physical "
                        "model.");
                break;
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

    using _space_time_function = initializer::SpaceTimeFunction<dimension>;

    /** @brief Whether a vector data field is given as a space/time function (vs a constant).
     * Used to drive a time-varying inflow (e.g. IMF turning). Reads the explicit boolean flag
     * "<key>_is_function" written by the Python layer, avoiding any cppdict variant
     * introspection (the pinned cppdict release has no Dict::is<T>()). */
    static bool isFunctionXYZ_(initializer::PHAREDict const& data, std::string const& key)
    {
        auto const flag = key + "_is_function";
        return data.contains(flag) && data[flag].template to<bool>();
    }

    /** @brief Lift a constant 3-vector to three constant-valued space-time functions. */
    static std::array<_space_time_function, 3> liftConst3_(std::array<double, 3> const& c)
    {
        return {inflow_compose::constFunction<dimension>(c[0]),
                inflow_compose::constFunction<dimension>(c[1]),
                inflow_compose::constFunction<dimension>(c[2])};
    }

    /** @brief A prescribable 3-vector as functions, whether given as constants or functions. */
    static std::array<_space_time_function, 3> vecAsFunctions_(initializer::PHAREDict const& data,
                                                               std::string const& key)
    {
        if (isFunctionXYZ_(data, key))
            return initializer::parseDimXYZType<_space_time_function, 3>(data, key);
        return liftConst3_(initializer::parseDimXYZType<double, 3>(data, key));
    }

    /** @brief Build the shared B sub-BC used by compound conditions (e.g.
     * TotalEnergyFromPressure) to provide ghost values of the prescribed total field.
     * It is a plain Dirichlet on the total field, used only inside the energy computation;
     * the magnetic field's own ghost condition is a separate divergence-free transverse
     * Dirichlet (see register_inflow_conditions_). */
    template<typename VecFieldT, typename VectorBcType>
    static std::shared_ptr<VectorBcType> make_inflow_B_bc_(std::array<double, 3> const& B)
    {
        return std::shared_ptr<VectorBcType>{
            FieldBoundaryConditionFactory::create<FieldBoundaryConditionType::Dirichlet, VecFieldT,
                                                  GridLayoutT>(B)};
    }

    /** @brief Time-varying overload: the prescribed total field B(t) is given as three
     * space-time functions, evaluated at @c ctx.time. */
    template<typename VecFieldT, typename VectorBcType>
    static std::shared_ptr<VectorBcType>
    make_inflow_B_bc_(std::array<_space_time_function, 3> const& Bfns)
    {
        return std::shared_ptr<VectorBcType>{
            FieldBoundaryConditionFactory::create<FieldBoundaryConditionType::Dirichlet, VecFieldT,
                                                  GridLayoutT>(Bfns)};
    }

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
                                                PHARE::initializer::PHAREDict const& data,
                                                _model_menu_type const& quantities)
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
                case (PhysicalQuantityT::Vector::B):
                    // Fill outside-domain B ghosts with a divergence-free transverse Neumann
                    // extrapolation of the interior field (Faraday runs on the interior box only,
                    // so the ghost B must be provided by this condition rather than CT).
                    boundary->template registerFieldCondition<
                        FieldBoundaryConditionType::DivergenceFreeTransverseNeumann>(quantity);
                    break;
                case (PhysicalQuantityT::Vector::J):
                    boundary->template registerFieldCondition<
                        FieldBoundaryConditionType::AntiSymmetric>(quantity);
                    break;
                case (PhysicalQuantityT::Vector::E):
                    boundary->template registerFieldCondition<
                        FieldBoundaryConditionType::AntiSymmetric>(quantity);
                    break;
                default:
                    boundary
                        ->template registerFieldCondition<FieldBoundaryConditionType::Symmetric>(
                            quantity);
                    break;
            }
        }
    }

    /** @brief Register boundary conditions to make an open boundary */
    static void register_open_conditions_(boundary_ptr_type& boundary,
                                          initializer::PHAREDict const& data,
                                          _model_menu_type const& quantities,
                                          std::shared_ptr<Thermo> thermo)
    {
        using VecFieldT    = VecField<FieldT, PhysicalQuantityT>;
        using ScalarBcType = IFieldBoundaryCondition<FieldT, GridLayoutT>;
        using VectorBcType = IFieldBoundaryCondition<VecFieldT, GridLayoutT>;

        auto rho_bc = std::shared_ptr<ScalarBcType>{
            FieldBoundaryConditionFactory::create<FieldBoundaryConditionType::Neumann, FieldT,
                                                  GridLayoutT>()};
        auto P_bc = std::shared_ptr<ScalarBcType>{
            FieldBoundaryConditionFactory::create<FieldBoundaryConditionType::Neumann, FieldT,
                                                  GridLayoutT>()};
        auto rhoV_bc = std::shared_ptr<VectorBcType>{
            FieldBoundaryConditionFactory::create<FieldBoundaryConditionType::Neumann, VecFieldT,
                                                  GridLayoutT>()};
        auto B_bc = std::shared_ptr<VectorBcType>{FieldBoundaryConditionFactory::create<
            FieldBoundaryConditionType::DivergenceFreeTransverseNeumann, VecFieldT, GridLayoutT>()};
        for (auto const quantity : quantities.scalars)
        {
            switch (quantity)
            {
                case (PhysicalQuantityT::Scalar::rho):
                    boundary->template registerFieldCondition<FieldBoundaryConditionType::Neumann>(
                        quantity);
                    break;
                case (PhysicalQuantityT::Scalar::Etot):
                    boundary->template registerFieldCondition<
                        FieldBoundaryConditionType::TotalEnergyFromPressure>(
                        quantity, rho_bc, rhoV_bc, B_bc, P_bc, thermo);
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
                case (PhysicalQuantityT::Vector::B):
                    boundary->template registerFieldCondition<
                        FieldBoundaryConditionType::DivergenceFreeTransverseNeumann>(quantity);
                    break;
                case (PhysicalQuantityT::Vector::E):
                    boundary->template registerFieldCondition<FieldBoundaryConditionType::None>(
                        quantity);
                    break;
                default:
                    boundary->template registerFieldCondition<FieldBoundaryConditionType::Neumann>(
                        quantity);
                    break;
            }
        }
    }

    /** @brief Shared builder for the two inflow variants (super-magnetofast and free-pressure).
     *
     *  They are identical except for the pressure sub-BC used by the energy reconstruction,
     *  which the caller supplies as @p P_bc: a prescribed Dirichlet pressure
     *  (super-magnetofast) or a zero-gradient Neumann (free-pressure). Everything else is
     *  common: Dirichlet ρ, the momentum composite ρv = ρ·v, the prescribed-B energy sub-BC,
     *  and a divergence-free transverse Dirichlet B ghost condition (E is left None). Each of
     *  density/velocity/B may independently be a constant or a space-time function; the
     *  all-constant fast path is the constant-lifted case of the space-time-function path.
     *
     *  Only available for physical quantity types carrying conserved MHD variables. */
    static void
    register_inflow_conditions_(boundary_ptr_type& boundary, initializer::PHAREDict const& data,
                                _model_menu_type const& quantities, std::shared_ptr<Thermo> thermo,
                                std::shared_ptr<IFieldBoundaryCondition<FieldT, GridLayoutT>> P_bc,
                                std::string const& bcName)
        requires HasInflowQuantities<PhysicalQuantityT>
    {
        if (!thermo)
            throw std::runtime_error("BoundaryFactory: a Thermo object is required for " + bcName
                                     + " boundaries but none was provided.");

        if (!data.contains("B"))
            throw std::runtime_error("BoundaryFactory: " + bcName
                                     + " requires the magnetic field 'B'.");

        using VecFieldT    = VecField<FieldT, PhysicalQuantityT>;
        using ScalarBcType = IFieldBoundaryCondition<FieldT, GridLayoutT>;
        using VectorBcType = IFieldBoundaryCondition<VecFieldT, GridLayoutT>;
        using STF          = _space_time_function;

        bool const rhoIsFn = isFunctionXYZ_(data, "density");
        bool const vIsFn   = isFunctionXYZ_(data, "velocity");
        bool const bIsFn   = isFunctionXYZ_(data, "B");

        // density read once; reused by rho_bc, the momentum composite, and the direct rho
        // condition. rhoFn is always a valid function (the prescribed one, or the constant
        // lifted) so the composite path never re-reads the dict.
        double const rhoDbl = rhoIsFn ? 0.0 : data["density"].template to<double>();
        STF const rhoFn     = rhoIsFn ? data["density"].template to<STF>()
                                      : inflow_compose::constFunction<dimension>(rhoDbl);

        // --- scalar Dirichlet sub/main BC: rho ---
        auto rho_bc
            = rhoIsFn
                  ? std::shared_ptr<ScalarBcType>{FieldBoundaryConditionFactory::create<
                        FieldBoundaryConditionType::Dirichlet, FieldT, GridLayoutT>(rhoFn)}
                  : std::shared_ptr<ScalarBcType>{
                        FieldBoundaryConditionFactory::create<FieldBoundaryConditionType::Dirichlet,
                                                              FieldT, GridLayoutT>(rhoDbl)};

        // --- momentum rhoV = rho * v: constant fast-path if rho and v both constant, else composed
        std::shared_ptr<VectorBcType> rhoV_bc;
        if (!rhoIsFn && !vIsFn)
        {
            auto const v = initializer::parseDimXYZType<double, 3>(data, "velocity");
            rhoV_bc      = std::shared_ptr<VectorBcType>{
                FieldBoundaryConditionFactory::create<FieldBoundaryConditionType::Dirichlet,
                                                           VecFieldT, GridLayoutT>(vToRhoV(rhoDbl, v))};
        }
        else
        {
            auto const vFns = vecAsFunctions_(data, "velocity");
            std::array<STF, 3> const rhoVfns
                = {inflow_compose::mulFunction<dimension>(rhoFn, vFns[0]),
                   inflow_compose::mulFunction<dimension>(rhoFn, vFns[1]),
                   inflow_compose::mulFunction<dimension>(rhoFn, vFns[2])};
            rhoV_bc = std::shared_ptr<VectorBcType>{
                FieldBoundaryConditionFactory::create<FieldBoundaryConditionType::Dirichlet,
                                                      VecFieldT, GridLayoutT>(rhoVfns)};
        }

        // --- magnetic field B: two conditions built from the same prescribed inflow field ---
        //  * B_bc: plain-Dirichlet total field, used only as a sub-BC of the energy
        //    reconstruction (TotalEnergyFromPressure) to provide ghost B values.
        //  * B_main: the field's own ghost condition, a divergence-free transverse Dirichlet.
        //    It sets the transverse ghost components to the prescribed field while keeping the
        //    normal face divergence free. Being value-prescribed (not extrapolated from the
        //    interior) it also produces valid ghosts at regrid/init, when the fine interior is
        //    not yet available, so no separate regrid fallback condition is needed.
        std::shared_ptr<VectorBcType> B_bc;
        std::shared_ptr<VectorBcType> B_main;
        if (!bIsFn)
        {
            auto const B = initializer::parseDimXYZType<double, 3>(data, "B");
            B_bc         = make_inflow_B_bc_<VecFieldT, VectorBcType>(B);
            B_main       = std::shared_ptr<VectorBcType>{FieldBoundaryConditionFactory::create<
                      FieldBoundaryConditionType::DivergenceFreeTransverseDirichlet, VecFieldT,
                      GridLayoutT>(B)};
        }
        else
        {
            auto const Bfns = vecAsFunctions_(data, "B");
            B_bc            = make_inflow_B_bc_<VecFieldT, VectorBcType>(Bfns);
            B_main          = std::shared_ptr<VectorBcType>{FieldBoundaryConditionFactory::create<
                         FieldBoundaryConditionType::DivergenceFreeTransverseDirichlet, VecFieldT,
                         GridLayoutT>(Bfns)};
        }

        // energy sub-BC needs rho/rhoV/B as ghost providers; the same rho_bc / rhoV_bc objects
        // also serve as the quantities' own main Dirichlet conditions (registered by pointer).
        for (auto const quantity : quantities.scalars)
        {
            switch (quantity)
            {
                case (PhysicalQuantityT::Scalar::rho):
                    boundary->registerFieldCondition(quantity, rho_bc);
                    break;
                case (PhysicalQuantityT::Scalar::Etot):
                    boundary->template registerFieldCondition<
                        FieldBoundaryConditionType::TotalEnergyFromPressure>(
                        quantity, rho_bc, rhoV_bc, B_bc, P_bc, thermo);
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
                case (PhysicalQuantityT::Vector::rhoV):
                    boundary->registerFieldCondition(quantity, rhoV_bc);
                    break;
                case (PhysicalQuantityT::Vector::B):
                    boundary->registerFieldCondition(quantity, B_main);
                    break;
                default:
                    boundary->template registerFieldCondition<FieldBoundaryConditionType::None>(
                        quantity);
                    break;
            }
        }
    }

    /** @brief Register boundary conditions to make a super-magnetofast inflow boundary.
     *
     *  Delegates to @c register_inflow_conditions_ with a prescribed (Dirichlet) pressure
     *  sub-BC: Etot is derived from @c data["pressure"] (constant or space-time function)
     *  through @c FieldTotalEnergyFromPressureBoundaryCondition. See that shared builder for
     *  the ρ / ρv / E / B handling. */
    static void register_super_magnetofast_inflow_conditions_(boundary_ptr_type& boundary,
                                                              initializer::PHAREDict const& data,
                                                              _model_menu_type const& quantities,
                                                              std::shared_ptr<Thermo> thermo)
        requires HasInflowQuantities<PhysicalQuantityT>
    {
        using ScalarBcType = IFieldBoundaryCondition<FieldT, GridLayoutT>;
        using STF          = _space_time_function;

        auto P_bc = isFunctionXYZ_(data, "pressure")
                        ? std::shared_ptr<ScalarBcType>{FieldBoundaryConditionFactory::create<
                              FieldBoundaryConditionType::Dirichlet, FieldT, GridLayoutT>(
                              data["pressure"].template to<STF>())}
                        : std::shared_ptr<ScalarBcType>{FieldBoundaryConditionFactory::create<
                              FieldBoundaryConditionType::Dirichlet, FieldT, GridLayoutT>(
                              data["pressure"].template to<double>())};

        register_inflow_conditions_(boundary, data, quantities, thermo, std::move(P_bc),
                                    "SuperMagnetofastInflow");
    }

    /**
     * @brief Register boundary conditions for a free-pressure inflow boundary.
     *
     * Delegates to @c register_inflow_conditions_ with a Neumann (zero-gradient) pressure
     * sub-BC instead of a prescribed value; the energy ghost values are derived from that
     * Neumann pressure via the EOS. See the shared builder for the ρ / ρv / E / B handling.
     */
    static void register_free_pressure_inflow_conditions_(boundary_ptr_type& boundary,
                                                          initializer::PHAREDict const& data,
                                                          _model_menu_type const& quantities,
                                                          std::shared_ptr<Thermo> thermo)
        requires HasInflowQuantities<PhysicalQuantityT>
    {
        using ScalarBcType = IFieldBoundaryCondition<FieldT, GridLayoutT>;

        auto P_bc = std::shared_ptr<ScalarBcType>{
            FieldBoundaryConditionFactory::create<FieldBoundaryConditionType::Neumann, FieldT,
                                                  GridLayoutT>()};

        register_inflow_conditions_(boundary, data, quantities, thermo, std::move(P_bc),
                                    "FreePressureInflow");
    }

    /**
     * @brief Register boundary conditions for a fixed-pressure outflow boundary.
     *
     * ρ, ρv, and B use Neumann (zero-gradient) conditions. The pressure uses a Dirichlet
     * condition and the total energy ghost values are derived from the resulting ghost
     * pressure via @c FieldTotalEnergyFromPressureBoundaryCondition.
     *
     * Only available for physical quantity types carrying conserved MHD variables.
     */
    static void register_fixed_pressure_outflow_conditions_(boundary_ptr_type& boundary,
                                                            initializer::PHAREDict const& data,
                                                            _model_menu_type const& quantities,
                                                            std::shared_ptr<Thermo> thermo)
        requires HasInflowQuantities<PhysicalQuantityT>
    {
        if (!thermo)
            throw std::runtime_error(
                "BoundaryFactory: a Thermo object is required for FixedPressureOutflow "
                "boundaries but none was provided.");

        double const pressure = data["pressure"].to<double>();

        using VecFieldT    = VecField<FieldT, PhysicalQuantityT>;
        using ScalarBcType = IFieldBoundaryCondition<FieldT, GridLayoutT>;
        using VectorBcType = IFieldBoundaryCondition<VecFieldT, GridLayoutT>;

        auto rho_bc = std::shared_ptr<ScalarBcType>{
            FieldBoundaryConditionFactory::create<FieldBoundaryConditionType::Neumann, FieldT,
                                                  GridLayoutT>()};
        auto P_bc = std::shared_ptr<ScalarBcType>{
            FieldBoundaryConditionFactory::create<FieldBoundaryConditionType::Dirichlet, FieldT,
                                                  GridLayoutT>(pressure)};
        auto rhoV_bc = std::shared_ptr<VectorBcType>{
            FieldBoundaryConditionFactory::create<FieldBoundaryConditionType::Neumann, VecFieldT,
                                                  GridLayoutT>()};
        auto B_bc = std::shared_ptr<VectorBcType>{FieldBoundaryConditionFactory::create<
            FieldBoundaryConditionType::DivergenceFreeTransverseNeumann, VecFieldT, GridLayoutT>()};

        for (auto const quantity : quantities.scalars)
        {
            switch (quantity)
            {
                case (PhysicalQuantityT::Scalar::rho):
                    boundary->template registerFieldCondition<FieldBoundaryConditionType::Neumann>(
                        quantity);
                    break;
                case (PhysicalQuantityT::Scalar::Etot):
                    boundary->template registerFieldCondition<
                        FieldBoundaryConditionType::TotalEnergyFromPressure>(
                        quantity, rho_bc, rhoV_bc, B_bc, P_bc, thermo);
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
                case (PhysicalQuantityT::Vector::rhoV):
                    boundary->template registerFieldCondition<FieldBoundaryConditionType::Neumann>(
                        quantity);
                    break;
                case (PhysicalQuantityT::Vector::B):
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
};

} // namespace PHARE::core

#endif // PHARE_CORE_BOUNDARY_BOUNDARY_FACTORY
