#ifndef PHARE_CORE_INNER_BOUNDARY_INNER_BOUNDARY_MANAGER_HPP
#define PHARE_CORE_INNER_BOUNDARY_INNER_BOUNDARY_MANAGER_HPP

#include "core/data/field/field_traits.hpp"
#include "core/data/grid/gridlayout_traits.hpp"
#include "core/data/vecfield/vecfield.hpp"
#include "core/inner_boundary/inner_bc_context.hpp"
#include "core/inner_boundary/inner_boundary_condition_factory.hpp"
#include "core/inner_boundary/inner_boundary_factory.hpp"
#include "core/inner_boundary/inner_boundary_geometry.hpp"
#include "core/inner_boundary/inner_boundary_mesh_classifier.hpp"
#include "core/inner_boundary/inner_boundary_mesh_data.hpp"
#include "core/inner_boundary/inner_boundary_defs.hpp"

#include "core/numerics/thermo/thermo.hpp"
#include "core/numerics/primite_conservative_converter/conversion_utils.hpp"

#include "initializer/data_provider.hpp"

#include <algorithm>
#include <array>
#include <functional>
#include <memory>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

namespace PHARE::core
{

/**
 * @brief Manages the lifecycle and application of inner boundary conditions.
 *
 * Owns the embedded boundary geometry, the per-patch mesh classification data,
 * and the per-quantity BC objects. Exposes:
 *   - classify(layout)                   — fill mesh data (call once per patch setup)
 *   - applyBC(field, layout, ctx)        — enforce BC on a scalar quantity
 *   - applyBC(vecfield, layout, ctx)     — enforce BC on a vector quantity
 *   - ResourcesUser interface            — expose mesh data to SAMRAI resource manager
 *
 * @tparam PhysicalQuantityT  Quantity traits (MHDQuantity, HybridQuantity, …).
 * @tparam FieldT             Scalar field type.
 * @tparam GridLayoutT        Grid layout type.
 * @tparam PhysicalStateT     Physical state type forwarded to BC appliers.
 */
template<typename PhysicalQuantityT, IsField FieldT, IsGridLayout GridLayoutT,
         typename PhysicalStateT>
class InnerBoundaryManager
{
public:
    static constexpr std::size_t dimension = GridLayoutT::dimension;

    using vecfield_type   = VecField<FieldT, PhysicalQuantityT>;
    using scalar_qty      = PhysicalQuantityT::Scalar;
    using vector_qty      = PhysicalQuantityT::Vector;
    using geometry_type   = InnerBoundaryGeometry<dimension>;
    using mesh_data_type  = InnerBoundaryMeshData<dimension, PhysicalQuantityT>;
    using classifier_type = InnerBoundaryMeshClassifier<dimension, GridLayoutT, PhysicalQuantityT>;
    using context_type    = InnerBCContext<PhysicalStateT>;
    using scalar_bc_type  = FieldInnerBoundaryCondition<FieldT, GridLayoutT, PhysicalStateT>;
    using vector_bc_type  = FieldInnerBoundaryCondition<vecfield_type, GridLayoutT, PhysicalStateT>;
    using scalar_bc_map_type = std::unordered_map<scalar_qty, std::unique_ptr<scalar_bc_type>>;
    using vector_bc_map_type = std::unordered_map<vector_qty, std::unique_ptr<vector_bc_type>>;
    using factory_type
        = InnerBoundaryConditionFactory<PhysicalQuantityT, FieldT, GridLayoutT, PhysicalStateT>;

    InnerBoundaryManager() = delete;

    /**
     * @brief Prescribed "safe state" pinned into inactive (inside-body) cells.
     *
     * Primitive values; momentum (rhoV) and total energy (Etot) are derived from these.
     * Defaults reproduce the historical hardcoded reset (rho=1, P=1, V=0, B=0).
     */
    struct SafeStateValues
    {
        double density{1.};
        double pressure{1.};
        std::array<double, 3> velocity{0., 0., 0.};
        std::array<double, 3> B{0., 0., 0.};
    };

    /**
     * @brief Construct from pre-built geometry and BC condition type.
     *
     * @param geometry         Embedded boundary geometry (transferred ownership).
     * @param conditionType    Physics preset for per-quantity BC assignment.
     * @param scalarQuantities Scalar quantities that need a BC.
     * @param vectorQuantities Vector quantities that need a BC.
     * @param safe             Prescribed safe state for inactive cells (see setSafeState).
     */
    InnerBoundaryManager(std::unique_ptr<geometry_type> geometry,
                         InnerBoundaryConditionType conditionType,
                         std::vector<scalar_qty> const& scalarQuantities,
                         std::vector<vector_qty> const& vectorQuantities,
                         std::shared_ptr<Thermo> thermo, double prescribedDensity = 0.,
                         double prescribedPressure = 0., SafeStateValues const& safe = {})
        : geometry_{std::move(geometry)}
        , meshData_{geometry_->name()}
        , thermo_{std::move(thermo)}
    {
        factory_type::create(conditionType, scalarQuantities, vectorQuantities, scalarBCs_,
                             vectorBCs_, thermo_, prescribedDensity, prescribedPressure,
                             geometry_.get());
        buildSafeState_(safe);
    }

    /**
     * @brief Static factory: create a manager from a simulation dict, or return nullptr.
     *
     * Returns nullptr when the dict contains no "inner_boundary" key, matching
     * the behaviour of InnerBoundaryFactory::create().
     */
    static std::unique_ptr<InnerBoundaryManager>
    create(initializer::PHAREDict const& dict, std::vector<scalar_qty> const& scalarQuantities,
           std::vector<vector_qty> const& vectorQuantities, std::shared_ptr<Thermo> thermo)
    {
        if (!dict.contains("inner_boundary"))
            return nullptr;

        auto geometry = InnerBoundaryFactory<dimension>::create(dict);
        if (!geometry)
            return nullptr;

        auto const& ibDict  = dict["inner_boundary"];
        auto const typeName = ibDict["condition_type"].template to<std::string>();
        auto const condType = getInnerBoundaryConditionTypeFromString(typeName);

        // prescribed reservoir values (used by types that impose Dirichlet moments, e.g.
        // ionospheric-convection); absent for types that do not need them.
        auto const prescribedDensity  = cppdict::get_value(ibDict, "density", 0.);
        auto const prescribedPressure = cppdict::get_value(ibDict, "pressure", 0.);

        // safe state pinned into inactive cells; optional, defaults reproduce the historical
        // hardcoded reset (rho=1, P=1, V=0, B=0).
        SafeStateValues safe{};
        if (ibDict.contains("inactive_safe_state"))
        {
            auto const& ss = ibDict["inactive_safe_state"];
            safe.density   = cppdict::get_value(ss, "density", 1.0);
            safe.pressure  = cppdict::get_value(ss, "pressure", 1.0);
            safe.velocity  = {cppdict::get_value(ss, "velocity/x", 0.0),
                              cppdict::get_value(ss, "velocity/y", 0.0),
                              cppdict::get_value(ss, "velocity/z", 0.0)};
            safe.B         = {cppdict::get_value(ss, "B/x", 0.0),
                              cppdict::get_value(ss, "B/y", 0.0),
                              cppdict::get_value(ss, "B/z", 0.0)};
        }

        return std::make_unique<InnerBoundaryManager>(std::move(geometry), condType,
                                                      scalarQuantities, vectorQuantities,
                                                      std::move(thermo), prescribedDensity,
                                                      prescribedPressure, safe);
    }

    // -------------------------------------------------------------------------
    //  Classification
    // -------------------------------------------------------------------------

    /**
     * @brief Classify mesh elements relative to the embedded boundary.
     *
     * Must be called after setBuffer() has been invoked on the mesh data
     * (i.e., after SAMRAI has allocated the patch data).
     *
     * @param layout Grid layout of the current patch.
     */
    void classify(GridLayoutT const& layout)
    {
        auto classifier = classifier_type::withDefaults(*geometry_, layout);
        classifier(layout, meshData_);
    }

    // -------------------------------------------------------------------------
    //  BC application
    // -------------------------------------------------------------------------

    /**
     * @brief Apply the inner BC for a scalar quantity.
     *
     * @param field  Field whose ghost cells are updated in place.
     * @param layout Grid layout of the current patch.
     * @param ctx    State context (statenew, state, time, dt).
     */
    void applyBC(FieldT& field, GridLayoutT const& layout, context_type const& ctx)
    {
        auto it = scalarBCs_.find(field.physicalQuantity());
        if (it == scalarBCs_.end())
            return; // quantity not registered — no-op
        it->second->apply(field, layout, meshData_, ctx);
    }

    /**
     * @brief Apply the inner BC for a vector quantity.
     *
     * @param vecfield VecField whose ghost cells are updated in place.
     * @param layout   Grid layout of the current patch.
     * @param ctx      State context (statenew, state, time, dt).
     */
    void applyBC(vecfield_type& vecfield, GridLayoutT const& layout, context_type const& ctx)
    {
        auto it = vectorBCs_.find(vecfield.physicalQuantity());
        if (it == vectorBCs_.end())
            return; // quantity not registered — no-op
        it->second->apply(vecfield, layout, meshData_, ctx);
    }

    /**
     * @brief Apply the inner BCs of the moment / energy quantities in priority order.
     *
     * Replaces the previously hardcoded rhoV → rho → Etot sequence duplicated at the call sites.
     * The order is defined in one place: each registered BC reports a priority() (default 0; the
     * total-energy-from-pressure BC reports a higher value so it runs last, after the ghost rho /
     * rhoV it depends on are filled). B is excluded — its inner BC is a no-op (B is enforced via
     * E + constrained transport, never written here) — and E is applied separately in ComputeFluxes.
     *
     * Fields are taken from the (mutable) ctx.statenew, which is the state being updated.
     */
    void applyToMoments(GridLayoutT const& layout, context_type const& ctx)
    {
        auto& state = ctx.statenew;
        std::array<std::pair<int, std::function<void()>>, 3> tasks{{
            scalarTask_(state.rho, layout, ctx),
            vectorTask_(state.rhoV, layout, ctx),
            scalarTask_(state.Etot, layout, ctx),
        }};

        std::stable_sort(tasks.begin(), tasks.end(),
                         [](auto const& a, auto const& b) { return a.first < b.first; });

        for (auto const& [prio, run] : tasks)
            if (run)
                run();
    }

    // -------------------------------------------------------------------------
    //  Safe-state enforcement
    // -------------------------------------------------------------------------

    /**
     * @brief Pin a scalar field to its prescribed safe value on elements tagged @p status.
     *
     * The safe value comes from the per-quantity map built at construction (rho, P from the
     * config; Etot derived). Throws if the field's quantity has no registered safe value.
     */
    void setSafeState(FieldT& field, GridLayoutT const& layout,
                      ElemStatus status = ElemStatus::Inactive)
    {
        auto it = safeScalars_.find(field.physicalQuantity());
        if (it == safeScalars_.end())
            throw std::runtime_error(
                "InnerBoundaryManager::setSafeState: no safe value for this scalar quantity");

        double const value = it->second;
        auto const centering = layout.centering(field.physicalQuantity());
        auto& statusField    = meshData_.getStatusFieldFromCentering(centering);
        layout.evalOnGhostBox(field, [&](auto&... args) {
            auto const idx = MeshIndex<dimension>{args...};
            if (statusField(idx) == toDouble(status))
                field(idx) = value;
        });
    }

    /**
     * @brief Pin a vector field to its prescribed safe value on elements tagged @p status.
     *
     * The safe value comes from the per-quantity map (velocity, B from the config; rhoV
     * derived). Throws if the vecfield's quantity has no registered safe value.
     */
    void setSafeState(vecfield_type& vecfield, GridLayoutT const& layout,
                      ElemStatus status = ElemStatus::Inactive)
    {
        auto it = safeVectors_.find(vecfield.physicalQuantity());
        if (it == safeVectors_.end())
            throw std::runtime_error(
                "InnerBoundaryManager::setSafeState: no safe value for this vector quantity");

        auto const& values = it->second;
        auto comps         = vecfield.components();
        for_N<3>([&](auto ci) {
            constexpr auto c     = ci();
            auto& comp           = std::get<c>(comps);
            auto const centering = layout.centering(comp.physicalQuantity());
            auto& statusField    = meshData_.getStatusFieldFromCentering(centering);
            layout.evalOnGhostBox(comp, [&](auto&... args) {
                auto const idx = MeshIndex<dimension>{args...};
                if (statusField(idx) == toDouble(status))
                    comp(idx) = values[c];
            });
        });
    }

    // -------------------------------------------------------------------------
    //  Mesh data accessor
    // -------------------------------------------------------------------------

    NO_DISCARD mesh_data_type& getMeshData() { return meshData_; }

    NO_DISCARD mesh_data_type const& getMeshData() const { return meshData_; }

    NO_DISCARD geometry_type const& getGeometry() const { return *geometry_; }

    // -------------------------------------------------------------------------
    //  ResourcesUser interface (exposes mesh data to SAMRAI resource manager)
    // -------------------------------------------------------------------------

    NO_DISCARD bool isUsable() const { return meshData_.isUsable(); }

    NO_DISCARD bool isSettable() const { return meshData_.isSettable(); }

    NO_DISCARD auto getCompileTimeResourcesViewList() const
    {
        return std::forward_as_tuple(meshData_);
    }

    NO_DISCARD auto getCompileTimeResourcesViewList() { return std::forward_as_tuple(meshData_); }

private:
    /// Build a (priority, applier) task for a scalar field, or an empty applier if unregistered.
    std::pair<int, std::function<void()>>
    scalarTask_(FieldT& field, GridLayoutT const& layout, context_type const& ctx)
    {
        auto it = scalarBCs_.find(field.physicalQuantity());
        if (it == scalarBCs_.end())
            return {0, {}};
        auto* bc = it->second.get();
        return {bc->priority(),
                [this, bc, &field, &layout, &ctx]() { bc->apply(field, layout, meshData_, ctx); }};
    }

    /// Build a (priority, applier) task for a vector field, or an empty applier if unregistered.
    std::pair<int, std::function<void()>>
    vectorTask_(vecfield_type& vecfield, GridLayoutT const& layout, context_type const& ctx)
    {
        auto it = vectorBCs_.find(vecfield.physicalQuantity());
        if (it == vectorBCs_.end())
            return {0, {}};
        auto* bc = it->second.get();
        return {bc->priority(), [this, bc, &vecfield, &layout, &ctx]() {
                    bc->apply(vecfield, layout, meshData_, ctx);
                }};
    }

    /// Build the per-quantity safe-value maps from the prescribed primitives, deriving the
    /// conserved momentum and total energy. Requires thermo_ to be set.
    void buildSafeState_(SafeStateValues const& s)
    {
        using Scalar = scalar_qty;
        using Vector = vector_qty;

        safeScalars_[Scalar::rho] = s.density;
        safeScalars_[Scalar::P]   = s.pressure;

        // total energy from the safe primitives.
        // Needs the equation of state — skipped if no thermo (e.g. unit tests that never
        // call setSafeState(Etot)); setSafeState would then throw for Etot.
        if (thermo_)
        {
            thermo_->setState_DP(s.density, s.pressure);
            auto const e_int            = s.density * thermo_->internalEnergy();
            safeScalars_[Scalar::Etot] = totalEnergyFromInternalEnergy(
                e_int, s.density, s.velocity[0], s.velocity[1], s.velocity[2], s.B[0], s.B[1],
                s.B[2]);
        }

        safeVectors_[Vector::V]    = s.velocity;
        safeVectors_[Vector::rhoV] = {s.density * s.velocity[0], s.density * s.velocity[1],
                                      s.density * s.velocity[2]};
        safeVectors_[Vector::B]    = s.B;
    }

    std::unique_ptr<geometry_type> geometry_;
    mesh_data_type meshData_;
    std::shared_ptr<Thermo> thermo_;
    scalar_bc_map_type scalarBCs_;
    vector_bc_map_type vectorBCs_;
    std::unordered_map<scalar_qty, double> safeScalars_;
    std::unordered_map<vector_qty, std::array<double, 3>> safeVectors_;
};

} // namespace PHARE::core

#endif // PHARE_CORE_INNER_BOUNDARY_INNER_BOUNDARY_MANAGER_HPP
