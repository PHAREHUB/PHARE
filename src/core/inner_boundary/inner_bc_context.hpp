#ifndef PHARE_CORE_INNER_BOUNDARY_INNER_BC_CONTEXT_HPP
#define PHARE_CORE_INNER_BOUNDARY_INNER_BC_CONTEXT_HPP

namespace PHARE::core
{

/**
 * @brief Bundle of state information passed to inner boundary condition appliers.
 *
 * Naming mirrors EulerUsingComputedFlux: @p statenew is the state being actively
 * updated (primary), @p state is the previous state available for future
 * complex boundary conditions (e.g. Robin-type or zone-dependent BCs).
 *
 * When there is no distinct previous state (e.g. in ComputeFluxes), the caller
 * simply passes the same reference for both @p statenew and @p state.
 *
 * @note @p statenew is a NON-const reference: most BCs only read it, but some
 *       (e.g. FieldTotalEnergyFromPressureInnerBoundaryCondition) need to recompute
 *       a derived field such as pressure in place before reconstructing their own
 *       ghost values. Accessing a reference member through a const InnerBCContext
 *       still yields the non-const referent, so this stays mutable inside apply().
 *
 * @tparam PhysicalStateT  Physical state type (MHDState, HybridState, …).
 */
template<typename PhysicalStateT>
struct InnerBCContext
{
    PhysicalStateT&       statenew; ///< Updated state (primary reference, mutable for derived fields).
    PhysicalStateT const& state;    ///< Previous state, may alias statenew when not applicable.
    double                time{0.}; ///< Current simulation time.
    double                dt{0.};   ///< Time-step size (zero until needed).
};

} // namespace PHARE::core

#endif // PHARE_CORE_INNER_BOUNDARY_INNER_BC_CONTEXT_HPP
