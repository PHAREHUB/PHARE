#ifndef PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_FIELD_ADAPTIVE_OUTFLOW_PRESSURE_BOUNDARY_CONDITION_HPP
#define PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_FIELD_ADAPTIVE_OUTFLOW_PRESSURE_BOUNDARY_CONDITION_HPP

#include "core/boundary/boundary_defs.hpp"
#include "core/data/grid/gridlayout.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/data/vecfield/vecfield.hpp"
#include "core/numerics/boundary_condition/field_boundary_condition.hpp"
#include "core/numerics/primite_conservative_converter/conversion_utils.hpp"
#include "core/numerics/riemann_solvers/mhd_speeds.hpp"
#include "core/numerics/thermo/thermo.hpp"
#include "core/utilities/point/point.hpp"

#include <array>
#include <map>
#include <memory>

namespace PHARE::core
{

/**
 * @brief Pressure boundary condition for an MHD outflow that adapts, per tangential
 * ghost-column, between a zero-gradient (Neumann) and a fixed-target (Dirichlet) pressure.
 *
 * Intended to be used as the @c P_bc of a @c FieldTotalEnergyFromPressureBoundaryCondition:
 * the energy BC reconstructs Etot from whatever ghost pressure this condition writes, so the
 * whole super-magnetofast ↔ fixed-pressure switch reduces to choosing the pressure ghost.
 *
 * For each ghost column the outward normal velocity @c u_n and the boundary-normal fast
 * magnetosonic speed @c c_f are evaluated at the first interior cell (previous-substage state):
 *   * @c u_n >= c_f  → super-magnetofast: zero-gradient pressure (Neumann). With ρ, ρv, B also
 *                      zero-gradient, the reconstructed energy equals the interior energy — i.e.
 *                      the @c open / super-magnetofast-outflow behavior.
 *   * @c u_n <  c_f  → sub-fast: the pressure is pinned to @c p_target (Dirichlet) — the
 *                      @c fixed-pressure-outflow behavior.
 *
 * @tparam FieldT       Scalar field type (must satisfy IsField).
 * @tparam GridLayoutT  Grid layout type (must satisfy IsGridLayout).
 */
template<typename FieldT, typename GridLayoutT>
class FieldAdaptiveOutflowPressureBoundaryCondition
    : public IFieldBoundaryCondition<FieldT, GridLayoutT>
{
public:
    using Super                  = IFieldBoundaryCondition<FieldT, GridLayoutT>;
    using field_type             = typename Super::field_type;
    using physical_quantity_type = typename GridLayoutT::Quantity;
    using scalar_quantity_type   = typename physical_quantity_type::Scalar;
    using vector_quantity_type   = typename physical_quantity_type::Vector;

    static constexpr std::size_t dimension = Super::dimension;
    static constexpr std::size_t N         = Super::N;
    static_assert(N == 1,
                  "FieldAdaptiveOutflowPressureBoundaryCondition only applies to scalar fields.");

    FieldAdaptiveOutflowPressureBoundaryCondition(double p_target, std::shared_ptr<Thermo> thermo)
        : p_target_{p_target}
        , thermo_{std::move(thermo)}
    {
    }

    virtual ~FieldAdaptiveOutflowPressureBoundaryCondition() = default;

    FieldBoundaryConditionType getType() const override
    {
        return FieldBoundaryConditionType::AdaptiveOutflowPressure;
    }

    void apply(FieldT& PField, BoundaryLocation const boundaryLocation,
               Box<std::uint32_t, dimension> const& localGhostBox, GridLayoutT const& gridLayout,
               typename Super::boundary_condition_context_type const& ctx) override
    {
        Direction const direction    = getDirection(boundaryLocation);
        Side const side              = getSide(boundaryLocation);
        std::size_t const dir_n      = static_cast<std::size_t>(direction);
        double const sign_n          = (side == Side::Upper) ? +1.0 : -1.0;
        QtyCentering const centering = GridLayoutT::centering(PField.physicalQuantity())[dir_n];

        std::uint32_t const interior_n
            = (side == Side::Upper) ? gridLayout.physicalEndIndex(QtyCentering::dual, direction)
                                    : gridLayout.physicalStartIndex(QtyCentering::dual, direction);

        // B is read from the current state (CT-updated); projected to cell center for c_f.
        auto BField = ctx.accessor_new.getVecField(vector_quantity_type::B);
        auto BComps = BField.components();
        auto& Bx    = std::get<0>(BComps);
        auto& By    = std::get<1>(BComps);
        auto& Bz    = std::get<2>(BComps);

        auto projectB = [&](auto const& idx) {
            double const bx
                = GridLayoutT::template project<GridLayoutT::faceXToCellCenter>(Bx, idx);
            double const by
                = GridLayoutT::template project<GridLayoutT::faceYToCellCenter>(By, idx);
            double const bz
                = GridLayoutT::template project<GridLayoutT::faceZToCellCenter>(Bz, idx);
            return std::array<double, 3>{bx, by, bz};
        };

        auto const pFieldBox = gridLayout.toFieldBox(localGhostBox, PField.physicalQuantity());

        std::map<TangentialKey, bool> superfastCache;

        for (auto const& index : pFieldBox)
        {
            bool const superfast
                = regimeForSlice_(superfastCache, index, dir_n, sign_n, interior_n, projectB, ctx);
            auto const mirrorIdx = gridLayout.boundaryMirrored(direction, side, centering, index);
            if (superfast)
                PField(index) = PField(mirrorIdx); // Neumann (zero-gradient)
            else
                PField(index) = 2.0 * p_target_ - PField(mirrorIdx); // Dirichlet(p_target)
        }
    }

private:
    using GhostIdx      = Point<std::uint32_t, dimension>;
    using TangentialKey = std::array<std::uint32_t, (dimension == 0 ? 1 : dimension)>;

    /// Pack the tangential coordinates of @p idx (excluding the normal axis) into a map key —
    /// groups the ghost cells of one column so the regime is decided once per column.
    static TangentialKey tangKey_(GhostIdx const& idx, std::size_t dir_n)
    {
        TangentialKey k{};
        std::size_t dir_t = 0;
        for (std::size_t i = 0; i < dimension; ++i)
            if (i != dir_n)
            {
                k[dir_t] = idx[i];
                ++dir_t;
            }
        return k;
    }

    /// Decide, for the column through @p ghostIdx, whether the first interior cell is
    /// super-magnetofast in the outward-normal direction (u_n >= c_f). Reads ρ, ρv and Etot
    /// from the previous-substage state (interior is valid and substage-stable); B from the
    /// current state. Cached per tangential column.
    template<typename ProjectB>
    bool regimeForSlice_(std::map<TangentialKey, bool>& cache, GhostIdx const& ghostIdx,
                         std::size_t dir_n, double sign_n, std::uint32_t interior_n,
                         ProjectB const& projectB,
                         typename Super::boundary_condition_context_type const& ctx)
    {
        auto const key = tangKey_(ghostIdx, dir_n);
        if (auto it = cache.find(key); it != cache.end())
            return it->second;

        GhostIdx interiorIdx = ghostIdx;
        interiorIdx[dir_n]   = interior_n;

        auto& rho_old  = ctx.accessor_old.getField(scalar_quantity_type::rho);
        auto& Etot_old = ctx.accessor_old.getField(scalar_quantity_type::Etot);
        auto rhoV_old  = ctx.accessor_old.getVecField(vector_quantity_type::rhoV);
        auto rhoVc     = rhoV_old.components();

        double const rho_i = rho_old(interiorIdx);
        double const vx    = std::get<0>(rhoVc)(interiorIdx) / rho_i;
        double const vy    = std::get<1>(rhoVc)(interiorIdx) / rho_i;
        double const vz    = std::get<2>(rhoVc)(interiorIdx) / rho_i;
        auto const Bv      = projectB(interiorIdx);

        double const e_int = internalEnergyFromTotalEnergy(Etot_old(interiorIdx), rho_i, vx, vy, vz,
                                                           Bv[0], Bv[1], Bv[2]);
        thermo_->setState_DU(rho_i, e_int / rho_i);
        double const cs = thermo_->soundSpeed();

        std::array<double, 3> const V{vx, vy, vz};
        double const u_n   = sign_n * V[dir_n];
        double const Bn    = Bv[dir_n];
        double const BdotB = Bv[0] * Bv[0] + Bv[1] * Bv[1] + Bv[2] * Bv[2];
        double const c_f   = compute_fast_magnetosonic_from_cs_(cs, rho_i, Bn, BdotB);

        bool const superfast = (u_n >= c_f);
        return cache.emplace(key, superfast).first->second;
    }

    double p_target_;
    std::shared_ptr<Thermo> thermo_;
};

} // namespace PHARE::core
#endif // PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_FIELD_ADAPTIVE_OUTFLOW_PRESSURE_BOUNDARY_CONDITION_HPP
