#ifndef PHARE_CORE_INNER_BOUNDARY_FIELD_TOTAL_ENERGY_FROM_PRESSURE_INNER_BOUNDARY_CONDITION_HPP
#define PHARE_CORE_INNER_BOUNDARY_FIELD_TOTAL_ENERGY_FROM_PRESSURE_INNER_BOUNDARY_CONDITION_HPP

#include "core/inner_boundary/field_inner_boundary_condition.hpp"
#include "core/inner_boundary/inner_boundary_defs.hpp"
#include "core/numerics/primite_conservative_converter/conversion_utils.hpp"
#include "core/numerics/thermo/thermo.hpp"
#include "core/utilities/index/index.hpp"

#include <cmath>
#include <limits>
#include <memory>
#include <tuple>

namespace PHARE::core
{
/**
 * @brief Inner-boundary condition for the total energy Etot derived from a zero-gradient
 *        condition on PRESSURE rather than on the conserved energy itself.
 *
 * Applying a plain Neumann condition directly on Etot mirrors the *conserved* energy, which
 * lumps internal + kinetic + magnetic energy together. When a strong gradient (shock / contact)
 * crosses the embedded boundary this produces spurious oscillations and over/under-heating at the
 * surface — the well-documented ghost-cell failure mode (Modified Ghost Fluid Method; inverse
 * Lax–Wendroff / hybrid-reconstruction literature).
 *
 * Instead, this BC imposes the zero-gradient condition on the *pressure* and reconstructs Etot in
 * the ghost cells from that pressure plus the already-filled ghost rho, rhoV, B. It mirrors the
 * outer-boundary FieldTotalEnergyFromPressureBoundaryCondition but uses the inner-boundary
 * machinery (continuous mirror points + linear point interpolator instead of integer mirror cells).
 *
 * apply():
 *   1. Reconstruct interior pressure on fluid/cut cells from the current conservatives, so the
 *      pressure sub-BC has valid interior P to interpolate at the continuous mirror points.
 *   2. Fill ghost pressure via the owned pressure sub-BC (a Neumann condition by default).
 *   3. Reconstruct ghost Etot from the freshly filled ghost P + ghost rho/rhoV/B.
 *
 * @note rho, rhoV, P and Etot are co-located (all-dual / cell-centered), so the cell-centered
 *       ghostElem.mirrorIsInterpolable flag already covers rho and rhoV — no extra interpolability
 *       check is needed for them. B is never written directly into ghost cells (its boundary
 *       condition is enforced via the electric field and constrained transport to keep div B = 0);
 *       the CT-maintained B is simply read here.
 *
 * @note Only meaningful for a scalar field (Etot). The vector instantiation is a no-op.
 */
template<typename ScalarOrTensorFieldT, typename GridLayoutT, typename PhysicalStateT>
class FieldTotalEnergyFromPressureInnerBoundaryCondition
    : public FieldInnerBoundaryCondition<ScalarOrTensorFieldT, GridLayoutT, PhysicalStateT>
{
public:
    using Super = FieldInnerBoundaryCondition<ScalarOrTensorFieldT, GridLayoutT, PhysicalStateT>;

    using field_type                    = Super::field_type;
    using inner_boundary_mesh_data_type = Super::inner_boundary_mesh_data_type;
    using ghost_elem_data_type          = Super::ghost_elem_data_type;
    using interpolator_type             = Super::interpolator_type;
    using context_type                  = Super::context_type;

    /// Base type of the owned scalar (pressure) sub-BC.
    using scalar_inner_bc_type
        = FieldInnerBoundaryCondition<field_type, GridLayoutT, PhysicalStateT>;

    static constexpr size_t dimension = Super::dimension;
    static constexpr bool is_scalar   = Super::is_scalar;

    /// Priority ensuring this BC runs after rho / rhoV ghosts are filled (see priority()).
    static constexpr int late_priority = 10;


    FieldTotalEnergyFromPressureInnerBoundaryCondition(
        std::unique_ptr<scalar_inner_bc_type> pressureBC, std::shared_ptr<Thermo> thermo)
        : pressureBC_{std::move(pressureBC)}
        , thermo_{std::move(thermo)}
    {
    }


    FieldInnerBoundaryConditionType getType() const override
    {
        return FieldInnerBoundaryConditionType::Neumann;
    }

    int priority() const override { return late_priority; }

    void apply(ScalarOrTensorFieldT& scalarOrTensorField, GridLayoutT const& layout,
               inner_boundary_mesh_data_type const& boundaryMeshData,
               context_type const& ctx) override
    {
        // Etot is a scalar quantity; the tensor instantiation is a meaningless no-op.
        if constexpr (!is_scalar)
        {
            return;
        }
        else
        {
            auto& EtotField = scalarOrTensorField;
            auto& P          = ctx.statenew.P;
            auto& rho        = ctx.statenew.rho;

            auto rhoVc  = ctx.statenew.rhoV.components();
            auto& rhoVx = std::get<0>(rhoVc);
            auto& rhoVy = std::get<1>(rhoVc);
            auto& rhoVz = std::get<2>(rhoVc);

            auto Bc  = ctx.statenew.B.components();
            auto& Bx = std::get<0>(Bc);
            auto& By = std::get<1>(Bc);
            auto& Bz = std::get<2>(Bc);

            auto const& cellStatus = boundaryMeshData.cellStatusField();

            // -- Step 1: reconstruct interior pressure on fluid/cut cells from conservatives, so
            //    the pressure sub-BC has valid interior P to interpolate at the mirror points.
            layout.evalOnGhostBox(P, [&](auto&... args) {
                auto const idx = core::MeshIndex<dimension>{args...};
                auto const st  = cellStatus(idx);
                if (st != toDouble(ElemStatus::Fluid) && st != toDouble(ElemStatus::Cut))
                    return;

                auto const rho_m = rho(idx);
                if (rho_m <= min_density_)
                    return;

                auto const vx = rhoVx(idx) / rho_m;
                auto const vy = rhoVy(idx) / rho_m;
                auto const vz = rhoVz(idx) / rho_m;

                auto const bx
                    = GridLayoutT::template project<GridLayoutT::faceXToCellCenter>(Bx, idx);
                auto const by
                    = GridLayoutT::template project<GridLayoutT::faceYToCellCenter>(By, idx);
                auto const bz
                    = GridLayoutT::template project<GridLayoutT::faceZToCellCenter>(Bz, idx);

                auto const e_int = internalEnergyFromTotalEnergy(EtotField(idx), rho_m, vx, vy, vz,
                                                                 bx, by, bz);
                thermo_->setState_DU(rho_m, e_int / rho_m);
                P(idx) = thermo_->pressure();
            });

            // -- Step 2: fill ghost pressure via the owned (Neumann) pressure sub-BC.
            pressureBC_->apply(P, layout, boundaryMeshData, ctx);

            // -- Step 3: reconstruct ghost Etot from the freshly filled ghost P + ghost
            //    rho/rhoV/B. rho/rhoV are co-located with Etot, so mirrorIsInterpolable covers
            //    them; B is read as maintained by constrained transport.
            auto const centering   = GridLayoutT::centering(EtotField);
            auto const& ghostElems = boundaryMeshData.getGhostDataFromCentering(centering);

            // When the pressure sub-BC fills every ghost (constant Dirichlet), step 2 set P on the
            // non-interpolable ghosts too; rho/rhoV are likewise filled there (co-located, applied
            // earlier by priority), so Etot can be reconstructed on those ghosts as well.
            bool const fillAll = pressureBC_->fillsNonInterpolableGhosts();

            for (ghost_elem_data_type const& ghostElem : ghostElems)
            {
                // Mirror not interpolable and a sub-BC that left this P ghost untouched: skip
                // (matches the "leave untouched" convention of every inner BC).
                if (!fillAll && !ghostElem.mirrorIsInterpolable)
                    continue;

                auto const idx = toMeshIndex_(ghostElem.index);

                auto const rho_g = rho(idx);
                if (rho_g <= min_density_)
                    continue;

                auto const vx = rhoVx(idx) / rho_g;
                auto const vy = rhoVy(idx) / rho_g;
                auto const vz = rhoVz(idx) / rho_g;

                auto const bx
                    = GridLayoutT::template project<GridLayoutT::faceXToCellCenter>(Bx, idx);
                auto const by
                    = GridLayoutT::template project<GridLayoutT::faceYToCellCenter>(By, idx);
                auto const bz
                    = GridLayoutT::template project<GridLayoutT::faceZToCellCenter>(Bz, idx);

                thermo_->setState_DP(rho_g, P(idx));
                auto const e_int = thermo_->internalEnergy() * rho_g;
                EtotField(idx)
                    = totalEnergyFromInternalEnergy(e_int, rho_g, vx, vy, vz, bx, by, bz);
            }
        }
    }

private:
    /// Convert a ghost-element index (Point<uint32_t, dim>) to a MeshIndex for project<>().
    static auto toMeshIndex_(Point<std::uint32_t, dimension> const& p)
    {
        if constexpr (dimension == 1)
            return core::MeshIndex<1>{p[0]};
        else if constexpr (dimension == 2)
            return core::MeshIndex<2>{p[0], p[1]};
        else
            return core::MeshIndex<3>{p[0], p[1], p[2]};
    }

    /// Density floor guarding the V = rhoV/rho division; mirrors min_value in to_primitive_converter.
    static double const min_density_;

    std::unique_ptr<scalar_inner_bc_type> pressureBC_;
    std::shared_ptr<Thermo>               thermo_;
};

template<typename ScalarOrTensorFieldT, typename GridLayoutT, typename PhysicalStateT>
double const FieldTotalEnergyFromPressureInnerBoundaryCondition<ScalarOrTensorFieldT, GridLayoutT,
                                                                PhysicalStateT>::min_density_
    = std::sqrt(1024 * std::numeric_limits<double>::min());

} // namespace PHARE::core
#endif // PHARE_CORE_INNER_BOUNDARY_FIELD_TOTAL_ENERGY_FROM_PRESSURE_INNER_BOUNDARY_CONDITION_HPP
