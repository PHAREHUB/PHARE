#include "gtest/gtest.h"

#include <cstdint>
#include <memory>
#include <vector>

#include "core/data/field/field.hpp"
#include "core/data/grid/grid.hpp"
#include "core/data/grid/gridlayout.hpp"
#include "core/data/grid/gridlayoutimplyee_mhd.hpp"
#include "core/data/ndarray/ndarray_vector.hpp"
#include "core/inner_boundary/field_total_energy_from_pressure_inner_boundary_condition.hpp"
#include "core/inner_boundary/field_neumann_inner_boundary_condition.hpp"
#include "core/inner_boundary/field_dirichlet_inner_boundary_condition.hpp"
#include "core/inner_boundary/inner_boundary_mesh_classifier.hpp"
#include "core/inner_boundary/plane_inner_boundary.hpp"
#include "core/mhd/mhd_quantities.hpp"
#include "core/numerics/thermo/ideal_gas_thermo.hpp"
#include "core/utilities/box/box.hpp"
#include "tests/core/data/vecfield/test_vecfield_fixtures_mhd.hpp"

namespace
{
using namespace PHARE::core;

constexpr double eps = 1e-12;

using GridLayoutImpl = GridLayoutImplYeeMHD<2, 2, 1>;
using GridLayout     = PHARE::core::GridLayout<GridLayoutImpl>;
using Classifier     = InnerBoundaryMeshClassifier<2, GridLayout, MHDQuantity>;
using MeshData       = InnerBoundaryMeshData<2, MHDQuantity>;
using ScalarField    = Field<2, MHDQuantity::Scalar, double>;

constexpr std::array<QtyCentering, 2> kCellC = {QtyCentering::dual, QtyCentering::dual};

constexpr double gamma = 5.0 / 3.0;

/// MHD-like state holding the conservative fields the energy BC needs. Scalars are Field views
/// over owned NdArrayVector storage (matching the real MHDState, whose members are field views);
/// vectors are self-contained UsableVecFieldMHD (a TensorField, so .components() works).
struct MhdState
{
    NdArrayVector<2, double> rhoStore;
    NdArrayVector<2, double> PStore;
    NdArrayVector<2, double> EtotStore;

    ScalarField rho;
    ScalarField P;
    ScalarField Etot;

    UsableVecFieldMHD<2> rhoV;
    UsableVecFieldMHD<2> B;

    explicit MhdState(GridLayout const& layout)
        : rhoStore{layout.allocSize(MHDQuantity::Scalar::rho)}
        , PStore{layout.allocSize(MHDQuantity::Scalar::P)}
        , EtotStore{layout.allocSize(MHDQuantity::Scalar::Etot)}
        , rho{"rho", MHDQuantity::Scalar::rho, rhoStore.data(), rhoStore.shape()}
        , P{"P", MHDQuantity::Scalar::P, PStore.data(), PStore.shape()}
        , Etot{"Etot", MHDQuantity::Scalar::Etot, EtotStore.data(), EtotStore.shape()}
        , rhoV{"rhoV", layout, MHDQuantity::Vector::rhoV}
        , B{"B", layout, MHDQuantity::Vector::B}
    {
    }
};

/// Allocate and wire all mesh-data buffers (same pattern as the per-BC inner tests).
struct MeshDataBuffers
{
    static constexpr char const* BOUNDARY_NAME = "test";

    explicit MeshDataBuffers(GridLayout const& layout)
        : phi_storage{layout.allocSize(MHDQuantity::Scalar::NodeCentered)}
        , tags{BOUNDARY_NAME}
    {
        ScalarField phi_field{std::string(BOUNDARY_NAME) + "_signed_distance",
                              MHDQuantity::Scalar::NodeCentered, phi_storage.data(),
                              phi_storage.shape()};
        tags.signedDistanceAtNodes.setBuffer(&phi_field);

        elem_storages.reserve(MeshData::num_elem_types);
        for (std::size_t i = 0; i < MeshData::num_elem_types; ++i)
        {
            auto const c   = MeshData::idxToCentering(i);
            auto const qty = MeshData::scalarFromCentering(c);
            elem_storages.emplace_back(layout.allocSize(qty));
            ScalarField tmp{tags.elemStatus[i].name(), qty, elem_storages[i].data(),
                            elem_storages[i].shape()};
            tags.elemStatus[i].setBuffer(&tmp);
        }
        tags.ghostElemsData._data    = &ghost_array;
        tags.degradedElemsData._data = &degraded_array;
    }

    NdArrayVector<2, double> phi_storage;
    std::vector<NdArrayVector<2, double>> elem_storages;
    GhostElemPack<2>::ghost_elem_array_type ghost_array{};
    GhostElemPack<2>::ghost_elem_array_type degraded_array{};
    MeshData tags;
};

/// Plane at x = 0, normal (1,0). Same 4×2 geometry as the other inner-BC tests.
struct Fixture
{
    PlaneInnerBoundary<2> plane{"plane", {0.0, 0.0}, {1.0, 0.0}};
    Box<int, 2> amr_box{{-2, 0}, {1, 1}};
    GridLayout layout{{1.0, 1.0}, {4u, 2u}, {0.0, 0.0}, amr_box};

    MeshDataBuffers buffers{layout};
    MhdState state{layout};

    Fixture()
    {
        Classifier::Overrides ov;
        ov.cut_eps      = 1e-12;
        ov.inactive_eps = 0.0;
        auto classifier = Classifier::withDefaults(plane, layout, ov);
        classifier(layout, buffers.tags);
    }

    template<typename F>
    static void fillAll(F& field, double value)
    {
        for (auto i = 0u; i < field.shape()[0]; ++i)
            for (auto j = 0u; j < field.shape()[1]; ++j)
                field(i, j) = value;
    }

    auto makeBC()
    {
        auto pressureBC = std::make_unique<
            FieldNeumannInnerBoundaryCondition<ScalarField, GridLayout, MhdState>>();
        auto thermo = std::make_shared<IdealGasThermo>(gamma);
        return FieldTotalEnergyFromPressureInnerBoundaryCondition<ScalarField, GridLayout,
                                                                  MhdState>{std::move(pressureBC),
                                                                            thermo};
    }

    // Variant whose pressure sub-BC is a constant (0th-order) Dirichlet at p_bc: it fills every
    // ghost P, so the energy reconstruction also covers non-interpolable ghosts.
    auto makeConstantDirichletBC(double p_bc)
    {
        using Dirichlet = FieldDirichletInnerBoundaryCondition<ScalarField, GridLayout, MhdState>;
        auto pressureBC = std::make_unique<Dirichlet>(p_bc, Dirichlet::ExtrapolationType::Constant);
        auto thermo     = std::make_shared<IdealGasThermo>(gamma);
        return FieldTotalEnergyFromPressureInnerBoundaryCondition<ScalarField, GridLayout,
                                                                  MhdState>{std::move(pressureBC),
                                                                            thermo};
    }
};

double totalEnergy(double rho, double vx, double vy, double vz, double bx, double by, double bz,
                   double P)
{
    return P / (gamma - 1.0) + 0.5 * rho * (vx * vx + vy * vy + vz * vz)
           + 0.5 * (bx * bx + by * by + bz * bz);
}

} // namespace

// ---------------------------------------------------------------------------
// A spatially uniform state is a fixed point: ghost Etot == interior Etot.
// ---------------------------------------------------------------------------
TEST(FieldTotalEnergyFromPressureInner, uniformStateGhostEtotEqualsInterior)
{
    Fixture fix;
    auto& st     = fix.state;
    auto& layout = fix.layout;

    constexpr double rho = 2.0, vx = 0.3, vy = -0.4, vz = 0.1;
    constexpr double bx = 0.5, by = -0.2, bz = 0.7, P = 1.5;
    double const etot = totalEnergy(rho, vx, vy, vz, bx, by, bz, P);

    Fixture::fillAll(st.rho, rho);
    Fixture::fillAll(st.P, P);
    Fixture::fillAll(st.Etot, etot);
    Fixture::fillAll(st.rhoV[0], rho * vx);
    Fixture::fillAll(st.rhoV[1], rho * vy);
    Fixture::fillAll(st.rhoV[2], rho * vz);
    Fixture::fillAll(st.B[0], bx);
    Fixture::fillAll(st.B[1], by);
    Fixture::fillAll(st.B[2], bz);

    auto bc = fix.makeBC();
    bc.apply(st.Etot, layout, fix.buffers.tags, InnerBCContext<MhdState>{st, st, 0.0});

    bool foundInPatch = false;
    for (auto const& g : fix.buffers.tags.getGhostDataFromCentering(kCellC))
    {
        if (!g.mirrorIsInterpolable)
            continue;
        foundInPatch = true;
        EXPECT_NEAR(st.Etot(g.index), etot, 1e-10)
            << "ghost (" << g.index[0] << "," << g.index[1] << ")";
    }
    EXPECT_TRUE(foundInPatch) << "at least one in-patch ghost must exist";
}

// ---------------------------------------------------------------------------
// Ghost Etot is reconstructed from the (Neumann) ghost pressure and the ghost
// momentum/density already present — NOT by mirroring the conserved energy.
//
// Interior is uniform; ghost cells carry a DIFFERENT rho/rhoV (as if their own
// BCs had already run). With zero-gradient pressure, the ghost energy must use
// the ghost kinetic term, so it differs from the interior energy.
// ---------------------------------------------------------------------------
TEST(FieldTotalEnergyFromPressureInner, ghostEnergyUsesGhostMomentumAndNeumannPressure)
{
    Fixture fix;
    auto& st     = fix.state;
    auto& layout = fix.layout;

    // interior (fluid) primitive state
    constexpr double rho_i = 2.0, vx_i = 0.5, vy_i = 0.0, vz_i = 0.0;
    constexpr double bx = 0.4, by = 0.0, bz = 0.0, P0 = 1.0;
    double const etot_i = totalEnergy(rho_i, vx_i, vy_i, vz_i, bx, by, bz, P0);

    // distinct ghost density / velocity (mimics rho & rhoV inner BCs having run)
    constexpr double rho_g = 3.0, vx_g = -0.2, vy_g = 0.6, vz_g = 0.0;

    Fixture::fillAll(st.rho, rho_i);
    Fixture::fillAll(st.P, P0);
    Fixture::fillAll(st.Etot, etot_i);
    Fixture::fillAll(st.rhoV[0], rho_i * vx_i);
    Fixture::fillAll(st.rhoV[1], rho_i * vy_i);
    Fixture::fillAll(st.rhoV[2], rho_i * vz_i);
    Fixture::fillAll(st.B[0], bx);
    Fixture::fillAll(st.B[1], by);
    Fixture::fillAll(st.B[2], bz);

    // overwrite the ghost cells' density / momentum with the distinct ghost values
    for (auto const& g : fix.buffers.tags.getGhostDataFromCentering(kCellC))
    {
        if (!g.mirrorIsInterpolable)
            continue;
        st.rho(g.index)     = rho_g;
        st.rhoV[0](g.index) = rho_g * vx_g;
        st.rhoV[1](g.index) = rho_g * vy_g;
        st.rhoV[2](g.index) = rho_g * vz_g;
    }

    double const etot_expected = totalEnergy(rho_g, vx_g, vy_g, vz_g, bx, by, bz, P0);

    auto bc = fix.makeBC();
    bc.apply(st.Etot, layout, fix.buffers.tags, InnerBCContext<MhdState>{st, st, 0.0});

    bool foundInPatch = false;
    for (auto const& g : fix.buffers.tags.getGhostDataFromCentering(kCellC))
    {
        if (!g.mirrorIsInterpolable)
            continue;
        foundInPatch = true;
        EXPECT_NEAR(st.Etot(g.index), etot_expected, 1e-10)
            << "ghost (" << g.index[0] << "," << g.index[1] << ")";
        // and it must NOT be the mirrored interior energy
        EXPECT_GT(std::abs(st.Etot(g.index) - etot_i), 1e-3)
            << "ghost energy should differ from mirrored interior energy";
    }
    EXPECT_TRUE(foundInPatch) << "at least one in-patch ghost must exist";
}

// ---------------------------------------------------------------------------
// Non-interpolable ghosts are left untouched.
// ---------------------------------------------------------------------------
TEST(FieldTotalEnergyFromPressureInner, nonInterpolableGhostsUntouched)
{
    Fixture fix;
    auto& st     = fix.state;
    auto& layout = fix.layout;

    constexpr double rho = 1.0, P = 1.0, bx = 0.0, by = 0.0, bz = 0.0;
    constexpr double sentinel = -999.0;
    double const etot         = totalEnergy(rho, 0, 0, 0, bx, by, bz, P);

    Fixture::fillAll(st.rho, rho);
    Fixture::fillAll(st.P, P);
    Fixture::fillAll(st.Etot, etot);
    Fixture::fillAll(st.rhoV[0], 0.0);
    Fixture::fillAll(st.rhoV[1], 0.0);
    Fixture::fillAll(st.rhoV[2], 0.0);
    Fixture::fillAll(st.B[0], bx);
    Fixture::fillAll(st.B[1], by);
    Fixture::fillAll(st.B[2], bz);

    // mark every non-interpolable ghost with a sentinel
    for (auto const& g : fix.buffers.tags.getGhostDataFromCentering(kCellC))
        if (!g.mirrorIsInterpolable)
            st.Etot(g.index) = sentinel;

    auto bc = fix.makeBC();
    bc.apply(st.Etot, layout, fix.buffers.tags, InnerBCContext<MhdState>{st, st, 0.0});

    for (auto const& g : fix.buffers.tags.getGhostDataFromCentering(kCellC))
    {
        if (!g.mirrorIsInterpolable)
        {
            EXPECT_NEAR(st.Etot(g.index), sentinel, eps)
                << "non-interpolable ghost (" << g.index[0] << "," << g.index[1] << ") changed";
        }
    }
}

// ---------------------------------------------------------------------------
// With a CONSTANT Dirichlet pressure sub-BC, ghost P is filled on every ghost
// (interpolable or not), so Etot is reconstructed even on non-interpolable
// ghosts — using the ghost rho/rhoV/B and the prescribed pressure.
// ---------------------------------------------------------------------------
TEST(FieldTotalEnergyFromPressureInner, constantPressureSubBCFillsNonInterpolableGhosts)
{
    Fixture fix;
    auto& st     = fix.state;
    auto& layout = fix.layout;

    constexpr double rho_i = 2.0, P_bc = 1.7;
    constexpr double bx = 0.3, by = -0.1, bz = 0.2;
    // distinct ghost state, as if rho/rhoV inner BCs had already filled every ghost
    constexpr double rho_g = 3.0, vx_g = -0.2, vy_g = 0.6, vz_g = 0.1;
    constexpr double sentinel = -999.0;

    Fixture::fillAll(st.rho, rho_i);
    Fixture::fillAll(st.P, P_bc);
    Fixture::fillAll(st.Etot, totalEnergy(rho_i, 0, 0, 0, bx, by, bz, P_bc));
    Fixture::fillAll(st.rhoV[0], 0.0);
    Fixture::fillAll(st.rhoV[1], 0.0);
    Fixture::fillAll(st.rhoV[2], 0.0);
    Fixture::fillAll(st.B[0], bx);
    Fixture::fillAll(st.B[1], by);
    Fixture::fillAll(st.B[2], bz);

    // every ghost carries the distinct ghost rho/rhoV; non-interpolable ones get a sentinel Etot
    bool sawNonInterpolable = false;
    for (auto const& g : fix.buffers.tags.getGhostDataFromCentering(kCellC))
    {
        st.rho(g.index)     = rho_g;
        st.rhoV[0](g.index) = rho_g * vx_g;
        st.rhoV[1](g.index) = rho_g * vy_g;
        st.rhoV[2](g.index) = rho_g * vz_g;
        if (!g.mirrorIsInterpolable)
        {
            sawNonInterpolable = true;
            st.Etot(g.index)  = sentinel;
        }
    }
    ASSERT_TRUE(sawNonInterpolable)
        << "fixture should expose at least one non-interpolable ghost to exercise the fix";

    double const etot_expected = totalEnergy(rho_g, vx_g, vy_g, vz_g, bx, by, bz, P_bc);

    auto bc = fix.makeConstantDirichletBC(P_bc);
    bc.apply(st.Etot, layout, fix.buffers.tags, InnerBCContext<MhdState>{st, st, 0.0});

    for (auto const& g : fix.buffers.tags.getGhostDataFromCentering(kCellC))
        if (!g.mirrorIsInterpolable)
        {
            EXPECT_NEAR(st.Etot(g.index), etot_expected, 1e-10)
                << "non-interpolable ghost (" << g.index[0] << "," << g.index[1]
                << ") not reconstructed with constant pressure sub-BC";
        }
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
