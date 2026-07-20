#include "gtest/gtest.h"

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include "core/data/field/field.hpp"
#include "core/data/grid/gridlayout.hpp"
#include "core/data/grid/gridlayoutimplyee_mhd.hpp"
#include "core/data/ndarray/ndarray_vector.hpp"
#include "core/data/vecfield/vecfield.hpp"
#include "core/inner_boundary/inner_boundary_manager.hpp"
#include "core/inner_boundary/inner_boundary_mesh_classifier.hpp"
#include "core/numerics/thermo/ideal_gas_thermo.hpp"
#include "core/inner_boundary/plane_inner_boundary.hpp"
#include "core/mhd/mhd_quantities.hpp"
#include "core/numerics/interpolator/field_at_point.hpp"
#include "core/utilities/box/box.hpp"
#include "initializer/data_provider.hpp"
#include "tests/core/data/vecfield/test_vecfield_fixtures_mhd.hpp"

namespace
{
constexpr double eps = 1e-10;

using GridLayoutImpl = PHARE::core::GridLayoutImplYeeMHD<2, 2, 1>;
using GridLayout     = PHARE::core::GridLayout<GridLayoutImpl>;
using ScalarField    = PHARE::core::Field<2, PHARE::core::MHDQuantity::Scalar, double>;
using VecFieldMHD2   = PHARE::core::UsableVecFieldMHD<2>;
using MeshData       = PHARE::core::InnerBoundaryMeshData<2, PHARE::core::MHDQuantity>;
using Classifier
    = PHARE::core::InnerBoundaryMeshClassifier<2, GridLayout, PHARE::core::MHDQuantity>;

/// Minimal physical-state stub. Holds MHD-like fields (with null buffers) so the
/// total-energy-from-pressure BC — instantiated by the Reflective factory — type-checks.
/// These are never dereferenced here: the BC is only exercised against a real state.
struct DummyState
{
    using VecF = PHARE::core::VecField<ScalarField, PHARE::core::MHDQuantity>;

    ScalarField rho{"rho", PHARE::core::MHDQuantity::Scalar::rho};
    ScalarField P{"P", PHARE::core::MHDQuantity::Scalar::P};
    ScalarField Etot{"Etot", PHARE::core::MHDQuantity::Scalar::Etot};
    VecF        rhoV{"rhoV", PHARE::core::MHDQuantity::Vector::rhoV};
    VecF        B{"B", PHARE::core::MHDQuantity::Vector::B};

    // mirror MHDState's accessor so criterion-based BCs (adaptive Dirichlet/Neumann) compile
    VecF& getVector(PHARE::core::MHDQuantity::Vector q)
    {
        if (q == PHARE::core::MHDQuantity::Vector::rhoV)
            return rhoV;
        if (q == PHARE::core::MHDQuantity::Vector::B)
            return B;
        throw std::runtime_error("DummyState::getVector: unsupported quantity");
    }
};

using Manager = PHARE::core::InnerBoundaryManager<PHARE::core::MHDQuantity, ScalarField, GridLayout,
                                                  DummyState>;

// ---------------------------------------------------------------------------
//  ManagerFixture — constructs a Manager directly (no dict needed)
//
//  Plane at x = 0, normal (1, 0).  Grid: 4 × 2 physical cells, dx=dy=1,
//  AMR box {{-2,0},{1,1}}.  Same geometry as all per-BC tests.
//
//  Storage arrays are owned by the fixture and wired into the manager's
//  MeshData via setBuffer(), exactly as MeshDataBuffers does for the
//  per-BC test fixtures.
// ---------------------------------------------------------------------------
struct ManagerFixture
{
    static constexpr char const* BOUNDARY_NAME = "test";

    PHARE::core::Box<int, 2> amr_box{{-2, 0}, {1, 1}};
    GridLayout layout{{1.0, 1.0}, {4u, 2u}, {0.0, 0.0}, amr_box};

    // Storage arrays — must be declared before manager so they outlive it.
    PHARE::core::NdArrayVector<2, double> phi_storage{
        layout.allocSize(PHARE::core::MHDQuantity::Scalar::NodeCentered)};
    std::vector<PHARE::core::NdArrayVector<2, double>> elem_storages;
    PHARE::core::GhostElemPack<2>::ghost_elem_array_type ghost_array{};
    PHARE::core::GhostElemPack<2>::ghost_elem_array_type degraded_array{};

    std::vector<PHARE::core::MHDQuantity::Scalar> scalarQtys{
        PHARE::core::MHDQuantity::Scalar::rho,
        PHARE::core::MHDQuantity::Scalar::Etot,
    };
    std::vector<PHARE::core::MHDQuantity::Vector> vectorQtys{
        PHARE::core::MHDQuantity::Vector::B,
        PHARE::core::MHDQuantity::Vector::rhoV,
        PHARE::core::MHDQuantity::Vector::E,
    };

    Manager manager;

    ManagerFixture(std::shared_ptr<PHARE::core::Thermo> thermo = nullptr,
                   Manager::SafeStateValues safe = {})
        : manager{std::make_unique<PHARE::core::PlaneInnerBoundary<2>>(
                      BOUNDARY_NAME, PHARE::core::Point<double, 2>{0.0, 0.0},
                      PHARE::core::Point<double, 2>{1.0, 0.0}),
                  PHARE::core::InnerBoundaryConditionType::Reflective, scalarQtys, vectorQtys,
                  std::move(thermo), 0., 0., safe}
    {
        // Wire storage into the manager's mesh data via temporary Field objects.
        // setBuffer() only retains the raw data pointer, so the temporaries can be
        // destroyed after the call; the NdArrayVector storage members keep the memory alive.
        auto& md = manager.getMeshData();

        std::string const bn{BOUNDARY_NAME};

        ScalarField phi_field{bn + "_signed_distance",
                              PHARE::core::MHDQuantity::Scalar::NodeCentered, phi_storage.data(),
                              phi_storage.shape()};
        md.signedDistanceAtNodes.setBuffer(&phi_field);

        elem_storages.reserve(MeshData::num_elem_types);
        for (std::size_t i = 0; i < MeshData::num_elem_types; ++i)
        {
            auto const c   = MeshData::idxToCentering(i);
            auto const qty = MeshData::scalarFromCentering(c);
            elem_storages.emplace_back(layout.allocSize(qty));
            ScalarField tmp{md.elemStatus[i].name(), qty, elem_storages[i].data(),
                            elem_storages[i].shape()};
            md.elemStatus[i].setBuffer(&tmp);
        }

        md.ghostElemsData._data    = &ghost_array;
        md.degradedElemsData._data = &degraded_array;

        // Classify mesh elements relative to the plane at x = 0.
        PHARE::core::PlaneInnerBoundary<2> plane{BOUNDARY_NAME, {0.0, 0.0}, {1.0, 0.0}};
        Classifier::Overrides ov;
        ov.cut_eps      = 1e-12;
        ov.inactive_eps = 0.0;
        auto classifier = Classifier::withDefaults(plane, layout, ov);
        classifier(layout, md);
    }
};

} // namespace


// ---------------------------------------------------------------------------
// Static factory: returns nullptr when no inner_boundary key in dict
// ---------------------------------------------------------------------------

TEST(InnerBoundaryManager, createReturnsNullptrWithoutInnerBoundaryKey)
{
    PHARE::initializer::PHAREDict dict;
    // dict has no "inner_boundary" key

    auto mgr = Manager::create(dict, {}, {}, nullptr);
    EXPECT_EQ(mgr, nullptr);
}


// ---------------------------------------------------------------------------
// Factory wiring: Reflective assigns the right BC type to each quantity
//
// We verify observable behaviour rather than implementation internals:
//   rho  (scalar) → Neumann   : ghost receives mirror value (constant field unchanged)
//   rhoV (vector) → Symmetric : normal component reversed, tangential preserved
//   E    (vector) → Antisymmetric : tangential component reversed, normal preserved
// ---------------------------------------------------------------------------

/**
 * @brief Scalar rho under Reflective gets Neumann: a constant field stays constant.
 */
TEST(InnerBoundaryManager, reflectiveScalarIsNeumann)
{
    ManagerFixture fix;
    auto const& layout = fix.layout;

    constexpr double C = 5.0;

    PHARE::core::NdArrayVector<2, double> storage{
        layout.allocSize(PHARE::core::MHDQuantity::Scalar::rho)};
    ScalarField rho{"rho", PHARE::core::MHDQuantity::Scalar::rho, storage.data(), storage.shape()};

    for (auto i = 0u; i < rho.shape()[0]; ++i)
        for (auto j = 0u; j < rho.shape()[1]; ++j)
            rho(i, j) = C;

    DummyState state;
    PHARE::core::InnerBCContext<DummyState> ctx{state, state, 0.0};
    fix.manager.applyBC(rho, layout, ctx);

    for (auto i = 0u; i < rho.shape()[0]; ++i)
        for (auto j = 0u; j < rho.shape()[1]; ++j)
            EXPECT_NEAR(rho(i, j), C, eps) << "cell (" << i << ", " << j << ") changed";
}

/**
 * @brief Vector rhoV under Reflective gets Symmetric: normal component is reversed.
 *
 * Constant field (Cx, Cy, Cz), plane normal n=(1,0,0).
 * Expected ghost: (-Cx, Cy, Cz).
 */
TEST(InnerBoundaryManager, reflectiveRhoVIsSymmetric_normalComponentFlips)
{
    ManagerFixture fix;
    auto const& layout   = fix.layout;
    auto const& meshData = fix.manager.getMeshData();

    constexpr double Cx = 3.0, Cy = 4.0, Cz = 5.0;

    VecFieldMHD2 rhoV{"rhoV", layout, PHARE::core::MHDQuantity::Vector::rhoV};
    auto& rhoVx = rhoV.getComponent(PHARE::core::Component::X);
    auto& rhoVy = rhoV.getComponent(PHARE::core::Component::Y);
    auto& rhoVz = rhoV.getComponent(PHARE::core::Component::Z);

    for (auto i = 0u; i < rhoVx.shape()[0]; ++i)
        for (auto j = 0u; j < rhoVx.shape()[1]; ++j)
        {
            rhoVx(i, j) = Cx;
            rhoVy(i, j) = Cy;
            rhoVz(i, j) = Cz;
        }

    DummyState state;
    PHARE::core::InnerBCContext<DummyState> ctx{state, state, 0.0};
    fix.manager.applyBC(rhoV, layout, ctx);

    constexpr std::array<PHARE::core::QtyCentering, 2> kCellC
        = {PHARE::core::QtyCentering::dual, PHARE::core::QtyCentering::dual};

    bool foundInPatch = false;
    for (auto const& g : meshData.getGhostDataFromCentering(kCellC))
    {
        if (!g.mirrorIsInterpolable)
            continue;
        foundInPatch = true;

        EXPECT_NEAR(rhoVx(g.index), -Cx, eps)
            << "rhoVx ghost at (" << g.index[0] << "," << g.index[1] << ")";
        EXPECT_NEAR(rhoVy(g.index), Cy, eps)
            << "rhoVy ghost at (" << g.index[0] << "," << g.index[1] << ")";
        EXPECT_NEAR(rhoVz(g.index), Cz, eps)
            << "rhoVz ghost at (" << g.index[0] << "," << g.index[1] << ")";
    }
    EXPECT_TRUE(foundInPatch) << "at least one in-patch ghost must exist";
}

/**
 * @brief Vector E under Reflective gets Antisymmetric: tangential component is reversed.
 *
 * Constant field (Ex, Ey, Ez), plane normal n=(1,0,0).
 * Antisymmetric: keep normal (Ex), reverse tangential → ghost = (Ex, -Ey, -Ez).
 *
 * But E components are edge-centred, so the ghost elements come from ghostEdgesData,
 * not ghostCellsData.  We verify the x-component (EdgeCenteredX centering = (dual,primal,primal)
 * → ghostEdgesData[0]) separately from y/z (EdgeCenteredY/Z).
 */
TEST(InnerBoundaryManager, reflectiveEIsAntisymmetric_tangentialComponentFlips)
{
    ManagerFixture fix;
    auto const& layout   = fix.layout;
    auto const& meshData = fix.manager.getMeshData();

    constexpr double Ex = 2.0, Ey = 3.0, Ez = -1.0;

    VecFieldMHD2 E{"E", layout, PHARE::core::MHDQuantity::Vector::E};
    auto& Ex_field = E.getComponent(PHARE::core::Component::X);
    auto& Ey_field = E.getComponent(PHARE::core::Component::Y);
    auto& Ez_field = E.getComponent(PHARE::core::Component::Z);

    for (auto i = 0u; i < Ex_field.shape()[0]; ++i)
        for (auto j = 0u; j < Ex_field.shape()[1]; ++j)
            Ex_field(i, j) = Ex;

    for (auto i = 0u; i < Ey_field.shape()[0]; ++i)
        for (auto j = 0u; j < Ey_field.shape()[1]; ++j)
            Ey_field(i, j) = Ey;

    // Ez is the out-of-plane component; fill it but its ghost list has no separate entry in 2D.
    for (auto i = 0u; i < Ez_field.shape()[0]; ++i)
        for (auto j = 0u; j < Ez_field.shape()[1]; ++j)
            Ez_field(i, j) = Ez;

    DummyState state;
    PHARE::core::InnerBCContext<DummyState> ctx{state, state, 0.0};
    fix.manager.applyBC(E, layout, ctx);

    // In 2D, EdgeCenteredX has centering (dual,primal) and EdgeCenteredY has (primal,dual).
    // EdgeCenteredZ collapses to all-primal = node centering, shared with Ez.
    // Normal n = (1,0,0):
    //   Ex (dual,primal): antisymmetric keeps normal component → Ex preserved.
    //   Ey (primal,dual): antisymmetric reverses tangential → Ey negated.
    constexpr std::array<PHARE::core::QtyCentering, 2> kEdgeXC // EdgeCenteredX
        = {PHARE::core::QtyCentering::dual, PHARE::core::QtyCentering::primal};
    constexpr std::array<PHARE::core::QtyCentering, 2> kEdgeYC // EdgeCenteredY
        = {PHARE::core::QtyCentering::primal, PHARE::core::QtyCentering::dual};
    constexpr std::array<PHARE::core::QtyCentering, 2> kEdgeZC // EdgeCenteredZ (2D: all-primal)
        = {PHARE::core::QtyCentering::primal, PHARE::core::QtyCentering::primal};

    // Antisymmetric BC skips a ghost when any sibling centering fails interpolability at
    // the mirror point (see field_antisymmetric_inner_boundary_condition.hpp). Mirror that
    // logic here: only assert BC behaviour on ghosts where all three sibling centerings
    // accept the mirror.
    auto siblings_ok = [&](auto const& g, auto const& sibA, auto const& sibB) {
        return PHARE::core::FieldAtPoint<2, 1>::pointIsInterpolable(layout, g.mirrorPoint, sibA)
               && PHARE::core::FieldAtPoint<2, 1>::pointIsInterpolable(layout, g.mirrorPoint,
                                                                       sibB);
    };

    bool foundAnyEdge = false;
    for (auto const& g : meshData.getGhostDataFromCentering(kEdgeXC))
    {
        if (!g.mirrorIsInterpolable)
            continue;
        if (!siblings_ok(g, kEdgeYC, kEdgeZC))
            continue;
        foundAnyEdge = true;
        EXPECT_NEAR(Ex_field(g.index), Ex, eps) << "Ex (normal component) should be preserved";
    }
    for (auto const& g : meshData.getGhostDataFromCentering(kEdgeYC))
    {
        if (!g.mirrorIsInterpolable)
            continue;
        if (!siblings_ok(g, kEdgeXC, kEdgeZC))
            continue;
        foundAnyEdge = true;
        EXPECT_NEAR(Ey_field(g.index), -Ey, eps) << "Ey (tangential component) should be negated";
    }

    EXPECT_TRUE(foundAnyEdge) << "no in-patch edge ghost found — check classifier / geometry";
}

/**
 * @brief applyBC is a no-op for a quantity that was not registered.
 */
TEST(InnerBoundaryManager, applyBCIsNoopForUnregisteredScalar)
{
    ManagerFixture fix;
    auto const& layout = fix.layout;

    PHARE::core::NdArrayVector<2, double> storage{
        layout.allocSize(PHARE::core::MHDQuantity::Scalar::P)};
    ScalarField field{"P", PHARE::core::MHDQuantity::Scalar::P, storage.data(), storage.shape()};

    constexpr double C = 42.0;
    for (auto i = 0u; i < field.shape()[0]; ++i)
        for (auto j = 0u; j < field.shape()[1]; ++j)
            field(i, j) = C;

    DummyState state;
    PHARE::core::InnerBCContext<DummyState> ctx{state, state, 0.0};
    // P was not registered — field.physicalQuantity() lookup misses, applyBC must no-op.
    fix.manager.applyBC(field, layout, ctx);

    for (auto i = 0u; i < field.shape()[0]; ++i)
        for (auto j = 0u; j < field.shape()[1]; ++j)
            EXPECT_NEAR(field(i, j), C, eps) << "unregistered scalar field was modified";
}

// ---------------------------------------------------------------------------
// setSafeState — pins inactive cells to the prescribed/derived safe state
// ---------------------------------------------------------------------------

namespace
{
// Count inactive cells (and assert at least one exists in this geometry).
template<typename MeshDataT, typename FieldT>
std::size_t expectHasInactive(MeshDataT const& md, FieldT const& shaped)
{
    auto const& cellStatus = md.cellStatusField();
    std::size_t n          = 0;
    for (auto i = 0u; i < shaped.shape()[0]; ++i)
        for (auto j = 0u; j < shaped.shape()[1]; ++j)
            if (cellStatus(i, j) == PHARE::core::toDouble(PHARE::core::ElemStatus::Inactive))
                ++n;
    EXPECT_GT(n, 0u) << "geometry must produce at least one inactive cell";
    return n;
}
} // namespace

/**
 * @brief setSafeState pins inactive cells to the prescribed scalar value, leaving others.
 */
TEST(InnerBoundaryManager, setSafeStateScalarPinsInactiveOnly)
{
    auto thermo = std::make_shared<PHARE::core::IdealGasThermo>(5.0 / 3.0);
    Manager::SafeStateValues safe;
    safe.density = 1.5;
    ManagerFixture fix{thermo, safe};
    auto const& layout   = fix.layout;
    auto const& meshData = fix.manager.getMeshData();

    constexpr double sentinel = 7.0;
    PHARE::core::NdArrayVector<2, double> storage{
        layout.allocSize(PHARE::core::MHDQuantity::Scalar::rho)};
    ScalarField rho{"rho", PHARE::core::MHDQuantity::Scalar::rho, storage.data(), storage.shape()};
    for (auto i = 0u; i < rho.shape()[0]; ++i)
        for (auto j = 0u; j < rho.shape()[1]; ++j)
            rho(i, j) = sentinel;

    expectHasInactive(meshData, rho);
    fix.manager.setSafeState(rho, layout);

    auto const& cellStatus = meshData.cellStatusField();
    for (auto i = 0u; i < rho.shape()[0]; ++i)
        for (auto j = 0u; j < rho.shape()[1]; ++j)
        {
            bool inactive
                = cellStatus(i, j) == PHARE::core::toDouble(PHARE::core::ElemStatus::Inactive);
            EXPECT_NEAR(rho(i, j), inactive ? 1.5 : sentinel, eps)
                << "cell (" << i << "," << j << ") inactive=" << inactive;
        }
}

/**
 * @brief setSafeState derives momentum (rhoV = density * velocity) for the vector path.
 */
TEST(InnerBoundaryManager, setSafeStateVectorDerivesMomentum)
{
    auto thermo = std::make_shared<PHARE::core::IdealGasThermo>(5.0 / 3.0);
    Manager::SafeStateValues safe;
    safe.density  = 1.5;
    safe.velocity = {2.0, 0.0, 0.0}; // rhoV safe = {3, 0, 0}
    ManagerFixture fix{thermo, safe};
    auto const& layout   = fix.layout;
    auto const& meshData = fix.manager.getMeshData();

    constexpr double sentinel = 9.0;
    VecFieldMHD2 rhoV{"rhoV", layout, PHARE::core::MHDQuantity::Vector::rhoV};
    for (auto comp : {PHARE::core::Component::X, PHARE::core::Component::Y, PHARE::core::Component::Z})
    {
        auto& c = rhoV.getComponent(comp);
        for (auto i = 0u; i < c.shape()[0]; ++i)
            for (auto j = 0u; j < c.shape()[1]; ++j)
                c(i, j) = sentinel;
    }

    fix.manager.setSafeState(rhoV, layout);

    auto const& cellStatus = meshData.cellStatusField();
    auto& rhoVx            = rhoV.getComponent(PHARE::core::Component::X);
    auto& rhoVy            = rhoV.getComponent(PHARE::core::Component::Y);
    for (auto i = 0u; i < rhoVx.shape()[0]; ++i)
        for (auto j = 0u; j < rhoVx.shape()[1]; ++j)
        {
            bool inactive
                = cellStatus(i, j) == PHARE::core::toDouble(PHARE::core::ElemStatus::Inactive);
            EXPECT_NEAR(rhoVx(i, j), inactive ? 3.0 : sentinel, eps);
            EXPECT_NEAR(rhoVy(i, j), inactive ? 0.0 : sentinel, eps);
        }
}

/**
 * @brief setSafeState derives Etot from the safe primitives via the equation of state.
 *
 * Ideal gas, gamma=5/3, V=0, B=0 ⇒ Etot = P/(gamma-1) = 2 / (2/3) = 3.
 */
TEST(InnerBoundaryManager, setSafeStateEtotDerived)
{
    auto thermo = std::make_shared<PHARE::core::IdealGasThermo>(5.0 / 3.0);
    Manager::SafeStateValues safe;
    safe.density  = 1.5;
    safe.pressure = 2.0;
    ManagerFixture fix{thermo, safe};
    auto const& layout   = fix.layout;
    auto const& meshData = fix.manager.getMeshData();

    PHARE::core::NdArrayVector<2, double> storage{
        layout.allocSize(PHARE::core::MHDQuantity::Scalar::Etot)};
    ScalarField etot{"Etot", PHARE::core::MHDQuantity::Scalar::Etot, storage.data(),
                     storage.shape()};
    for (auto i = 0u; i < etot.shape()[0]; ++i)
        for (auto j = 0u; j < etot.shape()[1]; ++j)
            etot(i, j) = 9.0;

    fix.manager.setSafeState(etot, layout);

    auto const& cellStatus = meshData.cellStatusField();
    for (auto i = 0u; i < etot.shape()[0]; ++i)
        for (auto j = 0u; j < etot.shape()[1]; ++j)
        {
            if (cellStatus(i, j) == PHARE::core::toDouble(PHARE::core::ElemStatus::Inactive))
            {
                EXPECT_NEAR(etot(i, j), 3.0, eps)
                    << "derived Etot mismatch at (" << i << "," << j << ")";
            }
        }
}

/**
 * @brief setSafeState throws for a quantity that has no registered safe value.
 */
TEST(InnerBoundaryManager, setSafeStateThrowsForUnknownQuantity)
{
    ManagerFixture fix; // null thermo, default safe state
    auto const& layout = fix.layout;

    // Etot (total energy) is not among the safe scalars (rho, P, Etot).
    PHARE::core::NdArrayVector<2, double> storage{
        layout.allocSize(PHARE::core::MHDQuantity::Scalar::Etot)};
    ScalarField etot{"Etot", PHARE::core::MHDQuantity::Scalar::Etot, storage.data(),
                     storage.shape()};
    EXPECT_THROW(fix.manager.setSafeState(etot, layout), std::runtime_error);

    // Etot has no safe value when no thermo was supplied (derivation skipped).
    PHARE::core::NdArrayVector<2, double> storage2{
        layout.allocSize(PHARE::core::MHDQuantity::Scalar::Etot)};
    ScalarField etot1{"Etot", PHARE::core::MHDQuantity::Scalar::Etot, storage2.data(),
                      storage2.shape()};
    EXPECT_THROW(fix.manager.setSafeState(etot1, layout), std::runtime_error);
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
