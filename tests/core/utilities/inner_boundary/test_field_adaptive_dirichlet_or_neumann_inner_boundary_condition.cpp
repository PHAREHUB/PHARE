#include "gtest/gtest.h"

#include <cstdint>
#include <vector>

#include "core/data/field/field.hpp"
#include "core/data/grid/gridlayout.hpp"
#include "core/data/grid/gridlayoutimplyee_mhd.hpp"
#include "core/data/ndarray/ndarray_vector.hpp"
#include "core/inner_boundary/field_adaptive_dirichlet_or_neumann_inner_boundary_condition.hpp"
#include "core/inner_boundary/inner_boundary_mesh_classifier.hpp"
#include "core/inner_boundary/plane_inner_boundary.hpp"
#include "core/mhd/mhd_quantities.hpp"
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

/// Minimal MHD-like state: a scalar field to fill (rho) + a vector criterion (rhoV). Exposes the
/// getVector(MHDQuantity::Vector) accessor the adaptive BC uses to resolve its criterion.
struct MhdState
{
    NdArrayVector<2, double> rhoStore;
    ScalarField              rho;
    UsableVecFieldMHD<2>     rhoV;

    explicit MhdState(GridLayout const& layout)
        : rhoStore{layout.allocSize(MHDQuantity::Scalar::rho)}
        , rho{"rho", MHDQuantity::Scalar::rho, rhoStore.data(), rhoStore.shape()}
        , rhoV{"rhoV", layout, MHDQuantity::Vector::rhoV}
    {
    }

    auto& getVector(MHDQuantity::Vector q)
    {
        if (q != MHDQuantity::Vector::rhoV)
            throw std::runtime_error("test MhdState::getVector only supports rhoV");
        return rhoV;
    }
};

/// Allocate and wire all mesh-data buffers (same pattern as the other inner-BC tests).
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

    NdArrayVector<2, double>                  phi_storage;
    std::vector<NdArrayVector<2, double>>     elem_storages;
    GhostElemPack<2>::ghost_elem_array_type   ghost_array{};
    GhostElemPack<2>::ghost_elem_array_type   degraded_array{};
    MeshData                                  tags;
};

/// Plane at x = 0, normal (1,0): fluid is x>0, ghosts are x<0. Same 4×2 geometry as the
/// other inner-BC tests. The per-ghost outward normal (ghost→mirror) is therefore (1,0),
/// so criterion·n == rhoVx at the boundary.
struct Fixture
{
    PlaneInnerBoundary<2> plane{"plane", {0.0, 0.0}, {1.0, 0.0}};
    Box<int, 2>           amr_box{{-2, 0}, {1, 1}};
    GridLayout            layout{{1.0, 1.0}, {4u, 2u}, {0.0, 0.0}, amr_box};

    MeshDataBuffers buffers{layout};
    MhdState        state{layout};

    Fixture()
    {
        Classifier::Overrides ov;
        ov.cut_eps      = 1e-12;
        ov.inactive_eps = 0.0;
        auto classifier = Classifier::withDefaults(plane, layout, ov);
        classifier(layout, buffers.tags);
    }

    // fill rho everywhere with the linear field f(x,y) = x + y, so linear interpolation at any
    // point reproduces it exactly and the mirror value is predictable.
    void fillRhoLinear()
    {
        auto& rho = state.rho;
        for (auto i = 0u; i < rho.shape()[0]; ++i)
            for (auto j = 0u; j < rho.shape()[1]; ++j)
            {
                auto amr_pos = layout.localToAMR(PHARE::core::Point<std::uint32_t, 2>{i, j});
                auto amr_ij  = PHARE::core::Point<int, 2>{static_cast<int>(amr_pos[0]),
                                                          static_cast<int>(amr_pos[1])};
                auto pos     = layout.fieldNodeCoordinates(rho, amr_ij);
                rho(i, j)    = pos[0] + pos[1];
            }
    }

    template<typename F>
    static void fillAll(F& field, double value)
    {
        for (auto i = 0u; i < field.shape()[0]; ++i)
            for (auto j = 0u; j < field.shape()[1]; ++j)
                field(i, j) = value;
    }

    auto makeBC(double dirichletValue)
    {
        return FieldAdaptiveDirichletOrNeumannInnerBoundaryCondition<ScalarField, GridLayout,
                                                                     MhdState>{
            MHDQuantity::Vector::rhoV, dirichletValue};
    }
};

} // namespace

// ---------------------------------------------------------------------------
// rhoV·n > 0 (flux into the fluid domain) selects Dirichlet on every ghost.
// ---------------------------------------------------------------------------
TEST(FieldAdaptiveDirichletOrNeumannInner, positiveNormalFluxGivesDirichlet)
{
    Fixture fix;
    auto&   st     = fix.state;
    auto&   layout = fix.layout;

    constexpr double boundaryValue = 5.0;

    fix.fillRhoLinear();
    Fixture::fillAll(st.rhoV[0], +1.0); // rhoVx > 0  → rhoV·n > 0
    Fixture::fillAll(st.rhoV[1], 0.0);
    Fixture::fillAll(st.rhoV[2], 0.0);

    auto bc = fix.makeBC(boundaryValue);
    bc.apply(st.rho, layout, fix.buffers.tags, InnerBCContext<MhdState>{st, st, 0.0});

    bool foundInPatch = false;
    for (auto const& g : fix.buffers.tags.getGhostDataFromCentering(kCellC))
    {
        if (!g.mirrorIsInterpolable)
            continue;
        foundInPatch          = true;
        double const mirror   = g.mirrorPoint[0] + g.mirrorPoint[1];
        double const expected = 2.0 * boundaryValue - mirror; // Dirichlet
        EXPECT_NEAR(st.rho(g.index), expected, 1e-10)
            << "ghost (" << g.index[0] << "," << g.index[1] << ")";
    }
    EXPECT_TRUE(foundInPatch) << "at least one in-patch ghost must exist";
}

// ---------------------------------------------------------------------------
// rhoV·n < 0 (flux into the body) selects Neumann (zero-gradient) on every ghost.
// ---------------------------------------------------------------------------
TEST(FieldAdaptiveDirichletOrNeumannInner, negativeNormalFluxGivesNeumann)
{
    Fixture fix;
    auto&   st     = fix.state;
    auto&   layout = fix.layout;

    fix.fillRhoLinear();
    Fixture::fillAll(st.rhoV[0], -1.0); // rhoVx < 0  → rhoV·n < 0
    Fixture::fillAll(st.rhoV[1], 0.0);
    Fixture::fillAll(st.rhoV[2], 0.0);

    auto bc = fix.makeBC(5.0);
    bc.apply(st.rho, layout, fix.buffers.tags, InnerBCContext<MhdState>{st, st, 0.0});

    bool foundInPatch = false;
    for (auto const& g : fix.buffers.tags.getGhostDataFromCentering(kCellC))
    {
        if (!g.mirrorIsInterpolable)
            continue;
        foundInPatch          = true;
        double const expected = g.mirrorPoint[0] + g.mirrorPoint[1]; // Neumann == mirror value
        EXPECT_NEAR(st.rho(g.index), expected, 1e-10)
            << "ghost (" << g.index[0] << "," << g.index[1] << ")";
    }
    EXPECT_TRUE(foundInPatch) << "at least one in-patch ghost must exist";
}

// ---------------------------------------------------------------------------
// The switch is per ghost element: a criterion whose sign varies in y must give
// Dirichlet on some ghosts and Neumann on others, decided by the local rhoVx sign.
// rhoVx is uniform in x, so the value interpolated at the surface equals the value
// stored at the ghost's own row — which we read back to derive the expectation.
// ---------------------------------------------------------------------------
TEST(FieldAdaptiveDirichletOrNeumannInner, perElementSwitchOnCriterionSign)
{
    Fixture fix;
    auto&   st     = fix.state;
    auto&   layout = fix.layout;

    constexpr double boundaryValue = 5.0;

    fix.fillRhoLinear();
    Fixture::fillAll(st.rhoV[1], 0.0);
    Fixture::fillAll(st.rhoV[2], 0.0);
    // row j=0 negative (Neumann), row j=1 positive (Dirichlet); uniform in x.
    auto& rhoVx = st.rhoV[0];
    for (auto i = 0u; i < rhoVx.shape()[0]; ++i)
    {
        rhoVx(i, 0u) = -1.0;
        rhoVx(i, 1u) = +1.0;
    }

    auto bc = fix.makeBC(boundaryValue);
    bc.apply(st.rho, layout, fix.buffers.tags, InnerBCContext<MhdState>{st, st, 0.0});

    bool sawDirichlet = false, sawNeumann = false;
    for (auto const& g : fix.buffers.tags.getGhostDataFromCentering(kCellC))
    {
        if (!g.mirrorIsInterpolable)
            continue;

        double const crit   = rhoVx(g.index); // x-uniform ⇒ equals the interpolated surface value
        double const mirror = g.mirrorPoint[0] + g.mirrorPoint[1];
        if (crit > 0.0)
        {
            sawDirichlet = true;
            EXPECT_NEAR(st.rho(g.index), 2.0 * boundaryValue - mirror, 1e-10)
                << "ghost (" << g.index[0] << "," << g.index[1] << ") expected Dirichlet";
        }
        else
        {
            sawNeumann = true;
            EXPECT_NEAR(st.rho(g.index), mirror, 1e-10)
                << "ghost (" << g.index[0] << "," << g.index[1] << ") expected Neumann";
        }
    }
    EXPECT_TRUE(sawDirichlet) << "expected at least one Dirichlet (rhoVx>0) ghost";
    EXPECT_TRUE(sawNeumann) << "expected at least one Neumann (rhoVx<0) ghost";
}

// ---------------------------------------------------------------------------
// Non-interpolable ghosts are left untouched.
// ---------------------------------------------------------------------------
TEST(FieldAdaptiveDirichletOrNeumannInner, nonInterpolableGhostsUntouched)
{
    Fixture fix;
    auto&   st     = fix.state;
    auto&   layout = fix.layout;

    constexpr double sentinel = -999.0;

    fix.fillRhoLinear();
    Fixture::fillAll(st.rhoV[0], +1.0);
    Fixture::fillAll(st.rhoV[1], 0.0);
    Fixture::fillAll(st.rhoV[2], 0.0);

    for (auto const& g : fix.buffers.tags.getGhostDataFromCentering(kCellC))
        if (!g.mirrorIsInterpolable)
            st.rho(g.index) = sentinel;

    auto bc = fix.makeBC(5.0);
    bc.apply(st.rho, layout, fix.buffers.tags, InnerBCContext<MhdState>{st, st, 0.0});

    for (auto const& g : fix.buffers.tags.getGhostDataFromCentering(kCellC))
    {
        if (!g.mirrorIsInterpolable)
        {
            EXPECT_NEAR(st.rho(g.index), sentinel, eps)
                << "non-interpolable ghost (" << g.index[0] << "," << g.index[1] << ") changed";
        }
    }
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
