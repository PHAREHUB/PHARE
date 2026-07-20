#include "gtest/gtest.h"

#include <array>
#include <cmath>
#include <cstdint>
#include <vector>

#include "core/data/field/field.hpp"
#include "core/data/grid/gridlayout.hpp"
#include "core/data/grid/gridlayoutimplyee_mhd.hpp"
#include "core/data/ndarray/ndarray_vector.hpp"
#include "core/inner_boundary/field_ionospheric_convection_momentum_inner_boundary_condition.hpp"
#include "core/inner_boundary/inner_boundary_mesh_classifier.hpp"
#include "core/inner_boundary/sphere_inner_boundary.hpp"
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
using VecFieldMHD2   = UsableVecFieldMHD<2>;
using VecF           = VecField<ScalarField, MHDQuantity>;

constexpr std::array<QtyCentering, 2> kCellC = {QtyCentering::dual, QtyCentering::dual};

/// Minimal MHD-like state exposing the conserved momentum + total-field components the
/// ionospheric-convection momentum BC reads (rhoV, B).
struct MhdState
{
    VecFieldMHD2 rhoV;
    VecFieldMHD2 B;

    explicit MhdState(GridLayout const& layout)
        : rhoV{"rhoV", layout, MHDQuantity::Vector::rhoV}
        , B{"B", layout, MHDQuantity::Vector::B}
    {
    }
};

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

    NdArrayVector<2, double>                phi_storage;
    std::vector<NdArrayVector<2, double>>   elem_storages;
    GhostElemPack<2>::ghost_elem_array_type ghost_array{};
    GhostElemPack<2>::ghost_elem_array_type degraded_array{};
    MeshData                                tags;
};

/// Sphere centred in a 16x16 grid (dx=1) so the classifier produces a ring of cell ghosts
/// with varied radial normals.
struct Fixture
{
    // off any node (integer) and cell-centre (half-integer) so geometry.normal is never queried
    // at the sphere centre during classification.
    PHARE::core::Point<double, 2> center{8.3, 8.3};
    double                        radius{4.0};
    SphereInnerBoundary<2>        sphere{"sphere", center, radius};
    Box<int, 2>                   amr_box{{0, 0}, {15, 15}};
    GridLayout                    layout{{1.0, 1.0}, {16u, 16u}, {0.0, 0.0}, amr_box};

    MeshDataBuffers buffers{layout};
    MhdState        state{layout};

    Fixture()
    {
        Classifier::Overrides ov;
        ov.cut_eps      = 1e-12;
        ov.inactive_eps = 0.0;
        auto classifier = Classifier::withDefaults(sphere, layout, ov);
        classifier(layout, buffers.tags);
    }

    static void fillAll(ScalarField& f, double v)
    {
        for (auto i = 0u; i < f.shape()[0]; ++i)
            for (auto j = 0u; j < f.shape()[1]; ++j)
                f(i, j) = v;
    }

    auto makeBC()
    {
        return FieldIonosphericConvectionMomentumInnerBoundaryCondition<VecF, GridLayout,
                                                                         MhdState>{&sphere};
    }

    PHARE::core::Point<double, 2> ghostCoord(GhostElemData<2> const& g) const
    {
        auto amr = layout.localToAMR(g.index);
        PHARE::core::Point<int, 2> idx{static_cast<int>(amr[0]), static_cast<int>(amr[1])};
        return layout.fieldNodeCoordinates(state.rhoV[0], idx);
    }
};

} // namespace

// ---------------------------------------------------------------------------
// Field-aligned + r^2 rhoV_n conservation (constant rhoV, uniform B).
// ---------------------------------------------------------------------------
TEST(FieldIonosphericConvectionMomentumInner, fieldAlignedRadialScaling)
{
    Fixture fix;
    auto&   st     = fix.state;
    auto&   layout = fix.layout;

    std::array<double, 3> const rhoV0{2.0, -1.0, 0.5};
    std::array<double, 3> const Bvec{1.0, 2.0, 0.0}; // b̂ in-plane, generally not ⟂ to normals
    double const                Bnorm = std::sqrt(Bvec[0] * Bvec[0] + Bvec[1] * Bvec[1]);
    std::array<double, 3> const bhat{Bvec[0] / Bnorm, Bvec[1] / Bnorm, 0.0};

    for (int c = 0; c < 3; ++c)
    {
        Fixture::fillAll(st.rhoV[c], rhoV0[c]);
        Fixture::fillAll(st.B[c], Bvec[c]);
    }

    auto bc = fix.makeBC();
    bc.apply(st.rhoV, layout, fix.buffers.tags, InnerBCContext<MhdState>{st, st, 0.0});

    bool found = false;
    for (auto const& g : fix.buffers.tags.getGhostDataFromCentering(kCellC))
    {
        if (!g.mirrorIsInterpolable)
            continue;
        // skip near-tangent ghosts: the BC caps the 1/(b̂·n) Tanaka gain by falling back to the
        // symmetric condition for |b̂·n| < eps_bn (=0.1, see the BC); match that threshold here.
        double bn = bhat[0] * g.normal[0] + bhat[1] * g.normal[1];
        if (std::abs(bn) < 0.1)
            continue;

        auto gc        = fix.ghostCoord(g);
        double r_g     = std::hypot(gc[0] - fix.center[0], gc[1] - fix.center[1]);
        double r_m     = std::hypot(g.mirrorPoint[0] - fix.center[0],
                                    g.mirrorPoint[1] - fix.center[1]);
        double rhoVn_m = rhoV0[0] * g.normal[0] + rhoV0[1] * g.normal[1]; // n is in-plane
        double ratio   = r_m / r_g;
        // rhoV_g = (1 + (r_m/r_g)^2) (rhoV_m·n / (b̂·n)) b̂ - rhoV_m
        double cexp = (1.0 + ratio * ratio) * (rhoVn_m / bn);

        found = true;
        for (int c = 0; c < 3; ++c)
            EXPECT_NEAR(st.rhoV[c](g.index), cexp * bhat[c] - rhoV0[c], 1e-9)
                << "comp " << c << " ghost (" << g.index[0] << "," << g.index[1] << ")";
    }
    EXPECT_TRUE(found) << "expected at least one non-tangent cell ghost";
}

// ---------------------------------------------------------------------------
// Negligible field -> symmetric fallback (reverse normal component, keep tangential).
// ---------------------------------------------------------------------------
TEST(FieldIonosphericConvectionMomentumInner, zeroFieldFallsBackToSymmetric)
{
    Fixture fix;
    auto&   st     = fix.state;
    auto&   layout = fix.layout;

    std::array<double, 3> const rhoV0{2.0, -1.0, 0.5};
    for (int c = 0; c < 3; ++c)
    {
        Fixture::fillAll(st.rhoV[c], rhoV0[c]);
        Fixture::fillAll(st.B[c], 0.0);
    }

    auto bc = fix.makeBC();
    bc.apply(st.rhoV, layout, fix.buffers.tags, InnerBCContext<MhdState>{st, st, 0.0});

    bool found = false;
    for (auto const& g : fix.buffers.tags.getGhostDataFromCentering(kCellC))
    {
        if (!g.mirrorIsInterpolable)
            continue;
        double rhoVn = rhoV0[0] * g.normal[0] + rhoV0[1] * g.normal[1];
        // symmetric: ghost = rhoV - 2 (rhoV·n) n   (n in-plane, z untouched)
        std::array<double, 3> expected{rhoV0[0] - 2.0 * rhoVn * g.normal[0],
                                       rhoV0[1] - 2.0 * rhoVn * g.normal[1], rhoV0[2]};
        found = true;
        for (int c = 0; c < 3; ++c)
            EXPECT_NEAR(st.rhoV[c](g.index), expected[c], 1e-9) << "comp " << c;
    }
    EXPECT_TRUE(found);
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
