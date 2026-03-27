#include "gtest/gtest.h"

#include "core/inner_boundary/inner_boundary_factory.hpp"
#include "core/inner_boundary/plane_inner_boundary.hpp"
#include "core/inner_boundary/sphere_inner_boundary.hpp"

#include "initializer/data_provider.hpp"

namespace
{
constexpr double eps = 1e-12;
}

using namespace PHARE;

TEST(InnerBoundarySphere, basicGeometry)
{
    core::SphereInnerBoundary<2> sphere{"sphere", {0., 0.}, 1.};

    auto const p = core::Point<double, 2>{2., 0.};
    EXPECT_NEAR(sphere.signedDistance(p), 1., eps);

    auto const n = sphere.normal(p);
    EXPECT_NEAR(n[0], 1., eps);
    EXPECT_NEAR(n[1], 0., eps);

    auto const proj = sphere.project(p);
    EXPECT_NEAR(proj[0], 1., eps);
    EXPECT_NEAR(proj[1], 0., eps);

    auto const sym = sphere.symmetric(p);
    EXPECT_NEAR(sym[0], 0., eps);
    EXPECT_NEAR(sym[1], 0., eps);
}

TEST(InnerBoundaryPlane, basicGeometry)
{
    core::PlaneInnerBoundary<3> plane{"plane", {0., 0., 0.}, {0., 0., 2.}};

    auto const p = core::Point<double, 3>{1., 2., 3.};
    EXPECT_NEAR(plane.signedDistance(p), 3., eps);

    auto const n = plane.normal(p);
    EXPECT_NEAR(n[0], 0., eps);
    EXPECT_NEAR(n[1], 0., eps);
    EXPECT_NEAR(n[2], 1., eps);

    auto const proj = plane.project(p);
    EXPECT_NEAR(proj[0], 1., eps);
    EXPECT_NEAR(proj[1], 2., eps);
    EXPECT_NEAR(proj[2], 0., eps);

    auto const sym = plane.symmetric(p);
    EXPECT_NEAR(sym[0], 1., eps);
    EXPECT_NEAR(sym[1], 2., eps);
    EXPECT_NEAR(sym[2], -3., eps);
}

TEST(InnerBoundaryParser, createsSphere)
{
    initializer::PHAREDict dict;
    dict["inner_boundary"]["name"]   = std::string{"my_sphere"};
    dict["inner_boundary"]["shape"]  = std::string{"sphere"};
    dict["inner_boundary"]["center"] = std::vector<double>{1., 2.};
    dict["inner_boundary"]["radius"] = 2.;

    auto boundary = core::InnerBoundaryFactory<2>::create(dict);
    ASSERT_NE(boundary, nullptr);

    auto const d = boundary->signedDistance(core::Point<double, 2>{1., 6.});
    EXPECT_NEAR(d, 2., eps);
}

TEST(InnerBoundaryParser, createsPlane)
{
    initializer::PHAREDict dict;
    dict["inner_boundary"]["name"]   = std::string{"my_plane"};
    dict["inner_boundary"]["shape"]  = std::string{"plane"};
    dict["inner_boundary"]["point"]  = std::vector<double>{0., 0., 0.};
    dict["inner_boundary"]["normal"] = std::vector<double>{0., 1., 0.};

    auto boundary = core::InnerBoundaryFactory<3>::create(dict);
    ASSERT_NE(boundary, nullptr);

    auto const d = boundary->signedDistance(core::Point<double, 3>{0., -4., 0.});
    EXPECT_NEAR(d, -4., eps);
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
