#include "gtest/gtest.h"

#include "core/boundary/boundary_inflow_compose.hpp"
#include "initializer/data_provider.hpp"

#include <memory>
#include <vector>

using namespace PHARE::core::inflow_compose;
using PHARE::core::Span;
using PHARE::core::VectorSpan;
using PHARE::initializer::SpaceTimeFunction;

namespace
{
// A 1D space-time function f(x, t) built from a plain lambda over the batch.
SpaceTimeFunction<1> make1D(std::function<double(double, double)> f)
{
    return [f](std::vector<double> const& x, double t) -> std::shared_ptr<Span<double>> {
        std::vector<double> out(x.size());
        for (std::size_t k = 0; k < x.size(); ++k)
            out[k] = f(x[k], t);
        return std::make_shared<VectorSpan<double>>(std::move(out));
    };
}
} // namespace

TEST(InflowCompose, ConstFunctionBroadcastsToNodeCount)
{
    auto c = constFunction<1>(2.5);
    std::vector<double> x{0.0, 1.0, 2.0, 3.0};
    auto s = c(x, 7.0); // time is ignored
    ASSERT_EQ(s->size(), x.size());
    for (std::size_t k = 0; k < x.size(); ++k)
        EXPECT_DOUBLE_EQ((*s)[k], 2.5);
}

TEST(InflowCompose, MulFunctionMultipliesElementwise)
{
    auto f = make1D([](double x, double) { return x; });
    auto g = make1D([](double, double t) { return t; });
    auto p = mulFunction<1>(f, g); // x * t
    std::vector<double> x{1.0, 2.0, 4.0};
    auto s = p(x, 3.0);
    ASSERT_EQ(s->size(), x.size());
    EXPECT_DOUBLE_EQ((*s)[0], 3.0);
    EXPECT_DOUBLE_EQ((*s)[1], 6.0);
    EXPECT_DOUBLE_EQ((*s)[2], 12.0);
}

TEST(InflowCompose, ProdComb2LinearlyCombinesTwoProducts)
{
    auto one = constFunction<1>(1.0);
    auto x   = make1D([](double x, double) { return x; });
    auto y2  = make1D([](double, double) { return 2.0; });
    // a*f1*g1 + b*f2*g2 = (-1)*x*1 + (3)*2*1 = -x + 6
    auto r = prodComb2<1>(-1.0, x, one, 3.0, y2, one);
    std::vector<double> xs{0.0, 1.0, 5.0};
    auto s = r(xs, 0.0);
    EXPECT_DOUBLE_EQ((*s)[0], 6.0);
    EXPECT_DOUBLE_EQ((*s)[1], 5.0);
    EXPECT_DOUBLE_EQ((*s)[2], 1.0);
}

TEST(InflowCompose, NegCrossMatchesMinusVCrossB)
{
    // Uniform V = (1,2,3), B = (4,5,6); E = -V x B = -(2*6-3*5, 3*4-1*6, 1*5-2*4)
    //                                             = -(-3, 6, -3) = (3, -6, 3)
    std::array<SpaceTimeFunction<1>, 3> V{constFunction<1>(1.0), constFunction<1>(2.0),
                                          constFunction<1>(3.0)};
    std::array<SpaceTimeFunction<1>, 3> B{constFunction<1>(4.0), constFunction<1>(5.0),
                                          constFunction<1>(6.0)};
    auto E = negCrossFunction<1>(V, B);
    std::vector<double> x{0.0, 1.0};
    auto ex = E[0](x, 0.0);
    auto ey = E[1](x, 0.0);
    auto ez = E[2](x, 0.0);
    EXPECT_DOUBLE_EQ((*ex)[0], 3.0);
    EXPECT_DOUBLE_EQ((*ey)[0], -6.0);
    EXPECT_DOUBLE_EQ((*ez)[0], 3.0);
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
