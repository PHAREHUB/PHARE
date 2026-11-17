#include "phare_core.hpp"
#include "simulator/simulator_def.hpp"

#include "core/numerics/faraday/faraday.hpp"
#include "core/numerics/marder/marder.hpp"

#include "tests/core/data/vecfield/test_vecfield_fixtures.hpp"

#include "gtest/gtest.h"

#include <array>
#include <cmath>
#include <cstring>
#include <numbers>
#include <random>
#include <vector>


using namespace PHARE::core;


template<std::size_t dim_, std::size_t interp_order_>
struct Case
{
    static constexpr std::size_t dim          = dim_;
    static constexpr std::size_t interp_order = interp_order_;
};


// even numbers of cells, for the Nyquist mode. Two grids: mesh sizes that differ in every
// direction, and powers of two, on which the arithmetic of the exact checks has no rounding.
template<typename Case_t>
class MarderTest : public ::testing::Test
{
protected:
    static constexpr std::size_t dim          = Case_t::dim;
    static constexpr std::size_t interp_order = Case_t::interp_order;
    using GridLayout_t = PHARE_Types<PHARE::SimOpts{dim, interp_order}>::Hybrid::GridLayout_t;
    using Grid_t       = Grid<NdArrayVector<dim>, HybridQuantity::Scalar>;

    static constexpr std::array<std::uint32_t, 3> allCells{16, 12, 10};
    static constexpr std::array<double, 3> generalMeshSizes{0.1, 0.13, 0.07};
    static constexpr std::array<double, 3> dyadicMeshSizes{0.5, 0.25, 0.125};
    static constexpr double c2 = 4.;

    static auto sub(auto const& all)
    {
        std::array<typename std::decay_t<decltype(all)>::value_type, dim> a;
        std::copy_n(all.begin(), dim, a.begin());
        return a;
    }

    // the electric field, the charge density and B (for Faraday) on one layout, whose AMR box
    // starts at `lower` (0 by default)
    struct State
    {
        explicit State(std::array<double, 3> const& meshSizes,
                       std::array<int, dim> const& lower = {})
            : layout{sub(meshSizes), sub(allCells), Point<double, dim>{}, amrBox(lower)}
        {
        }

        static Box<int, dim> amrBox(std::array<int, dim> const& lower)
        {
            auto upper = lower;
            for (std::size_t d = 0; d < dim; ++d)
                upper[d] += static_cast<int>(allCells[d]) - 1;
            return {Point<int, dim>{lower}, Point<int, dim>{upper}};
        }

        GridLayout_t layout;
        UsableVecField<dim> E{"E", layout, HybridQuantity::Vector::E};
        UsableVecField<dim> B{"B", layout, HybridQuantity::Vector::B};
        Grid_t rho{"rho", layout, HybridQuantity::Scalar::rho, 0.};

        // d such that d dt sum_i 2/dl_i^2 = s
        double coefficient(double const s, double const dt) const
        {
            return marderCoefficient(layout, s, dt);
        }
    };

    State general{generalMeshSizes};
    State dyadic{dyadicMeshSizes};
    static constexpr double dt = 0.01;


    static bool isPhysical(auto const& layout, auto const& f, auto const&... ijk)
    {
        std::array<std::uint32_t, dim> const index{static_cast<std::uint32_t>(ijk)...};
        for (std::size_t d = 0; d < dim; ++d)
        {
            auto const direction = static_cast<Direction>(d);
            if (index[d] < layout.physicalStartIndex(f, direction)
                or index[d] > layout.physicalEndIndex(f, direction))
                return false;
        }
        return true;
    }

    // ghosts of a periodic domain made of this one patch: copies of the physical nodes one
    // period away (primal and dual nodes alike start at the same index)
    static void fillPeriodicGhosts(auto const& layout, auto& f)
    {
        layout.evalOnGhostBox(f, [&](auto const&... ijk) {
            if (isPhysical(layout, f, ijk...))
                return;
            std::array<std::uint32_t, dim> source{static_cast<std::uint32_t>(ijk)...};
            for (std::size_t d = 0; d < dim; ++d)
            {
                auto const direction = static_cast<Direction>(d);
                if (source[d] < layout.physicalStartIndex(f, direction))
                    source[d] += allCells[d];
                else if (source[d] > layout.physicalEndIndex(f, direction))
                    source[d] -= allCells[d];
            }
            f(ijk...) = std::apply([&](auto... s) { return f(s...); }, source);
        });
    }

    static auto position(auto const& layout, auto const& f, auto const&... ijk)
    {
        return layout.fieldNodeCoordinates(
            f, layout.localToAMR(Point{static_cast<std::uint32_t>(ijk)...}.as_signed()));
    }

    // div E - c2 rho on the physical nodes of rho, written independently of the operator
    static auto gaussResidual(State const& s)
    {
        std::vector<double> F;
        s.layout.evalOnBox(*s.rho, [&](auto const&... ijk) {
            std::array<std::uint32_t, dim> const index{static_cast<std::uint32_t>(ijk)...};
            double div = 0.;
            for (std::size_t d = 0; d < dim; ++d)
            {
                auto left = index;
                --left[d]; // the dual node at i - 1/2 has index i - 1
                auto const& Ed = s.E[d];
                div += (std::apply([&](auto... i) { return Ed(i...); }, index)
                        - std::apply([&](auto... i) { return Ed(i...); }, left))
                       / s.layout.meshSize()[d];
            }
            F.push_back(div - c2 * (*s.rho)(ijk...));
        });
        return F;
    }

    static double laplacianEigenvalue(State const& s, std::array<int, dim> const& modes)
    {
        double lambda = 0.;
        for (std::size_t d = 0; d < dim; ++d)
        {
            auto const dl     = s.layout.meshSize()[d];
            auto const sinkdl = std::sin(std::numbers::pi * modes[d] / allCells[d]);
            lambda += 4. / (dl * dl) * sinkdl * sinkdl;
        }
        return lambda;
    }

    // rho = -prod_d cos(k_d x_d + phase_d) / c2 and E = 0, so that F = prod_d cos(...)
    static void seed(State& s, std::array<int, dim> const& modes)
    {
        for (auto& c : s.E)
            std::fill(c.data(), c.data() + c.size(), 0.);
        auto& rho = *s.rho;
        s.layout.evalOnGhostBox(rho, [&](auto const&... ijk) {
            auto const r = position(s.layout, rho, ijk...);
            double value = 1.;
            for (std::size_t d = 0; d < dim; ++d)
            {
                auto const length = allCells[d] * s.layout.meshSize()[d];
                value *= std::cos(2 * std::numbers::pi * modes[d] / length * r[d] + 0.3 + d);
            }
            rho(ijk...) = -value / c2;
        });
    }

    // apply the operator n times on a periodic domain, F after each application must be
    // F0 times the predicted amplification factor to the power of the step
    static void checkAmplification(State& s, Marder<GridLayout_t> const& marder, double const g,
                                   int const nsteps)
    {
        auto const F0 = gaussResidual(s);
        for (int step = 1; step <= nsteps; ++step)
        {
            marder(s.E, *s.rho, dt);
            for (auto& c : s.E)
                fillPeriodicGhosts(s.layout, c);
            auto const F      = gaussResidual(s);
            auto const factor = std::pow(g, step);
            for (std::size_t n = 0; n < F.size(); ++n)
                ASSERT_NEAR(F[n], factor * F0[n], 1e-12 * std::max(1., std::abs(factor)))
                    << "step " << step << " node " << n << " factor " << factor;
        }
    }
};

using Cases = ::testing::Types<Case<1, 1>, Case<1, 2>, Case<1, 3>, Case<2, 1>, Case<2, 2>,
                               Case<2, 3>, Case<3, 1>, Case<3, 2>, Case<3, 3>>;
TYPED_TEST_SUITE(MarderTest, Cases);



// (a) F = div E - c2 rho is diffused by the discrete Laplacian: a mode is multiplied each step
// by g = 1 - d dt lambda_k, lambda_k = sum_i 4/dl_i^2 sin^2(k_i dl_i / 2)
TYPED_TEST(MarderTest, seededResidualDecaysAtTheDiscreteLaplacianRate)
{
    constexpr auto dim = TestFixture::dim;
    auto& s            = this->general;
    double const d     = s.coefficient(0.4, this->dt);
    Marder<typename TestFixture::GridLayout_t> const marder{s.layout, this->c2, d};

    std::vector<std::array<int, dim>> allModes;
    for (int const mode : {1, 3})
    {
        std::array<int, dim> modes{};
        modes[0] = mode;
        allModes.push_back(modes);
    }
    std::array<int, dim> mixed;
    for (std::size_t i = 0; i < dim; ++i)
        mixed[i] = static_cast<int>(i) + 1;
    allModes.push_back(mixed);

    for (auto const& modes : allModes)
    {
        this->seed(s, modes);
        auto const g = 1. - d * this->dt * this->laplacianEigenvalue(s, modes);
        this->checkAmplification(s, marder, g, 5);
    }
}


// (b) the correction is a discrete gradient: its discrete curl is zero, Faraday sees nothing.
// Exactly zero on the power of two grid (no rounding), rounding level on the general one.
TYPED_TEST(MarderTest, correctionHasNoCurl)
{
    for (auto* s : {&this->dyadic, &this->general})
    {
        // on the power of two grid, d and dt are powers of two and rho integers, so that every
        // operation is exact
        bool const exact = s == &this->dyadic;
        double const dt  = exact ? 1. / 128 : this->dt;
        double const d   = exact ? 1. / 8 : s->coefficient(0.25, dt);
        Marder<typename TestFixture::GridLayout_t> const marder{s->layout, this->c2, d};

        std::mt19937 gen{12};
        std::uniform_int_distribution<int> integers{-8, 8};
        auto& rho = *s->rho;
        s->layout.evalOnGhostBox(rho, [&](auto const&... ijk) { rho(ijk...) = integers(gen); });

        marder(s->E, rho, dt); // E = 0 before: E is now the correction

        double maxE = 0.;
        for (auto& c : s->E)
            s->layout.evalOnGhostBox(c, [&](auto const&... ijk) {
                if (!TestFixture::isPhysical(s->layout, c, ijk...))
                {
                    EXPECT_EQ(c(ijk...), 0.) << c.name() << " ghost modified";
                }
                maxE = std::max(maxE, std::abs(c(ijk...)));
            });
        ASSERT_GT(maxE, 0.);

        Faraday<typename TestFixture::GridLayout_t>{s->layout}(s->B, s->E, s->B, 1.);
        double maxCurl = 0.;
        for (auto& c : s->B)
            s->layout.evalOnBox(
                c, [&](auto const&... ijk) { maxCurl = std::max(maxCurl, std::abs(c(ijk...))); });

        if (exact)
            EXPECT_EQ(maxCurl, 0.);
        else
        {
            auto const minDl
                = *std::min_element(s->layout.meshSize().begin(), s->layout.meshSize().end());
            EXPECT_LT(maxCurl, 1e-14 * maxE / minDl);
            std::ostringstream os;
            os << std::scientific << maxCurl / (maxE / minDl);
            ::testing::Test::RecordProperty("curl_over_E_dl_general_grid", os.str());
        }
    }
}


// (c) a state with div E = c2 rho exactly is not modified, ghosts included
TYPED_TEST(MarderTest, gaussExactStateIsUnchanged)
{
    constexpr auto dim = TestFixture::dim;
    auto& s            = this->dyadic;
    double const d     = s.coefficient(0.5, this->dt);
    Marder<typename TestFixture::GridLayout_t> const marder{s.layout, this->c2, d};

    // E_c = slope_c i_c + a function of the other indices, multiples of 1/64: div E is
    // sum_c slope_c / dl_c exactly, and so is c2 rho
    std::array<double, 3> const slopes{0.5, -0.25, 0.75};
    std::mt19937 gen{7};
    std::uniform_int_distribution<int> integers{-64, 64};
    double divE = 0.;
    for (std::size_t c = 0; c < 3; ++c)
    {
        auto& Ec = s.E[c];
        std::vector<double> transverse(Ec.size());
        for (auto& v : transverse)
            v = integers(gen) / 64.;
        s.layout.evalOnGhostBox(Ec, [&](auto const&... ijk) {
            std::array<std::uint32_t, dim> index{static_cast<std::uint32_t>(ijk)...};
            double along = 0.;
            if (c < dim)
            {
                along    = slopes[c] * index[c];
                index[c] = 0; // the transverse part does not depend on i_c
            }
            auto const flat = std::apply([&](auto... i) { return &Ec(i...) - Ec.data(); }, index);
            Ec(ijk...)      = along + transverse[flat];
        });
        if (c < dim)
            divE += slopes[c] / s.layout.meshSize()[c];
    }
    auto& rho = *s.rho;
    std::fill(rho.data(), rho.data() + rho.size(), divE / this->c2);

    std::array<std::vector<double>, 3> before;
    for (std::size_t c = 0; c < 3; ++c)
        before[c].assign(s.E[c].data(), s.E[c].data() + s.E[c].size());

    marder(s.E, rho, this->dt);

    for (std::size_t c = 0; c < 3; ++c)
        EXPECT_EQ(std::memcmp(before[c].data(), s.E[c].data(), before[c].size() * sizeof(double)),
                  0)
            << s.E[c].name() << " modified";
}


// (d) stability: the Nyquist mode in every direction has lambda = sum_i 4/dl_i^2, so
// g = 1 - 2 s with s = d dt sum_i 2/dl_i^2: s = 1 is marginal (g = -1), s > 1 grows
TYPED_TEST(MarderTest, stabilityLimit)
{
    constexpr auto dim = TestFixture::dim;
    auto& s            = this->general;

    std::array<int, dim> nyquist;
    for (std::size_t i = 0; i < dim; ++i)
        nyquist[i] = static_cast<int>(this->allCells[i] / 2);

    for (double const stability : {1., 1.05})
    {
        double const d = s.coefficient(stability, this->dt);
        Marder<typename TestFixture::GridLayout_t> const marder{s.layout, this->c2, d};
        this->seed(s, nyquist);
        auto const g = 1. - d * this->dt * this->laplacianEigenvalue(s, nyquist);
        EXPECT_NEAR(g, 1. - 2. * stability, 1e-14);
        this->checkAmplification(s, marder, g, 40);
    }
}



// (e) skipped nodes, as the solver skips the level boundary nodes: F = 0 there. A residual on
// them alone does not move E; with any state, the E nodes away from them are bitwise those of
// the cleaning without skipped nodes, and the two E nodes next to each skipped node (one per
// direction and side) differ. Skipped: both faces normal to x, and one interior node. The AMR box
// does not start at 0, as on a patch: the skipped boxes are AMR indexes, not local ones.
TYPED_TEST(MarderTest, skippedNodesHaveNoResidual)
{
    constexpr auto dim = TestFixture::dim;
    using Box_t        = Box<int, dim>;
    typename TestFixture::State s{TestFixture::dyadicMeshSizes, TestFixture::sub(std::array{5, -7, 3})};
    auto const& layout = s.layout;
    ASSERT_EQ(layout.AMRBox().lower[0], 5);
    double const dt    = 1. / 128;
    Marder<typename TestFixture::GridLayout_t> const marder{layout, this->c2, 1. / 8};

    // AMR indexes of the primal nodes: node i is the lower corner of cell i
    auto const cells = layout.AMRBox();
    auto lowFace     = Box_t{cells.lower, cells.upper};
    for (std::size_t d = 0; d < dim; ++d)
        lowFace.upper[d] += 1;
    auto highFace     = lowFace;
    lowFace.upper[0]  = lowFace.lower[0];
    highFace.lower[0] = highFace.upper[0];
    auto interior     = Box_t{cells.lower, cells.lower};
    for (std::size_t d = 0; d < dim; ++d)
        interior.lower[d] = interior.upper[d] = cells.lower[d] + 3;
    std::vector<Box_t> const skipped{lowFace, highFace, interior};

    auto const isSkipped = [&](std::array<std::uint32_t, dim> const& local) {
        auto const amr = layout.localToAMR(Point<std::uint32_t, dim>{local});
        return isIn(amr, skipped);
    };
    // E_d at local index i (dual node i + 1/2 along d) reads F on the nodes i and i + 1 along d
    auto const nextToSkipped = [&](std::size_t const d, auto const&... ijk) {
        std::array<std::uint32_t, dim> left{static_cast<std::uint32_t>(ijk)...};
        auto right = left;
        ++right[d];
        return isSkipped(left) or isSkipped(right);
    };

    auto snapshot = [](auto const& E) {
        std::array<std::vector<double>, 3> c;
        for (std::size_t i = 0; i < 3; ++i)
            c[i].assign(E[i].data(), E[i].data() + E[i].size());
        return c;
    };
    auto restore = [](auto& E, auto const& c) {
        for (std::size_t i = 0; i < 3; ++i)
            std::copy(c[i].begin(), c[i].end(), E[i].data());
    };

    auto& rho = *s.rho;
    std::mt19937 gen{3};
    std::uniform_int_distribution<int> integers{-64, 64};

    // a residual on the skipped nodes only: E = 0, rho on the skipped nodes
    for (auto& c : s.E)
        std::fill(c.data(), c.data() + c.size(), 0.);
    std::fill(rho.data(), rho.data() + rho.size(), 0.);
    layout.evalOnBox(rho, [&](auto const&... ijk) {
        if (isSkipped({static_cast<std::uint32_t>(ijk)...}))
            rho(ijk...) = integers(gen) | 1; // odd: never 0
    });
    auto const zero = snapshot(s.E);
    marder(s.E, rho, dt);
    EXPECT_NE(snapshot(s.E), zero) << "the residual must move E without skipped nodes";
    restore(s.E, zero);
    marder(s.E, rho, dt, skipped);
    EXPECT_EQ(snapshot(s.E), zero) << "a residual on skipped nodes moved E";

    // any state, multiples of 1/64 (exact arithmetic)
    for (auto& c : s.E)
        layout.evalOnGhostBox(c, [&](auto const&... ijk) { c(ijk...) = integers(gen) / 64.; });
    layout.evalOnGhostBox(rho, [&](auto const&... ijk) { rho(ijk...) = integers(gen) / 64.; });
    auto const before = snapshot(s.E);
    marder(s.E, rho, dt);
    auto const all = snapshot(s.E);
    restore(s.E, before);
    marder(s.E, rho, dt, skipped);

    std::size_t near = 0, differ = 0;
    for (std::size_t d = 0; d < dim; ++d)
    {
        auto const& Ed = s.E[d];
        layout.evalOnBox(Ed, [&](auto const&... ijk) {
            auto const flat = &Ed(ijk...) - Ed.data();
            if (nextToSkipped(d, ijk...))
            {
                ++near;
                differ += Ed.data()[flat] != all[d][flat];
            }
            else
                EXPECT_EQ(Ed.data()[flat], all[d][flat]) << Ed.name() << " node " << flat;
        });
    }
    ASSERT_GT(near, 0u);
    EXPECT_GT(differ, 0u);
    for (std::size_t c = dim; c < 3; ++c)
        EXPECT_EQ(std::vector<double>(s.E[c].data(), s.E[c].data() + s.E[c].size()), before[c]);
}



int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
