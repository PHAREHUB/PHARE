#include "phare_core.hpp"
#include "initializer/data_provider.hpp"

#include "core/data/grid/grid.hpp"
#include "core/data/vecfield/vecfield.hpp"
#include "core/data/vecfield/vecfield_initializer.hpp"
#include "core/numerics/curl/edge_to_face_curl.hpp"
#include "core/data/field/initializers/field_user_initializer.hpp"
#include "core/utilities/span.hpp"

#include "gtest/gtest.h"

#include <array>
#include <cmath>
#include <memory>
#include <string>
#include <vector>
#include <algorithm>

using namespace PHARE::core;
using namespace PHARE::initializer;


// MHD layouts need consistent MHD options, hybrid ones only an interpolation order
constexpr PHARE::SimOpts mhdOpts(std::size_t dim)
{
    return PHARE::SimOpts{dim,
                          0,
                          0,
                          PHARE::MHDOpts::ReconstructionType::WENOZ,
                          PHARE::MHDOpts::SlopeLimiterType::None,
                          PHARE::MHDOpts::RiemannSolverType::Rusanov};
}

template<std::size_t dim_, typename Layout_, typename Quantity_>
struct Case
{
    static constexpr std::size_t dim = dim_;
    using Layout                     = Layout_;
    using Quantity                   = Quantity_;
};

template<std::size_t dim>
using HybridCase
    = Case<dim, typename PHARE_Types<PHARE::SimOpts{dim, 1}>::Hybrid::GridLayout_t, HybridQuantity>;

template<std::size_t dim>
using MHDCase = Case<dim, typename PHARE_Types<mhdOpts(dim)>::MHD::GridLayout_t, MHDQuantity>;


// analytic, non periodic, with a linear part (as a uniform-field vector potential would have)
double ax_fn(double x, double y, double z)
{
    return 0.7 * std::cos(1.3 * y + 0.4 * z) + 0.2 * x * y;
}
double ay_fn(double x, double /*y*/, double z)
{
    return 1.1 * std::sin(0.9 * x - 0.5 * z) - 0.3 * z + 0.15 * x;
}
double az_fn(double x, double y, double /*z*/)
{
    return 2.0 * std::sin(1.7 * x) * std::cos(0.8 * y) + 0.5 * y + 3.0;
}
double bz_fn(double x, double y, double /*z*/)
{
    return 0.4 + 0.1 * x * y;
}
double bx_fn(double x, double y, double z)
{
    return 1.0 + 0.3 * std::sin(x + 2 * y + z);
}
double by_fn(double x, double y, double z)
{
    return -0.2 * std::cos(x * y) + 0.1 * z;
}


template<std::size_t dim, typename F>
InitFunction<dim> make_fn(F f)
{
    using Param  = std::vector<double> const&;
    using Return = std::shared_ptr<Span<double>>;

    if constexpr (dim == 2)
        return [f](Param x, Param y) -> Return {
            std::vector<double> v(x.size());
            for (std::size_t i = 0; i < x.size(); ++i)
                v[i] = f(x[i], y[i], 0.);
            return std::make_shared<VectorSpan<double>>(std::move(v));
        };
    else
        return [f](Param x, Param y, Param z) -> Return {
            std::vector<double> v(x.size());
            for (std::size_t i = 0; i < x.size(); ++i)
                v[i] = f(x[i], y[i], z[i]);
            return std::make_shared<VectorSpan<double>>(std::move(v));
        };
}


template<typename CaseT>
struct VecFieldInitializerTest : public ::testing::Test
{
    static constexpr auto dim = CaseT::dim;
    using Layout              = typename CaseT::Layout;
    using Quantity            = typename CaseT::Quantity;
    using Scalar              = typename Quantity::Scalar;
    using Grid_t              = Grid<NdArrayVector<dim, double>, Scalar>;
    using Field_t             = Field<dim, Scalar>;
    using VecField_t          = VecField<Field_t, Quantity>;

    static Layout make_layout()
    {
        // non-zero AMR offset so AMR <-> local index conversions are exercised
        if constexpr (dim == 2)
            return {{{0.1, 0.2}}, {{20, 16}}, {0.5, -0.3}, Box<int, 2>{{3, 5}, {22, 20}}};
        else
            return {{{0.1, 0.2, 0.3}},
                    {{12, 10, 8}},
                    {0.5, -0.3, 0.2},
                    Box<int, 3>{{3, 5, 1}, {14, 14, 8}}};
    }

    Layout layout = make_layout();

    std::array<Grid_t, 3> bgrids{Grid_t{"B_x", layout, Scalar::Bx, 0.},
                                 Grid_t{"B_y", layout, Scalar::By, 0.},
                                 Grid_t{"B_z", layout, Scalar::Bz, 0.}};
    VecField_t B{"B", Quantity::Vector::B};

    // A sampled independently of the code under test, for the reference stencil
    std::array<Grid_t, 3> agrids{Grid_t{"Ax", layout, Scalar::Ex, 0.},
                                 Grid_t{"Ay", layout, Scalar::Ey, 0.},
                                 Grid_t{"Az", layout, Scalar::Ez, 0.}};

    VecFieldInitializerTest()
    {
        for (std::size_t i = 0; i < 3; ++i)
            B[i].setBuffer(&bgrids[i]);
    }

    void sample_A(bool with_in_plane = true)
    {
        auto fx = with_in_plane ? make_fn<dim>(ax_fn) : make_fn<dim>([](auto...) { return 0.; });
        auto fy = with_in_plane ? make_fn<dim>(ay_fn) : make_fn<dim>([](auto...) { return 0.; });
        FieldUserFunctionInitializer::initialize(*&agrids[0], layout, fx);
        FieldUserFunctionInitializer::initialize(*&agrids[1], layout, fy);
        FieldUserFunctionInitializer::initialize(*&agrids[2], layout, make_fn<dim>(az_fn));
    }

    double A_scale() const
    {
        double m = 0;
        for (auto const& g : agrids)
            for (std::size_t i = 0; i < g.size(); ++i)
                m = std::max(m, std::abs(g.data()[i]));
        return m;
    }

    double min_dl() const
    {
        auto const dl = layout.meshSize();
        return *std::min_element(dl.begin(), dl.end());
    }

    auto range(Scalar qty, Direction dir) const
    {
        auto const [s, e] = layout.ghostStartToEnd(qty, dir);
        return std::array<std::uint32_t, 2>{s, e};
    }

    // hand-written edge -> face two-point curl, independent of GridLayout::deriv
    double ref_B(std::size_t c, std::uint32_t i, std::uint32_t j, std::uint32_t k = 0) const
    {
        auto const& Ax  = *&agrids[0];
        auto const& Ay  = *&agrids[1];
        auto const& Az  = *&agrids[2];
        auto const idx  = layout.inverseMeshSize();
        auto const idx_ = idx[0], idy = idx[1];
        if constexpr (dim == 2)
        {
            if (c == 0)
                return idy * (Az(i, j + 1) - Az(i, j));
            if (c == 1)
                return -(idx_ * (Az(i + 1, j) - Az(i, j)));
            return idx_ * (Ay(i + 1, j) - Ay(i, j)) - idy * (Ax(i, j + 1) - Ax(i, j));
        }
        else
        {
            auto const idz = idx[2];
            if (c == 0)
                return idy * (Az(i, j + 1, k) - Az(i, j, k))
                       - idz * (Ay(i, j, k + 1) - Ay(i, j, k));
            if (c == 1)
                return idz * (Ax(i, j, k + 1) - Ax(i, j, k))
                       - idx_ * (Az(i + 1, j, k) - Az(i, j, k));
            return idx_ * (Ay(i + 1, j, k) - Ay(i, j, k)) - idy * (Ax(i, j + 1, k) - Ax(i, j, k));
        }
    }

    // B component c equals the reference curl at every ghost-box node, returns nbr nodes checked
    std::size_t expect_B_is_curl(std::size_t c)
    {
        auto const qty   = std::array{Scalar::Bx, Scalar::By, Scalar::Bz}[c];
        auto const& Bc   = *&bgrids[c];
        auto const rx    = range(qty, Direction::X);
        auto const ry    = range(qty, Direction::Y);
        std::size_t n    = 0;
        auto const check = [&](double actual, double expected, auto... ijk) {
            ++n;
            if (actual != expected)
            {
                ADD_FAILURE() << "B" << c << " mismatch at (" << ((std::to_string(ijk) + ",") + ...)
                              << ") " << actual << " != " << expected;
            }
        };
        for (auto i = rx[0]; i <= rx[1]; ++i)
            for (auto j = ry[0]; j <= ry[1]; ++j)
            {
                if constexpr (dim == 2)
                {
                    check(Bc(i, j), ref_B(c, i, j), i, j);
                }
                else
                {
                    auto const rz = range(qty, Direction::Z);
                    for (auto k = rz[0]; k <= rz[1]; ++k)
                    {
                        check(Bc(i, j, k), ref_B(c, i, j, k), i, j, k);
                    }
                }
            }
        return n;
    }

    // max |div B| over every cell whose faces all lie in the ghost box
    double max_divB() const
    {
        auto const& Bx  = *&bgrids[0];
        auto const& By  = *&bgrids[1];
        auto const& Bz  = *&bgrids[2];
        auto const idx  = layout.inverseMeshSize();
        auto const cx   = range(Scalar::Bz, Direction::X); // dual in x
        auto const cy   = range(Scalar::Bz, Direction::Y); // dual in y
        double max_divB = 0;
        for (auto i = cx[0]; i <= cx[1]; ++i)
            for (auto j = cy[0]; j <= cy[1]; ++j)
            {
                if constexpr (dim == 2)
                {
                    auto const d
                        = idx[0] * (Bx(i + 1, j) - Bx(i, j)) + idx[1] * (By(i, j + 1) - By(i, j));
                    max_divB = std::max(max_divB, std::abs(d));
                }
                else
                {
                    auto const cz = range(Scalar::Bx, Direction::Z); // dual in z
                    for (auto k = cz[0]; k <= cz[1]; ++k)
                    {
                        auto const d = idx[0] * (Bx(i + 1, j, k) - Bx(i, j, k))
                                       + idx[1] * (By(i, j + 1, k) - By(i, j, k))
                                       + idx[2] * (Bz(i, j, k + 1) - Bz(i, j, k));
                        max_divB = std::max(max_divB, std::abs(d));
                    }
                }
            }
        return max_divB;
    }

    double divB_tolerance() const { return 1e-12 * A_scale() / (min_dl() * min_dl()); }

    // B component c equals the user function sampled at its own centring
    void expect_B_is_sampled(std::size_t c, InitFunction<dim> const& fn)
    {
        auto const qty = std::array{Scalar::Bx, Scalar::By, Scalar::Bz}[c];
        Grid_t expected{"expected", layout, qty, 0.};
        FieldUserFunctionInitializer::initialize(*&expected, layout, fn);
        auto const& actual = bgrids[c];
        ASSERT_EQ(actual.size(), expected.size());
        for (std::size_t i = 0; i < actual.size(); ++i)
        {
            if (actual.data()[i] != expected.data()[i])
            {
                ADD_FAILURE() << "B" << c << " flat index " << i << ": " << actual.data()[i]
                              << " != " << expected.data()[i];
            }
        }
    }

    PHAREDict vecpot_dict(bool with_in_plane = true) const
    {
        PHAREDict dict;
        auto zero                               = make_fn<dim>([](auto...) { return 0.; });
        dict["vector_potential"]["x_component"] = with_in_plane ? make_fn<dim>(ax_fn) : zero;
        dict["vector_potential"]["y_component"] = with_in_plane ? make_fn<dim>(ay_fn) : zero;
        dict["vector_potential"]["z_component"] = make_fn<dim>(az_fn);
        return dict;
    }
};


using Cases = ::testing::Types<HybridCase<2>, HybridCase<3>, MHDCase<2>, MHDCase<3>>;
TYPED_TEST_SUITE(VecFieldInitializerTest, Cases);


TYPED_TEST(VecFieldInitializerTest, curlOfEdgeFieldMatchesStencilOverGhostBox)
{
    this->sample_A();
    curl_edges_to_faces(
        this->layout, std::array{*&this->agrids[0], *&this->agrids[1], *&this->agrids[2]}, this->B);

    for (std::size_t c = 0; c < 3; ++c)
    {
        EXPECT_GT(this->expect_B_is_curl(c), 0u);
    }
}

TYPED_TEST(VecFieldInitializerTest, curlOfEdgeFieldIsDivergenceFreeOverGhostBox)
{
    this->sample_A();
    curl_edges_to_faces(
        this->layout, std::array{*&this->agrids[0], *&this->agrids[1], *&this->agrids[2]}, this->B);

    EXPECT_GT(this->A_scale(), 1.); // non-trivial A, so the bound is meaningful
    EXPECT_LE(this->max_divB(), this->divB_tolerance());
}

TYPED_TEST(VecFieldInitializerTest, vectorPotentialModeSetsBToCurlA)
{
    VecFieldInitializer<TestFixture::dim> init{this->vecpot_dict()};
    init.initialize(this->B, this->layout);

    this->sample_A();
    for (std::size_t c = 0; c < 3; ++c)
    {
        EXPECT_GT(this->expect_B_is_curl(c), 0u);
    }
    EXPECT_LE(this->max_divB(), this->divB_tolerance());
}

TYPED_TEST(VecFieldInitializerTest, directComponentOverridesCurlInVectorPotentialMode)
{
    constexpr auto dim = TestFixture::dim;
    if constexpr (dim == 2)
    {
        auto dict           = this->vecpot_dict(/*with_in_plane=*/false);
        dict["z_component"] = make_fn<dim>(bz_fn);

        VecFieldInitializer<dim> init{dict};
        init.initialize(this->B, this->layout);

        this->sample_A(/*with_in_plane=*/false);
        EXPECT_GT(this->expect_B_is_curl(0), 0u);
        EXPECT_GT(this->expect_B_is_curl(1), 0u);
        this->expect_B_is_sampled(2, make_fn<dim>(bz_fn));

        // Bz does not enter the 2D divergence
        EXPECT_LE(this->max_divB(), this->divB_tolerance());
    }
    else
    {
        GTEST_SKIP() << "out-of-plane override only exists in 2D";
    }
}

TYPED_TEST(VecFieldInitializerTest, componentModeStillSamplesUserFunctions)
{
    constexpr auto dim = TestFixture::dim;
    PHAREDict dict;
    dict["x_component"] = make_fn<dim>(bx_fn);
    dict["y_component"] = make_fn<dim>(by_fn);
    dict["z_component"] = make_fn<dim>(bz_fn);

    VecFieldInitializer<dim> init{dict};
    init.initialize(this->B, this->layout);

    this->expect_B_is_sampled(0, make_fn<dim>(bx_fn));
    this->expect_B_is_sampled(1, make_fn<dim>(by_fn));
    this->expect_B_is_sampled(2, make_fn<dim>(bz_fn));
}

TYPED_TEST(VecFieldInitializerTest, componentModeRequiresAllComponents)
{
    constexpr auto dim = TestFixture::dim;
    PHAREDict dict;
    dict["x_component"] = make_fn<dim>(bx_fn);
    dict["y_component"] = make_fn<dim>(by_fn);

    EXPECT_THROW(VecFieldInitializer<dim>{dict}, std::runtime_error);
}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
