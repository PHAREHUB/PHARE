

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <list>
#include <random>

#include "core/utilities/box/box.hpp"
#include "core/utilities/range/range.hpp"
#include "phare_core.hpp"
#include "core/data/electromag/electromag.hpp"
#include "core/data/field/field.hpp"
#include "core/data/grid/gridlayout.hpp"
#include "core/data/grid/gridlayout_impl.hpp"
#include "core/data/ndarray/ndarray_vector.hpp"
#include "core/data/particles/particle.hpp"
#include "core/data/particles/particle_array.hpp"
#include "core/data/vecfield/vecfield.hpp"
#include "core/hybrid/hybrid_quantities.hpp"
#include "core/numerics/interpolator/interpolator.hpp"

#include "tests/core/data/vecfield/test_vecfield_fixtures.hpp"

using namespace PHARE::core;



template<typename Weighter>
class AWeighter : public ::testing::Test
{
public:
    // All tests using this class are 1D and should remain that way or update this.
    auto static constexpr dimension    = 1;
    auto static constexpr interp_order = Weighter::interp_order;

    using GridLayout_t   = GridLayout<GridLayoutImplYee<dimension, interp_order>>;
    using Interpolator_t = Interpolator<dimension, interp_order>;


    AWeighter()
    {
        std::random_device rd;  // Will be used to obtain a seed for the random number engine
        std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
        std::uniform_real_distribution<> dis(5, 10);


        // generate the nbr_tests random (normalized) particle position
        std::generate(std::begin(normalizedPositions), std::end(normalizedPositions),
                      [&dis, &gen]() { return dis(gen); });


        // now for each random position, calculate
        // the start index and the N(interp_order) weights
        for (auto i = 0u; i < nbr_tests; ++i)
        {
            auto icell = static_cast<int>(normalizedPositions[i]);
            auto delta = normalizedPositions[i] - icell;
            auto startIndex
                = icell
                  - Interpolator_t::template computeStartLeftShift<QtyCentering,
                                                                   QtyCentering::primal>(delta);
            this->weighter.computeWeight(normalizedPositions[i], startIndex, weights_[i]);
        }


        // for each particle, we have N(interp_order) weights
        // the sum of these N(interp_order) weights for each particle is stored in weightsSum
        std::transform(std::begin(weights_), std::end(weights_), std::begin(weightsSums),
                       [](auto const& weight_list) {
                           return std::accumulate(std::begin(weight_list), std::end(weight_list),
                                                  0.);
                       });
    }

protected:
    Weighter weighter;
    static int const nbr_tests = 100000;
    std::array<double, nbr_tests> normalizedPositions;
    std::array<double, nbr_tests> weightsSums;
    std::array<std::array<double, nbrPointsSupport(Weighter::interp_order)>, nbr_tests> weights_;
};


using Weighters = ::testing::Types<Weighter<1>, Weighter<2>, Weighter<3>>;




TYPED_TEST_SUITE(AWeighter, Weighters);


TYPED_TEST(AWeighter, ComputesWeightThatSumIsOne)
{
    auto equalsOne = [](double sum) { return std::abs(sum - 1.) < 1e-10; };
    EXPECT_TRUE(std::all_of(std::begin(this->weightsSums), std::end(this->weightsSums), equalsOne));
}




TEST(Weights, NbrPointsInBSplineSupportIsCorrect)
{
    EXPECT_EQ(2, nbrPointsSupport(1));
    EXPECT_EQ(3, nbrPointsSupport(2));
    EXPECT_EQ(4, nbrPointsSupport(3));
}




struct BSpline
{
    std::vector<std::vector<double>> weights;
    std::vector<std::vector<int>> indexes;
};


// this function reads the indices and weights
// from the file generated by the python script interpolator_test.py
BSpline readFromFile(std::string centering, std::size_t order)
{
    std::ostringstream oss;
    oss << "bsplines_" << order << "_" << centering << ".dat";

    std::ifstream file{oss.str()};
    if (!file)
    {
        std::cout << "file not found\n";
    }

    BSpline bs;
    bs.indexes.resize(10);
    bs.weights.resize(10);

    for (auto ipos = 0u; ipos < 10; ++ipos)
    {
        bs.indexes[ipos].resize(order + 1);
        bs.weights[ipos].resize(order + 1);

        file.read(reinterpret_cast<char*>(bs.indexes[ipos].data()),
                  static_cast<std::streamsize>(bs.indexes[ipos].size() * sizeof(int)));

        file.read(reinterpret_cast<char*>(bs.weights[ipos].data()),
                  static_cast<std::streamsize>(bs.weights[ipos].size() * sizeof(double)));
    }
    return bs;
}

template<typename AWeighter_t, typename Centering, Centering centering, typename Weighter>
void check_bspline(Weighter& weighter, std::string centering_id)
{
    using Interpolator_t        = typename AWeighter_t::Interpolator_t;
    constexpr auto interp_order = AWeighter_t::interp_order;

    auto data = readFromFile(centering_id, interp_order);
    std::array<double, nbrPointsSupport(interp_order)> weights;

    // python file hard-codes 10 particle positions
    // that start at x = 3 every 0.1

    std::size_t const n_particles = 10;
    double const icell            = 3;
    double const dx               = 0.1;

    for (auto ipos = 0u; ipos < n_particles; ++ipos)
    {
        auto delta = static_cast<double>(ipos) * dx;

        auto startIndex
            = icell - Interpolator_t::template computeStartLeftShift<Centering, centering>(delta);

        double normalizedPosition = icell + delta;
        if constexpr (centering == QtyCentering::dual)
            normalizedPosition -= .5;

        weighter.computeWeight(normalizedPosition, startIndex, weights);

        for (auto inode = 0u; inode < weights.size(); ++inode)
        {
            std::cout << inode << " " << ipos << " " << data.weights[ipos][inode] << " =? "
                      << weights[inode] << "\n";
            EXPECT_DOUBLE_EQ(data.weights[ipos][inode], weights[inode]);
        }
    }
}


TYPED_TEST(AWeighter, computesPrimalBSplineWeightsForAnyParticlePosition)
{
    using AWeighter_t    = TestFixture;
    using Interpolator_t = typename AWeighter_t::Interpolator_t;
    using GridLayout_t   = typename AWeighter_t::GridLayout_t;

    static_assert(Interpolator_t::interp_order == GridLayout_t::interp_order);

    check_bspline<AWeighter_t, QtyCentering, QtyCentering::primal>(this->weighter, "primal");
}
TYPED_TEST(AWeighter, computesDualBSplineWeightsForAnyParticlePosition)
{
    using AWeighter_t    = TestFixture;
    using Interpolator_t = typename AWeighter_t::Interpolator_t;
    using GridLayout_t   = typename AWeighter_t::GridLayout_t;

    static_assert(Interpolator_t::interp_order == GridLayout_t::interp_order);

    check_bspline<AWeighter_t, QtyCentering, QtyCentering::dual>(this->weighter, "dual");
}


template<typename InterpolatorT>
class A1DInterpolator : public ::testing::Test
{
public:
    static constexpr auto dimension    = InterpolatorT::dimension;
    static constexpr auto interp_order = InterpolatorT::interp_order;
    // arbitrary number of cells
    static constexpr std::uint32_t nx = 50;

    using PHARE_TYPES      = PHARE::core::PHARE_Types<dimension, interp_order>;
    using GridLayout_t     = typename PHARE_TYPES::GridLayout_t;
    using ParticleArray_t  = typename PHARE_TYPES::ParticleArray_t;
    using Electromag_t     = typename PHARE_TYPES::Electromag_t;
    using UsableVecFieldND = UsableVecField<dimension>;

    Electromag_t em;
    GridLayout_t layout{{0.1}, {nx}, {0.}};
    ParticleArray_t particles;
    InterpolatorT interp;
    constexpr static auto safeLayer = static_cast<int>(1 + ghostWidthForParticles<interp_order>());

    UsableVecFieldND B, E;

    static constexpr double ex0 = 2.25;
    static constexpr double ey0 = 2.50;
    static constexpr double ez0 = 2.75;
    static constexpr double bx0 = 2.25;
    static constexpr double by0 = 2.50;
    static constexpr double bz0 = 2.75;

    A1DInterpolator()
        : em{"EM"}
        , particles{grow(layout.AMRBox(), safeLayer), 1}
        , B{"EM_B", layout, HybridQuantity::Vector::B}
        , E{"EM_E", layout, HybridQuantity::Vector::E}
    {
        for (auto ix = 0u; ix < nx; ++ix) // B & E are constant on their grid
        {
            B(Component::X)(ix) = bx0;
            B(Component::Y)(ix) = by0;
            B(Component::Z)(ix) = bz0;

            E(Component::X)(ix) = ex0;
            E(Component::Y)(ix) = ey0;
            E(Component::Z)(ix) = ez0;
        }

        for (auto& part : particles)
        {
            part.iCell[0] = 5;
            part.delta[0] = 0.32;
        }

        B.set_on(em.B);
        E.set_on(em.E);
    }
};




using Interpolators1D
    = ::testing::Types<Interpolator<1, 1>, Interpolator<1, 2>, Interpolator<1, 3>>;

TYPED_TEST_SUITE(A1DInterpolator, Interpolators1D);



TYPED_TEST(A1DInterpolator, canComputeAllEMfieldsAtParticle)
{
    for (auto const& part : this->particles)
    {
        auto const [E, B]        = this->interp(part, this->em, this->layout);
        auto const& [Ex, Ey, Ez] = E;
        auto const& [Bx, By, Bz] = B;
        EXPECT_NEAR(Ex, this->ex0, 1e-8);
        EXPECT_NEAR(Ey, this->ey0, 1e-8);
        EXPECT_NEAR(Ez, this->ez0, 1e-8);
        EXPECT_NEAR(Bx, this->bx0, 1e-8);
        EXPECT_NEAR(By, this->by0, 1e-8);
        EXPECT_NEAR(Bz, this->bz0, 1e-8);
    }
}



template<typename InterpolatorT>
class A2DInterpolator : public ::testing::Test
{
public:
    static constexpr auto dimension    = InterpolatorT::dimension;
    static constexpr auto interp_order = InterpolatorT::interp_order;
    // arbitrary number of cells
    static constexpr std::uint32_t nx = 50;
    static constexpr std::uint32_t ny = 50;

    using PHARE_TYPES      = PHARE::core::PHARE_Types<dimension, interp_order>;
    using GridLayoutImpl   = GridLayoutImplYee<dimension, interp_order>;
    using ParticleArray_t  = typename PHARE_TYPES::ParticleArray_t;
    using Electromag_t     = typename PHARE_TYPES::Electromag_t;
    using UsableVecFieldND = UsableVecField<dimension>;

    Electromag_t em;
    GridLayout<GridLayoutImpl> layout{{0.1, 0.1}, {nx, ny}, {0., 0.}};
    ParticleArray_t particles;
    InterpolatorT interp;
    constexpr static auto safeLayer = static_cast<int>(1 + ghostWidthForParticles<interp_order>());

    UsableVecFieldND B, E;

    static constexpr double ex0 = 2.25;
    static constexpr double ey0 = 2.50;
    static constexpr double ez0 = 2.75;
    static constexpr double bx0 = 2.25;
    static constexpr double by0 = 2.50;
    static constexpr double bz0 = 2.75;

    A2DInterpolator()
        : em{"EM"}
        , particles{grow(layout.AMRBox(), safeLayer), 1}
        , B{"EM_B", layout, HybridQuantity::Vector::B}
        , E{"EM_E", layout, HybridQuantity::Vector::E}
    {
        for (auto ix = 0u; ix < nx; ++ix)
        {
            for (auto iy = 0u; iy < ny; ++iy)
            {
                B(Component::X)(ix, iy) = bx0;
                B(Component::Y)(ix, iy) = by0;
                B(Component::Z)(ix, iy) = bz0;
                E(Component::X)(ix, iy) = ex0;
                E(Component::Y)(ix, iy) = ey0;
                E(Component::Z)(ix, iy) = ez0;
            }
        }

        for (auto& part : particles)
        {
            part.iCell[0] = 5;
            part.delta[0] = 0.32;
        }

        B.set_on(em.B);
        E.set_on(em.E);
    }
};




using Interpolators2D
    = ::testing::Types<Interpolator<2, 1>, Interpolator<2, 2>, Interpolator<2, 3>>;

TYPED_TEST_SUITE(A2DInterpolator, Interpolators2D);



TYPED_TEST(A2DInterpolator, canComputeAllEMfieldsAtParticle)
{
    for (auto const& part : this->particles)
    {
        auto const [E, B]        = this->interp(part, this->em, this->layout);
        auto const& [Ex, Ey, Ez] = E;
        auto const& [Bx, By, Bz] = B;
        EXPECT_NEAR(Ex, this->ex0, 1e-8);
        EXPECT_NEAR(Ey, this->ey0, 1e-8);
        EXPECT_NEAR(Ez, this->ez0, 1e-8);
        EXPECT_NEAR(Bx, this->bx0, 1e-8);
        EXPECT_NEAR(By, this->by0, 1e-8);
        EXPECT_NEAR(Bz, this->bz0, 1e-8);
    }
}




template<typename InterpolatorT>
class A3DInterpolator : public ::testing::Test
{
public:
    static constexpr auto dimension    = InterpolatorT::dimension;
    static constexpr auto interp_order = InterpolatorT::interp_order;
    // arbitrary number of cells
    static constexpr std::uint32_t nx = 50;
    static constexpr std::uint32_t ny = 50;
    static constexpr std::uint32_t nz = 50;

    using PHARE_TYPES      = PHARE::core::PHARE_Types<dimension, interp_order>;
    using GridLayout_t     = typename PHARE_TYPES::GridLayout_t;
    using ParticleArray_t  = typename PHARE_TYPES::ParticleArray_t;
    using Electromag_t     = typename PHARE_TYPES::Electromag_t;
    using UsableVecFieldND = UsableVecField<dimension>;

    Electromag_t em;
    GridLayout_t layout{{0.1, 0.1, 0.1}, {nx, ny, nz}, {0., 0., 0.}};
    ParticleArray_t particles;
    InterpolatorT interp;
    constexpr static auto safeLayer = static_cast<int>(1 + ghostWidthForParticles<interp_order>());

    UsableVecFieldND B, E;

    static constexpr double ex0 = 2.25;
    static constexpr double ey0 = 2.50;
    static constexpr double ez0 = 2.75;
    static constexpr double bx0 = 2.25;
    static constexpr double by0 = 2.50;
    static constexpr double bz0 = 2.75;

    A3DInterpolator()
        : em{"EM"}
        , particles{grow(layout.AMRBox(), safeLayer), 1}
        , B{"EM_B", layout, HybridQuantity::Vector::B}
        , E{"EM_E", layout, HybridQuantity::Vector::E}
    {
        for (auto ix = 0u; ix < nx; ++ix)
        {
            for (auto iy = 0u; iy < ny; ++iy)
            {
                for (auto iz = 0u; iz < nz; ++iz)
                {
                    B(Component::X)(ix, iy, iz) = bx0;
                    B(Component::Y)(ix, iy, iz) = by0;
                    B(Component::Z)(ix, iy, iz) = bz0;
                    E(Component::X)(ix, iy, iz) = ex0;
                    E(Component::Y)(ix, iy, iz) = ey0;
                    E(Component::Z)(ix, iy, iz) = ez0;
                }
            }
        }

        for (auto& part : particles)
        {
            part.iCell[0] = 5;
            part.delta[0] = 0.32;
        }

        B.set_on(em.B);
        E.set_on(em.E);
    }
};




using Interpolators3D
    = ::testing::Types<Interpolator<3, 1>, Interpolator<3, 2>, Interpolator<3, 3>>;

TYPED_TEST_SUITE(A3DInterpolator, Interpolators3D);



TYPED_TEST(A3DInterpolator, canComputeAllEMfieldsAtParticle)
{
    for (auto const& part : this->particles)
    {
        auto const [E, B]        = this->interp(part, this->em, this->layout);
        auto const& [Ex, Ey, Ez] = E;
        auto const& [Bx, By, Bz] = B;
        EXPECT_NEAR(Ex, this->ex0, 1e-8);
        EXPECT_NEAR(Ey, this->ey0, 1e-8);
        EXPECT_NEAR(Ez, this->ez0, 1e-8);
        EXPECT_NEAR(Bx, this->bx0, 1e-8);
        EXPECT_NEAR(By, this->by0, 1e-8);
        EXPECT_NEAR(Bz, this->bz0, 1e-8);
    }
}




// set a collection of particle (the number depending on interpOrder) so that
// their cumulative density equals 1 at index 20. idem for velocity components...




template<typename Interpolator>
class ACollectionOfParticles_1d : public ::testing::Test
{
    static constexpr auto dimension    = Interpolator::dimension;
    static constexpr auto interp_order = Interpolator::interp_order;

    using PHARE_TYPES      = PHARE::core::PHARE_Types<dimension, interp_order>;
    using ParticleArray_t  = typename PHARE_TYPES::ParticleArray_t;
    using GridLayout_t     = typename PHARE_TYPES::GridLayout_t;
    using Grid_t           = typename PHARE_TYPES::Grid_t;
    using Particle_t       = typename ParticleArray_t::Particle_t;
    using UsableVecFieldND = UsableVecField<dimension>;

public:
    static constexpr std::uint32_t nx        = 30;
    static constexpr std::uint32_t nbrPoints = nbrPointsSupport(Interpolator::interp_order);
    static constexpr std::uint32_t numOfPart = Interpolator::interp_order + 2;

    GridLayout_t layout{{0.1}, {nx}, {0.}};
    constexpr static auto safeLayer = static_cast<int>(1 + ghostWidthForParticles<interp_order>());

    Particle_t part;
    ParticleArray_t particles;

    Grid_t rho;
    Grid_t rho_c;

    UsableVecFieldND v;
    std::array<double, nbrPointsSupport(Interpolator::interp_order)> weights;


    // For the set of particles herebelow, whatever the interp_order,
    // we have a set of particles around idx=20, so that their particle
    // density is 1, their charge density is 2 and velocity is (2, -1, 1)
    // depending on the interp order, iCell, delta, weight, charge and v
    // are then correctly tuned

    ACollectionOfParticles_1d()
        : part{}
        , particles{grow(layout.AMRBox(), safeLayer)}
        , rho{"field", HybridQuantity::Scalar::rho, nx}
        , rho_c{"field", HybridQuantity::Scalar::rho, nx}
        , v{"v", layout, HybridQuantity::Vector::V}
    {
        if constexpr (Interpolator::interp_order == 1)
        {
            part.iCell[0] = 19; // AMR index
            part.delta[0] = 0.5;
            part.weight   = 1.0;
            part.charge   = 2.0;
            part.v[0]     = +2.;
            part.v[1]     = -1.;
            part.v[2]     = +1.;
            particles.push_back(part);

            part.iCell[0] = 20; // AMR index
            part.delta[0] = 0.5;
            part.weight   = 0.4;
            part.charge   = 1.85;
            part.v[0]     = +2.;
            part.v[1]     = -1.;
            part.v[2]     = +1.;
            particles.push_back(part);

            part.iCell[0] = 20; // AMR index
            part.delta[0] = 0.5;
            part.weight   = 0.6;
            part.charge   = 2.1;
            part.v[0]     = +2.;
            part.v[1]     = -1.;
            part.v[2]     = +1.;
            particles.push_back(part);
        }

        if constexpr (Interpolator::interp_order == 2)
        {
            part.iCell[0] = 19; // AMR index
            part.delta[0] = 0.0;
            part.weight   = 1.0;
            part.charge   = 2.0;
            part.v[0]     = +2.;
            part.v[1]     = -1.;
            part.v[2]     = +1.;
            particles.push_back(part);

            part.iCell[0] = 20; // AMR index
            part.delta[0] = 0.0;
            part.weight   = 0.2;
            part.charge   = 3.2;
            part.v[0]     = +2.;
            part.v[1]     = -1.;
            part.v[2]     = +1.;
            particles.push_back(part);

            part.iCell[0] = 20; // AMR index
            part.delta[0] = 0.0;
            part.weight   = 0.8;
            part.charge   = 1.7;
            part.v[0]     = +2.;
            part.v[1]     = -1.;
            part.v[2]     = +1.;
            particles.push_back(part);

            part.iCell[0] = 21; // AMR index
            part.delta[0] = 0.0;
            part.weight   = 1.0;
            part.charge   = 2.0;
            part.v[0]     = +2.;
            part.v[1]     = -1.;
            part.v[2]     = +1.;
            particles.push_back(part);
        }

        if constexpr (Interpolator::interp_order == 3)
        {
            part.iCell[0] = 18; // AMR index
            part.delta[0] = 0.5;
            part.weight   = 1.0;
            part.charge   = 2.0;
            part.v[0]     = +2.;
            part.v[1]     = -1.;
            part.v[2]     = +1.;
            particles.push_back(part);

            part.iCell[0] = 19; // AMR index
            part.delta[0] = 0.5;
            part.weight   = 1.0;
            part.charge   = 2.0;
            part.v[0]     = +2.;
            part.v[1]     = -1.;
            part.v[2]     = +1.;
            particles.push_back(part);

            part.iCell[0] = 20; // AMR index
            part.delta[0] = 0.5;
            part.weight   = 1.0;
            part.charge   = 2.0;
            part.v[0]     = +2.;
            part.v[1]     = -1.;
            part.v[2]     = +1.;
            particles.push_back(part);

            part.iCell[0] = 21; // AMR index
            part.delta[0] = 0.5;
            part.weight   = 0.1;
            part.charge   = 3.35;
            part.v[0]     = +2.;
            part.v[1]     = -1.;
            part.v[2]     = +1.;
            particles.push_back(part);

            part.iCell[0] = 21; // AMR index
            part.delta[0] = 0.5;
            part.weight   = 0.9;
            part.charge   = 1.85;
            part.v[0]     = +2.;
            part.v[1]     = -1.;
            part.v[2]     = +1.;
            particles.push_back(part);
        }
        interpolator(makeIndexRange(particles), rho, rho_c, v, layout);
    }



protected:
    Interpolator interpolator;
};
TYPED_TEST_SUITE_P(ACollectionOfParticles_1d);

TYPED_TEST_P(ACollectionOfParticles_1d, DepositCorrectlyTheirWeight_1d)
{
    constexpr auto interp = TypeParam::interp_order;

    auto idx = 20 + this->layout.nbrGhosts(QtyCentering::dual);

    auto const& [vx, vy, vz] = this->v();
    EXPECT_DOUBLE_EQ(this->rho(idx), 1.0);
    EXPECT_DOUBLE_EQ(this->rho_c(idx), 2.0);
    EXPECT_DOUBLE_EQ(vx(idx), 2.0);
    EXPECT_DOUBLE_EQ(vy(idx), -1.0);
    EXPECT_DOUBLE_EQ(vz(idx), 1.0);
}
REGISTER_TYPED_TEST_SUITE_P(ACollectionOfParticles_1d, DepositCorrectlyTheirWeight_1d);

using MyTypes = ::testing::Types<Interpolator<1, 1>, Interpolator<1, 2>, Interpolator<1, 3>>;
INSTANTIATE_TYPED_TEST_SUITE_P(testInterpolator, ACollectionOfParticles_1d, MyTypes);


template<typename Interpolator>
struct ACollectionOfParticles_2d : public ::testing::Test
{
    static constexpr auto interp_order = Interpolator::interp_order;
    static constexpr std::size_t dim   = 2;
    static constexpr std::uint32_t nx = 15, ny = 15;
    static constexpr int start = 0, end = 5;
    static constexpr auto safeLayer = static_cast<int>(1 + ghostWidthForParticles<interp_order>());

    using PHARE_TYPES      = PHARE::core::PHARE_Types<dim, interp_order>;
    using ParticleArray_t  = typename PHARE_TYPES::ParticleArray_t;
    using GridLayout_t     = typename PHARE_TYPES::GridLayout_t;
    using Grid_t           = typename PHARE_TYPES::Grid_t;
    using UsableVecFieldND = UsableVecField<dim>;

    GridLayout_t layout{ConstArray<double, dim>(.1), {nx, ny}, ConstArray<double, dim>(0)};
    ParticleArray_t particles;
    Grid_t rho;
    Grid_t rho_c;
    UsableVecFieldND v;
    Interpolator interpolator;

    ACollectionOfParticles_2d()
        : particles{grow(layout.AMRBox(), safeLayer)}
        , rho{"field", HybridQuantity::Scalar::rho, nx, ny}
        , rho_c{"field", HybridQuantity::Scalar::rho, nx, ny}
        , v{"v", layout, HybridQuantity::Vector::V}
    {
        for (int i = start; i < end; i++)
            for (int j = start; j < end; j++)
            {
                auto& part  = particles.emplace_back();
                part.iCell  = {i, j};
                part.delta  = ConstArray<double, dim>(.5);
                part.weight = 1.;
                part.v[0]   = +2.;
                part.v[1]   = -1.;
                part.v[2]   = +1.;
            }
        interpolator(makeIndexRange(particles), rho, rho_c, v, layout);
    }
};
TYPED_TEST_SUITE_P(ACollectionOfParticles_2d);


TYPED_TEST_P(ACollectionOfParticles_2d, DepositCorrectlyTheirWeight_2d)
{
    constexpr auto interp = TypeParam::interp_order;

    auto idx                 = 2 + this->layout.nbrGhosts(QtyCentering::dual);
    auto const& [vx, vy, vz] = this->v();
    EXPECT_DOUBLE_EQ(this->rho(idx, idx), 1.0);
    EXPECT_DOUBLE_EQ(vx(idx, idx), 2.0);
    EXPECT_DOUBLE_EQ(vy(idx, idx), -1.0);
    EXPECT_DOUBLE_EQ(vz(idx, idx), 1.0);
}
REGISTER_TYPED_TEST_SUITE_P(ACollectionOfParticles_2d, DepositCorrectlyTheirWeight_2d);


using My2dTypes = ::testing::Types<Interpolator<2, 1>, Interpolator<2, 2>, Interpolator<2, 3>>;
INSTANTIATE_TYPED_TEST_SUITE_P(testInterpolator, ACollectionOfParticles_2d, My2dTypes);


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
