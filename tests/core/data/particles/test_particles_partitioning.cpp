
#include "phare_core.hpp"

#include "core/utilities/geom.hpp"
// #include "amr/utilities/box/amr_box.hpp"
//  #include "amr/data/particles/particles_data.hpp"

#include "tests/core/data/vecfield/test_vecfield_fixtures.hpp"
#include "tests/core/data/electromag/test_electromag_fixtures.hpp"
#include "tests/core/data/tensorfield/test_tensorfield_fixtures.hpp"
#include "tests/core/data/ion_population/test_ion_population_fixtures.hpp"

#include "tests/core/utilities/box/test_box_fixtures.hpp"

#include "gtest/gtest.h"

using namespace PHARE::core;



std::size_t constexpr static dim = 3; // only 3d
using Box_t                      = PHARE::core::Box<int, dim>;
using FileParticles              = ParticleArray<ParticleArrayOptions{dim, LayoutMode::AoS}>;

auto file_name_for_box(Box_t const& box)
{
    std::stringstream ss;
    ss << "particles." << std::hex << std::hash<std::string>()(to_string(box)) << ".bin";
    return ss.str();
}

struct StaticParticlesDataFixture
{
    std::size_t constexpr static ghost_cells = 1;
    std::size_t constexpr static ppc         = 125;
    std::size_t constexpr static box_size    = 50;

    static StaticParticlesDataFixture instance()
    {
        static StaticParticlesDataFixture i;
        return i;
    }

    auto b() const { return true; }
    auto& middle() const { return std::get<0>(box_tuple); }
    auto& neighbours() const { return std::get<1>(box_tuple); }

    StaticParticlesDataFixture()
        : box_tuple{get_some_boxes(box_size)}
    {
        add_particles_in(particles, grow(middle(), ghost_cells), ppc);
        // disperse(particles, 13337);
    }

    std::tuple<Box_t, std::vector<Box_t>> const box_tuple;
    FileParticles particles;
};
static auto is_initialized = StaticParticlesDataFixture::instance().b();

template<typename Particles, std::size_t ghost_cells = 2>
struct ParticlesData
{
    auto static constexpr dim = Particles::dimension;
    ParticlesData(Box<int, dim> const& box_)
        : box{box_}
    {
    }
    Box<int, dim> box;
    Box<int, dim> ghost_box = grow(box, ghost_cells);
    Particles domainParticles;
};


template<typename Particles_, std::size_t impl_ = 0>
struct ParticlesDataFixture
{
    std::size_t constexpr static ghost_cells = StaticParticlesDataFixture::ghost_cells;
    std::size_t constexpr static ppc         = StaticParticlesDataFixture::ppc;
    std::size_t constexpr static box_size    = StaticParticlesDataFixture::box_size;
    std::size_t constexpr static impl        = impl_;

    using Particles_t     = Particles_;
    using ParticlesData_t = ParticlesData<Particles_t, ghost_cells>;

    ParticlesDataFixture()
        : middle{StaticParticlesDataFixture::instance().middle()}
        , neighbours{
              generate_from([&](auto& box) { return std::make_unique<ParticlesData_t>(box); },
                            StaticParticlesDataFixture::instance().neighbours())}
    {
        auto& particles = middle.domainParticles;
        for (auto const& particle : StaticParticlesDataFixture::instance().particles)
            particles.push_back(particle);

        // particles.norcell.resize(particles.size());
        // for (std::size_t i = 0; i < particles.size(); ++i)
        //     particles[i].norcell
        //         = normalized_icell_relative_to_closest_edge(middle.box, particles.iCell(i));
    }

    ParticlesData_t middle;
    std::vector<std::unique_ptr<ParticlesData_t>> neighbours;
};


template<typename ParticlesData>
struct ParticlesDataTest : public ::testing::Test, public ParticlesData
{
};


using ParticlesDatas = testing::Types< //
    ParticlesDataFixture<ParticleArray<ParticleArrayOptions{dim, LayoutMode::AoS}>>
    //, ParticlesDataFixture<ParticleArray< ParticleArrayOptions{dim, LayoutMode::AoS}>, 1>
    /*,
     ParticlesDataFixture<
         ParticleArray< ParticleArrayOptions{dim, LayoutMode::SoA, StorageMode::VECTOR,
                                                   PHARE::AllocatorMode::GPU_UNIFIED}>>,
     ParticlesDataFixture<
         ParticleArray< ParticleArrayOptions{dim, LayoutMode::AoS, StorageMode::VECTOR,
                                                   PHARE::AllocatorMode::GPU_UNIFIED}>>*/
    /*,
ParticlesDataFixture<ParticleArray< ParticleArrayOptions{dim, LayoutMode::AoSMapped}>>,
ParticlesDataFixture<ParticleArray< ParticleArrayOptions{dim, LayoutMode::SoA}>>,

*/

    >;

TYPED_TEST_SUITE(ParticlesDataTest, ParticlesDatas, );

TYPED_TEST(ParticlesDataTest, works)
{
    using Test            = TypeParam;
    using ParticleArray_t = typename TestFixture::Particles_t;
    using Partitioner     = ParticleArrayPartitioner<ParticleArray_t::alloc_mode, ParticleArray_t>;
    auto& middle          = this->middle;
    auto& particles       = middle.domainParticles;

    for (std::size_t i = 0; i < particles.size(); ++i)
    {
        EXPECT_TRUE(isIn(particles.iCell(i), middle.ghost_box));
        // if (isIn(particles.iCell(i), middle.box))
        // {
        //     EXPECT_TRUE(particles[i].norcell >= 0);
        // }
    }

    auto before_size = particles.size();
    auto start       = now_in_microseconds();
    auto iterator    = Partitioner{particles}(middle.box);
    auto end         = now_in_microseconds();

    PHARE_LOG_LINE_STR("time " << end - start);
    EXPECT_EQ(before_size, particles.size());
    EXPECT_EQ(iterator.size(), std::pow(Test::box_size, dim) * Test::ppc);

    for (std::size_t i = 0; i < iterator.size(); ++i)
        EXPECT_TRUE(isIn(particles.iCell(i), middle.box));

    for (std::size_t i = iterator.size(); i < particles.size(); ++i)
    {
        auto const& iCell = particles.iCell(i);
        EXPECT_TRUE(isIn(iCell, middle.ghost_box));
        EXPECT_TRUE(not isIn(iCell, middle.box));
    }
}


TYPED_TEST(ParticlesDataTest, worksVector)
{
    using Test            = TypeParam;
    using ParticleArray_t = typename TestFixture::Particles_t;
    using Partitioner     = ParticleArrayPartitioner<ParticleArray_t::alloc_mode, ParticleArray_t>;
    auto& middle          = this->middle;
    auto& particles       = middle.domainParticles;

    for (std::size_t i = 0; i < particles.size(); ++i)
    {
        EXPECT_TRUE(isIn(particles.iCell(i), middle.ghost_box));
    }

    auto before_size = particles.size();
    auto start       = now_in_microseconds();
    auto iterators   = Partitioner{particles}(std::vector{middle.box});
    auto end         = now_in_microseconds();

    auto iterator = iterators.back();

    PHARE_LOG_LINE_STR("time " << end - start);
    EXPECT_EQ(before_size, particles.size());
    EXPECT_EQ(iterator.size(), std::pow(Test::box_size, dim) * Test::ppc);

    for (std::size_t i = 0; i < iterator.size(); ++i)
        EXPECT_TRUE(isIn(particles.iCell(i), middle.box));

    for (std::size_t i = iterator.size(); i < particles.size(); ++i)
    {
        auto const& iCell = particles.iCell(i);
        EXPECT_TRUE(isIn(iCell, middle.ghost_box));
        EXPECT_TRUE(not isIn(iCell, middle.box));
    }
}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
