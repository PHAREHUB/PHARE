
#include "phare_core.hpp"

#include "core/utilities/geom.hpp"
#include "tests/core/utilities/box/test_box_fixtures.hpp"
#include "tests/core/data/particles/test_particles_fixtures.hpp"

#include "gtest/gtest.h"

using namespace PHARE::core;

std::size_t constexpr static dim = 3; // only 3d
using Box_t                      = PHARE::core::Box<int, dim>;
using FileParticles              = ParticleArray<ParticleArrayOptions{dim, LayoutMode::AoS}>;

struct StaticParticlesDataFixture
{
    std::size_t constexpr static ghost_cells = 2;
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
        disperse(particles, 12345);
        std::cout << "Allocated RAM for particles: " << ram_in_gbs(particles) << std::endl;
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
    Box<int, dim> ghostBox = grow(box, ghost_cells);
    Particles domainParticles;
};

// PP = per particle
template<typename Particles_>
struct ParticlesDataFixture
{
    std::size_t constexpr static ghost_cells = StaticParticlesDataFixture::ghost_cells;
    std::size_t constexpr static ppc         = StaticParticlesDataFixture::ppc;
    std::size_t constexpr static box_size    = StaticParticlesDataFixture::box_size;
    std::size_t const in_box                 = std::pow(box_size, Particles_::dimension) * ppc;

    using Particles_t     = Particles_;
    using ParticlesData_t = ParticlesData<Particles_t, ghost_cells>;

    ParticlesDataFixture()
        : middle{StaticParticlesDataFixture::instance().middle()}
        , neighbours{
              generate_from([&](auto& box) { return std::make_unique<ParticlesData_t>(box); },
                            StaticParticlesDataFixture::instance().neighbours())}
    {
        auto add_particles_in_ = [](auto& pdata) {
            // for (auto const& particle :
            //      read_raw_from_file<FileParticles>(file_name_for_box(pdata.box)))
            for (auto const& particle : StaticParticlesDataFixture::instance().particles)
                pdata.domainParticles.push_back(particle);
        };
        add_particles_in_(middle);
    }


    ParticlesData_t middle;
    std::vector<std::unique_ptr<ParticlesData_t>> neighbours;
};


template<typename ParticlesData>
struct GeomTest : public ::testing::Test, public ParticlesData
{
};

using ParticlesDatas = testing::Types< //
    ParticlesDataFixture<ParticleArray<ParticleArrayOptions{dim, LayoutMode::AoS}>>,
    ParticlesDataFixture<ParticleArray<ParticleArrayOptions{
        dim, LayoutMode::SoA, StorageMode::VECTOR, PHARE::AllocatorMode::GPU_UNIFIED}>>,
    ParticlesDataFixture<ParticleArray<ParticleArrayOptions{
        dim, LayoutMode::AoS, StorageMode::VECTOR, PHARE::AllocatorMode::GPU_UNIFIED}>>
    /*,
ParticlesDataFixture<ParticleArray< ParticleArrayOptions{dim, LayoutMode::AoSMapped}>>,
ParticlesDataFixture<ParticleArray< ParticleArrayOptions{dim, LayoutMode::SoA}>>,

*/

    >;

TYPED_TEST_SUITE(GeomTest, ParticlesDatas, );

template<std::size_t impl, typename Fixture>
auto count(Fixture const& fixture)
{
    auto& middle       = fixture.middle;
    auto& particles    = middle.domainParticles;
    std::size_t in_box = 0;
    auto start         = now_in_microseconds();
    PointInBoxChecker<Box_t> checker{middle.box};
    for (auto const& p : particles)
        if (checker(p.iCell()))
            ++in_box;
    auto end = now_in_microseconds();
    PHARE_LOG_LINE_STR("time " << end - start);
    return in_box;
}

TYPED_TEST(GeomTest, works0)
{
    EXPECT_EQ(count<0>(*this), this->in_box);
}
TYPED_TEST(GeomTest, works1)
{
    EXPECT_EQ(count<1>(*this), this->in_box);
}
TYPED_TEST(GeomTest, works2)
{
    EXPECT_EQ(count<2>(*this), this->in_box);
}
TYPED_TEST(GeomTest, works3)
{
    EXPECT_EQ(count<3>(*this), this->in_box);
}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
