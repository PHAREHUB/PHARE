#include "phare_core.hpp"
#include "core/data/particles/particle.hpp"
#include "core/data/particles/particle_array.hpp"
#include "core/data/particles/particle_array_sorter.hpp"
#include "core/data/particles/particle_array_converter.hpp"

#include "tests/core/data/gridlayout/test_gridlayout.hpp"

#include "gtest/gtest.h"

using namespace PHARE;
using namespace PHARE::core;


template<std::size_t _dim, auto _layout_mode, auto _alloc_mode = AllocatorMode::CPU>
struct TestParam
{
    static_assert(all_are<LayoutMode>(_layout_mode));
    static_assert(all_are<AllocatorMode>(_alloc_mode));
    auto constexpr static dim         = _dim;
    auto constexpr static layout_mode = _layout_mode;
    auto constexpr static alloc_mode  = _alloc_mode;

    using ParticleArray_t
        = ParticleArray<ParticleArrayOptions{dim, layout_mode, StorageMode::VECTOR, alloc_mode}>;
};



template<typename Param>
class ParticleArraySortingTest : public ::testing::Test
{
    static constexpr std::size_t dim = 3;
    auto constexpr static interp     = 1;
    auto constexpr static sim_opts
        = SimOpts{.dimension = dim, .interp_order = interp, .layout_mode = Param::layout_mode,
                 .alloc_mode = Param::alloc_mode};
    using Particle_t       = Particle<dim>;
    using box_t            = Box<int, dim>;
    using GridLayout_t     = PHARE_Types<sim_opts>::GridLayout_t;
    using TestGridLayout_t = TestGridLayout<GridLayout_t>;

public:
    using ParticleArray_t    = Param::ParticleArray_t;
    using RefParticleArray_t = AoSParticleArray<dim>;

    void add(std::array<int, dim> const& iCell)
    {
        particles.push_back(Particle_t{0.01, 1., iCell, {0.002f, 0.2f, 0.8f}, {1.8, 1.83, 2.28}});
    }



    void print()
    {
        for (auto const& particle : particles)
            std::cout << __LINE__ << " " << Point{particle.iCell()} << std::endl;
        std::cout << std::endl;
    }

    box_t domain_box{{0, 0, 0}, {3, 3, 3}};
    TestGridLayout_t layout{domain_box};
    ParticleArray_t particles{make_particles<ParticleArray_t>(layout)};
};


// clang-format off
using Permutations_t = testing::Types< // ! notice commas !

    TestParam<3, LayoutMode::SoA>
   ,TestParam<3, LayoutMode::AoS>
   ,TestParam<3, LayoutMode::AoSTS>

PHARE_WITH_GPU(

   ,TestParam<3, LayoutMode::SoA, AllocatorMode::GPU_UNIFIED>
   ,TestParam<3, LayoutMode::AoS, AllocatorMode::GPU_UNIFIED>
   ,TestParam<3, LayoutMode::AoSTS, AllocatorMode::GPU_UNIFIED>

)

>;
// clang-format on


TYPED_TEST_SUITE(ParticleArraySortingTest, Permutations_t, );

TYPED_TEST(ParticleArraySortingTest, _3d_sorting_test)
{
    using ParticleArray_t    = TestFixture::ParticleArray_t;
    using RefParticleArray_t = TestFixture::RefParticleArray_t;

    PHARE_LOG_LINE_SS(ParticleArray_t::type_id);

    auto& particles = this->particles;
    std::vector<std::array<int, 3>> iCells{{2, 2, 2}, {2, 1, 0}, {1, 1, 1}, //
                                           {1, 0, 1}, {0, 2, 0}, {0, 0, 0},
                                           {1, 1, 2}, {0, 1, 0}, {2, 0, 2}};
    for (auto const& iCell : iCells)
        this->add(iCell);

    std::vector<std::array<int, 3>> expected{{0, 0, 0}, {0, 1, 0}, {0, 2, 0}, //
                                             {1, 0, 1}, {1, 1, 1}, {1, 1, 2},
                                             {2, 0, 2}, {2, 1, 0}, {2, 2, 2}};
    ASSERT_EQ(particles.size(), expected.size());

    auto const& ps = convert_particles_and_sort<RefParticleArray_t>(particles, this->layout);

    for (std::size_t i = 0; i < particles.size(); ++i)
        ASSERT_EQ(ps.iCell(i), expected[i]);
}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
