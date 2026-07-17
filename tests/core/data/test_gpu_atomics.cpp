// test_updater_pp_main.cpp

#include "core/numerics/ion_updater/ion_updater.hpp"
#include "core/numerics/ion_updater/ion_updater_pc.hpp"
#include "core/numerics/ion_updater/ion_updater_per_particle.hpp"

#include "tests/core/data/electromag/test_electromag_fixtures.hpp"
#include "tests/core/data/ion_population/test_ion_population_fixtures.hpp"

#include "gtest/gtest.h"

#include <cooperative_groups.h> // nvidia/cuda
#include <unordered_set>

namespace PHARE::core
{

template<auto opts, typename internals>
void PrintTo(ParticleArray<opts, internals> const& arr, std::ostream* os)
{
    assert(arr.size());
    *os << arr;
}

auto static const ppc     = get_env_as("PHARE_PPC", std::size_t{64});
auto static const seed    = get_env_as("PHARE_SEED", std::size_t{8});
auto static const cells   = get_env_as("PHARE_CELLS", std::uint32_t{5});
auto static const dt      = get_env_as("PHARE_TIMESTEP", double{.001});
auto static const shufle  = get_env_as("PHARE_UNSORTED", std::size_t{0});
bool static const premain = []() {
    PHARE_WITH_PHLOP(          //
        using namespace PHARE; //
        using namespace std::literals;
        if (auto e = core::get_env("PHARE_SCOPE_TIMING", "false"); e == "1" || e == "true")
            phlop::ScopeTimerMan::INSTANCE()
                .file_name(".phare_times.0.csv")
                .force_strings()
                .headers("fn"s, "dim"s, "layout"s, "alloc"s, "storage"s, "time"s)
                .init(); //
    )
    return true;
}();




template<typename Particles_t, typename GridLayout_t>
auto make_ions(GridLayout_t const& layout)
{
    auto constexpr static alloc_mode = Particles_t::alloc_mode;
    auto constexpr static dim        = GridLayout_t::dimension;
    auto constexpr static interp     = GridLayout_t::interp_order;

    auto ions_p = std::make_shared<UsableIons_t<Particles_t, interp>>(layout, "protons");
    auto& ions  = *ions_p;

    EXPECT_EQ(ions.populations[0].particles.domain_particles.size(), 0);

    auto disperse = [](auto& particles) {
        // delta_disperse(particles.domain_particles, seed);
        // delta_disperse(particles.patch_ghost_particles, seed);
        if (shufle > 0)
        {
            // shuffle(particles.domain_particles, seed);
            // shuffle(particles.patch_ghost_particles, seed);
        }
    };

    auto add_particles = [&](auto& particles) {
        add_particles_in(particles.domain_particles, layout.AMRBox(), ppc);
        add_ghost_particles(particles.patch_ghost_particles, layout.AMRBox(), ppc,
                            GridLayout_t::nbrParticleGhosts());
    };

    add_particles(ions.populations[0].particles);
    disperse(ions.populations[0].particles);

    EXPECT_EQ(ions.populations[0].particles.domain_particles.size(), layout.AMRBox().size() * ppc);

    return ions_p;
}



template<std::size_t _dim, auto _layout_mode, auto _alloc_mode = PHARE::AllocatorMode::CPU>
struct TestParam
{
    static_assert(std::is_same_v<decltype(_layout_mode), LayoutMode>);
    static_assert(std::is_same_v<decltype(_alloc_mode), PHARE::AllocatorMode>);
    auto constexpr static dim         = _dim;
    auto constexpr static layout_mode = _layout_mode;
    auto constexpr static alloc_mode  = _alloc_mode;
};



template<typename Param>
struct AtomicsTest : public ::testing::Test
{
    auto constexpr static dim         = Param::dim;
    auto constexpr static interp      = 1;
    auto constexpr static layout_mode = Param::layout_mode;
    auto constexpr static alloc_mode  = Param::alloc_mode;

    using GridLayout_t = TestGridLayout<typename PHARE_Types<SimOpts{dim, interp}>::GridLayout_t>;
    using ParticleArray_t
        = ParticleArray<ParticleArrayOptions{dim, layout_mode, StorageMode::VECTOR, alloc_mode}>;

    using Ions_t = UsableIons_t<ParticleArray_t, interp>;

    GridLayout_t const layout{cells};
    std::shared_ptr<Ions_t> ions = make_ions<ParticleArray_t>(*layout);
};

// clang-format off
using Permutations_t = testing::Types< // ! notice commas !
     // TestParam<1, LayoutMode::AoS>                               // 0
    // ,TestParam<1, LayoutMode::AoSPC>                             // 1
PHARE_WITH_MKN_GPU(
    // ,TestParam<1, LayoutMode::AoS, AllocatorMode::GPU_UNIFIED>   // 3
    /*,*/TestParam<1, LayoutMode::AoSPC, AllocatorMode::GPU_UNIFIED> // 4
)
    // ,TestParam<2, LayoutMode::AoS>
    // ,TestParam<2, LayoutMode::AoSPC>
// PHARE_WITH_MKN_GPU(
//     ,TestParam<2, LayoutMode::SoA>
//     ,TestParam<2, LayoutMode::AoS, AllocatorMode::GPU_UNIFIED>
//     ,TestParam<2, LayoutMode::AoSPC, AllocatorMode::GPU_UNIFIED> // 10
//     ,TestParam<2, LayoutMode::SoA, AllocatorMode::GPU_UNIFIED>
// )
//     ,TestParam<3, LayoutMode::AoS>
//     ,TestParam<3, LayoutMode::AoSPC>
// PHARE_WITH_MKN_GPU(
//     ,TestParam<3, LayoutMode::SoA>
//     ,TestParam<3, LayoutMode::AoS, AllocatorMode::GPU_UNIFIED>
//     ,TestParam<3, LayoutMode::AoSPC, AllocatorMode::GPU_UNIFIED> // 16
//     ,TestParam<3, LayoutMode::SoA, AllocatorMode::GPU_UNIFIED>
// )
>;
// clang-format on

TYPED_TEST_SUITE(AtomicsTest, Permutations_t);

__device__ int atomicAggInc2(int* ctr)
{
    namespace cg = cooperative_groups;
    auto g       = cg::coalesced_threads();
    int warp_res;
    if (g.thread_rank() == 0)
        warp_res = atomicAdd(ctr, g.size());
    return g.shfl(warp_res, 0) + g.thread_rank();
}
__device__ int atomicAggInc1(int* ptr)
{
    namespace cg = cooperative_groups;
    auto g       = cg::coalesced_threads();
    int prev;
    // elect the first active thread to perform atomic add
    if (g.thread_rank() == 0)
    {
        prev = atomicAdd(ptr, g.size());
    }
    // broadcast previous value within the warp
    // and add each active thread’s rank to it
    prev = g.thread_rank() + g.shfl(prev, 0);
    return prev;
}
/*
TYPED_TEST(AtomicsTest, incrementWorksGlobal)
{
  using ManagedVector = std::vector<std::size_t, mkn::gpu::ManagedAllocator<std::size_t>>;

  PHARE_WITH_MKN_GPU(         //
      ManagedVector mv(ppc);
      auto& particles = this->ions->populations[0].particles.domain_particles; //

      auto view = *particles; //
      auto v = mv.data();

      auto needle = ConstArray<int, TypeParam::dim>(3);

      std::size_t count = 0;
      for(auto const& p : particles)
        if(array_equals(p.iCell(), needle)) ++count;
      EXPECT_EQ(count, ppc);

      mkn::gpu::GDLauncher{particles.size()}([=] _PHARE_ALL_FN_() mutable {
          __syncthreads();
          auto particle = view[mkn::gpu::idx()];
          auto b = array_equals(particle.iCell(), needle);
          printf("addr %d %d %lu %lu \n", b, particle.iCell()[0], &view.idx, mkn::gpu::idx());
          // __syncthreads();
          if(b) {
              auto idx = atomicAggInc1(&view.idx);
              v[idx] = mkn::gpu::idx();
          } //
          __syncthreads();
      }); //

      for(std::size_t i = 0; i < mv.size(); ++i){
        PHARE_LOG_LINE_STR(mv[i]);
      }
      std::size_t idx = mv[0];
      for(std::size_t i = 1; i < mv.size(); ++i){
          EXPECT_EQ(mv[i], ++idx);
      }
  )
}
*/

template<typename F, typename... Args>
__global__ static void my_kernel(F f)
{
    /*if (auto i = mkn::gpu::cuda::idx(); i < s)*/
    f();
}

// TYPED_TEST(AtomicsTest, incrementWorks)
// {
//   using ManagedVector = std::vector<std::size_t, mkn::gpu::ManagedAllocator<std::size_t>>;

//   PHARE_WITH_MKN_GPU(         //
//       ManagedVector mv(ppc);
//       std::vector<int, mkn::gpu::ManagedAllocator<int>> minfo(5);
//       auto& particles = this->ions->populations[0].particles.domain_particles; //

//       auto view = *particles; //
//       auto v = mv.data();
//       auto info = minfo.data();

//       auto needle = ConstArray<int, TypeParam::dim>(3);

//       // auto lambda = [=] _PHARE_ALL_FN_() mutable {
//       //     if(auto tid = mkn::gpu::idx(); view[mkn::gpu::idx()].iCell()[0] == needle[0]) {
//       //         auto idx = atomicAdd(&info[0], 1);
//       //         v[idx] = tid;
//       //     }
//       // };
//       // mkn::gpu::Launcher{512,1,512,1}(my_kernel<decltype(lambda)>, lambda); //

//           // __syncthreads();
//       mkn::gpu::GDLauncher{particles.size()}([=] _PHARE_ALL_FN_() mutable {

//           // if(auto tid = mkn::gpu::idx(); view[tid].iCell()[0] == needle[0]) {
//           //     auto idx = atomicAdd(&info[0], 1);
//           //     v[idx] = tid;
//           // }

//           auto tid = mkn::gpu::idx();
//           auto particle = view[tid];
//           if(/*view[tid].iCell()[0] == needle[0]*/
//           array_equals(particle.iCell(), needle)) {
//               auto idx = atomicAdd(&view.idx, 1);
//               v[idx] = tid;
//           }
//       }); //

//       std::size_t idx = mv[0];
//       for(std::size_t i = 1; i < mv.size(); ++i){
//           EXPECT_EQ(mv[i], ++idx);
//       }
//   )
// }



TYPED_TEST(AtomicsTest, icellChangeWorks)
{
    using ManagedVector = std::vector<std::size_t, mkn::gpu::ManagedAllocator<std::size_t>>;

    PHARE_WITH_MKN_GPU( //

        auto& particles = this->ions->populations[0].particles.domain_particles;
        auto view       = *particles;

        mkn::gpu::GDLauncher{particles.size()}([=] _PHARE_ALL_FN_() mutable {
            auto tid      = mkn::gpu::idx();
            auto particle = view[tid];

            if (particle.iCell()[0] == 3 and tid % 2 == 0)
            {
                auto newCell = particle.iCell();
                newCell[0] += 1;
                particle.icell_changer(newCell);
                particle.iCell() = newCell;
            }
        });

        ParticleArrayService::sync<0, ParticleType::Domain>(particles);

        std::unordered_set<std::size_t> ids;
        for (auto const& bix : particles.local_box()) for (auto const& p : particles(bix)) {
            if (ids.count(p.id) != 0)
            {
                PHARE_LOG_LINE_STR("DUPLICATE ID! " << p.id);
            }
            EXPECT_EQ(ids.count(p.id), 0);
            ids.emplace(p.id);
            EXPECT_EQ(ids.count(p.id), 1);
        }
        //
    )
}


} // namespace PHARE::core


int main(int argc, char** argv)
{
    PHARE_LOG_LINE_STR(PHARE::core::seed);
    assert(phlop::ScopeTimerMan::INSTANCE().active);
    ::testing::InitGoogleTest(&argc, argv);
    auto r = RUN_ALL_TESTS();
    PHARE_WITH_PHLOP(                  //
        phlop::ScopeTimerMan::reset(); //
    )
    return r;
}
