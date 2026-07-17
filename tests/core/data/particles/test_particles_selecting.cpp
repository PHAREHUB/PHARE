

#include "phare_core.hpp"
#include "core/utilities/types.hpp"
#include "core/data/particles/particle_array.hpp"
#include "core/data/particles/particle_array_service.hpp"
#include "core/data/particles/particle_array_selector.hpp"

#include "tests/core/data/particles/test_particles.hpp"
#include "tests/core/data/gridlayout/test_gridlayout.hpp"

#include "gtest/gtest.h"


namespace PHARE::core
{

auto static const bytes   = get_env_as("PHARE_GPU_BYTES", std::uint64_t{500000000});
auto static const cells   = get_env_as("PHARE_CELLS", std::uint32_t{4});
auto static const ppc     = get_env_as("PHARE_PPC", std::size_t{1});
bool static const premain = []() {
    //     PHARE_WITH_MKN_GPU({
    //         PHARE_LOG_LINE_STR("bytes: " << bytes); //
    //         mkn::gpu::setDevice(0);
    //         mkn::gpu::prinfo();
    //         mkn::gpu::setLimitMallocHeapSize(bytes);
    //         mkn::gpu::print_gpu_mem_used();
    //     })
    //     PHARE_WITH_PHLOP(                           //
    //         PHARE_LOG_LINE_STR("cells: " << cells); //
    //         PHARE_LOG_LINE_STR("ppc  : " << ppc);   //

    //         using namespace PHARE; //
    //         using namespace std::literals;
    //         if (auto e = core::get_env("PHARE_SCOPE_TIMING", "false"); e == "1" || e == "true")
    //             phlop::threaded::ScopeTimerMan::INSTANCE()
    //                 .file_name(".phare_times.0.txt")
    //                 // .force_strings()
    //                 // .headers("fn"s, "dim"s, "layout"s, "alloc"s, "storage"s, "time"s)
    //                 .init(); //
    //     )
    return true;
}();

template<std::size_t _dim, auto lm, auto am>
struct TestParam
{
    static_assert(all_are<LayoutMode>(lm));
    static_assert(all_are<AllocatorMode>(am));

    auto constexpr static dim         = _dim;
    auto constexpr static layout_mode = lm;
    auto constexpr static alloc_mode  = am;
};


template<typename Param>
struct AParticleArraySelectingTest : public ::testing::Test
{
    auto constexpr static dim         = Param::dim;
    auto constexpr static interp      = 1;
    auto constexpr static layout_mode = Param::layout_mode;
    auto constexpr static alloc_mode  = Param::alloc_mode;
    auto constexpr static sim_opts
        = SimOpts{.dimension = dim, .interp_order = interp, .layout_mode = Param::layout_mode,
                 .alloc_mode = Param::alloc_mode};

    using GridLayout_t     = PHARE_Types<sim_opts>::GridLayout_t;
    using TestGridLayout_t = TestGridLayout<GridLayout_t>;
    using Box_t            = Box<int, dim>;
    using ParticleArray_t
        = ParticleArray<ParticleArrayOptions{dim, layout_mode, StorageMode::VECTOR, alloc_mode}>;

    AParticleArraySelectingTest() = default;


    struct Patch
    {
        Patch(TestGridLayout_t const& _layout)
            : layout{_layout}
        {
            add_particles_in(domain, layout.AMRBox(), ppc);
            abort_if_not(domain.size() == ppc * layout.AMRBox().size());
        }


        GridLayout_t layout;
        ParticleArray_t domain = make_particles<ParticleArray_t>(layout);
        ParticleArray_t ghost  = make_particles<ParticleArray_t>(layout);
    };

    bool at_periodic_border(Box<int, dim> const& box) const
    {
        auto const ghostbox = grow(box, 1);
        return domainBox * ghostbox != ghostbox;
    }

    auto& neighbours_for(std::size_t const pid)
    {
        neighbours.clear();

        auto const ghostbox = grow(patches[pid].layout.AMRBox(), 1);
        auto const emplace  = [&](auto from, auto to) {
            for (std::size_t i = from; i < to; ++i)
                if (auto const overlap = ghostbox * patches[i].layout.AMRBox())
                    neighbours.emplace_back(&patches[i], *overlap);
        };

        emplace(0, pid);
        emplace(pid + 1, patches.size());

        return neighbours;
    }


    auto& periodic_neighbours_for(std::size_t const pid)
    {
        periodic_neighbours.clear();
        if (!at_periodic_border(patches[pid].layout.AMRBox()))
            return periodic_neighbours;

        auto const& patch = patches[pid];
        int const icells  = cells;
        int const mid     = icells * 3 / 2;
        auto shifts
            = for_N<7, for_N_R_mode::make_array>([&](auto i) { return Point<int, 3>{0, 0, 0}; });
        for_N<3>([&](auto i) {
            int const span  = icells * 3;
            int const shift = patch.layout.AMRBox().upper[i] < mid ? 1 : -1;
            shifts[i][i]    = span * shift;
        });

        shifts[3] = {shifts[0][0], shifts[1][1], 0};
        shifts[4] = {0, shifts[1][1], shifts[2][2]};
        shifts[5] = {shifts[0][0], 0, shifts[2][2]};
        shifts[6] = {shifts[0][0], shifts[1][1], shifts[2][2]};

        for (auto const& shifter : shifts)
        {
            auto const shift_box = shift(patch.layout.AMRBox(), shifter);
            auto const ghostbox  = grow(shift_box, 1);
            for (std::size_t i = 0; i < pid; ++i)
                if (auto const overlap = ghostbox * patches[i].layout.AMRBox())
                    periodic_neighbours.emplace_back(&patches[i], *overlap, shifter);
            for (std::size_t i = pid + 1; i < patches.size(); ++i)
                if (auto const overlap = ghostbox * patches[i].layout.AMRBox())
                    periodic_neighbours.emplace_back(&patches[i], *overlap, shifter);
        }

        return periodic_neighbours;
    }


    TestGridLayout_t layout{cells * 3};
    Box_t const domainBox = layout.AMRBox();
    Box_t const ghostBox  = grow(layout.AMRBox(), 1);

    std::vector<Patch> patches{};
    std::vector<std::tuple<Patch*, Box_t>> neighbours{};
    std::vector<std::tuple<Patch*, Box_t, Point<int, 3>>> periodic_neighbours{};
};


template<typename Param>
struct ParticleArraySelectingTestSymmetric : public AParticleArraySelectingTest<Param>
{
    using Super        = AParticleArraySelectingTest<Param>;
    using GridLayout_t = Super::GridLayout_t;
    using Super::patches;

    ParticleArraySelectingTestSymmetric()
    {
        // 27 patches of cells**3 in 3**3 config
        patches.reserve(3 * 3 * 3);
        int const off = cells - 1;
        for (std::uint8_t i = 0; i < 3; ++i)
            for (std::uint8_t j = 0; j < 3; ++j)
                for (std::uint8_t k = 0; k < 3; ++k)
                {
                    int const cellx = i * cells;
                    int const celly = j * cells;
                    int const cellz = k * cells;
                    patches.emplace_back(Box<int, Super::dim>{
                        Point{cellx, celly, cellz}, Point{cellx + off, celly + off, cellz + off}});
                }

        assert(!any_overlaps_in(patches, [](auto const& patch) { return patch.layout.AMRBox(); }));
    }
};


template<typename ParticleArraySelectingTest_t>
auto run(ParticleArraySelectingTest_t& self)
{
    abort_if(self.periodic_neighbours_for(13).size());

    for (std::size_t pid = 0; pid < self.patches.size(); ++pid)
    {
        auto& dst = self.patches[pid];

        for (auto const& [src, overlap] : self.neighbours_for(pid))
            select_particles(src->domain, dst.ghost, overlap);

        for (auto const& [src, overlap, shift] : self.periodic_neighbours_for(pid))
            select_particles(src->domain, dst.ghost, overlap, shift);
    }

    auto const expected = pow(cells + 2, 3) - pow(cells, 3); // ghost box layer

    for (auto& patch : self.patches)
        ParticleArrayService::template sync<2, ParticleType::Ghost>(patch.ghost);

    for (std::size_t pid = 0; pid < self.patches.size(); ++pid)
    {
        auto const& patch = self.patches[pid];

        EXPECT_EQ(expected, patch.ghost.size());
    }
}

// clang-format off
using Permutations_t = testing::Types< // ! notice commas !

    TestParam<3, LayoutMode::AoSMapped, AllocatorMode::CPU>
    // ,TestParam<3, LayoutMode::AoSTS, AllocatorMode::CPU>

// PHARE_WITH_THRUST(
//     ,TestParam<3, LayoutMode::SoAPC, AllocatorMode::CPU>
// )

// PHARE_WITH_GPU( // GPU implies WITH_THRUST
//     ,TestParam<3, LayoutMode::AoSPC, AllocatorMode::GPU_UNIFIED>
//     ,TestParam<3, LayoutMode::SoAPC, AllocatorMode::GPU_UNIFIED>
// )

>;
// clang-format on

TYPED_TEST_SUITE(ParticleArraySelectingTestSymmetric, Permutations_t, );

TYPED_TEST(ParticleArraySelectingTestSymmetric, patch_ghost_copy)
{
    run(*this);
}



} // namespace PHARE::core


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
