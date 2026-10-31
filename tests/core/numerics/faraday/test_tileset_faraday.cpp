
#include "core/utilities/types.hpp"
#include "core/numerics/faraday/faraday.hpp"

#include "simulator/simulator_def.hpp"

#include "tests/core/data/gridlayout/test_gridlayout.hpp"
#include "tests/core/data/vecfield/test_vecfield_fixtures.hpp"
#include "tests/core/data/electromag/test_electromag_fixtures.hpp"

#include "gtest/gtest.h"

#include <cstddef>

#include "BS_thread_pool.hpp"

namespace PHARE::core
{

auto static const cells = get_env_as("PHARE_CELLS", std::uint32_t{4});
auto static const dt    = get_env_as("PHARE_TIMESTEP", double{.001});

static ::BS::thread_pool pool{4};

template<typename Patch>
void ref_do(Patch& patch)
{
    using GridLayout_t = Patch::GridLayout_t;

    Faraday<GridLayout_t>{patch.layout}(patch.em.B, patch.em.E, patch.emNew.B, dt);
}

template<typename Patch>
void cmp_do(Patch& patch)
{
    using GridLayout_t = Patch::GridLayout_t;

    auto const n_tiles = patch.em.B[0]().size();
    for (std::uint16_t ti = 0; ti < n_tiles; ++ti)
    {
        auto const B = tile_at(patch.em.B, ti);
        auto const E = tile_at(patch.em.E, ti);
        auto BNew    = tile_at(patch.emNew.B, ti);

        pool.detach_task([=, layout = patch.em.B[0][ti].layout()]() mutable {
            Faraday<GridLayout_t>{layout}(B, E, BNew, dt);
        });
    }
    pool.wait();
}



template<auto alloc_mode, typename GridLayout_t, typename R, typename C>
void compare(GridLayout_t const& layout, R& ref, C& cmp)
{
    using enum LayoutMode;
    using enum AllocatorMode;
    auto constexpr static n_components = 3;

    double diff = 1e-15;

    auto const& patch_box = layout.AMRBox();

    for (std::size_t c = 0; c < n_components; ++c)
    {
        auto const eq = compare_fields(ref.emNew.B[c], cmp.emNew.B[c]);
        PHARE_LOG_LINE_SS(eq.why());
        EXPECT_TRUE(eq) << "Failure for B New: " << eq.why();
    }
}


template<std::size_t _dim, auto _layout_mode, auto _alloc_mode>
struct TestParam
{
    static_assert(all_are<LayoutMode>(_layout_mode));
    static_assert(all_are<AllocatorMode>(_alloc_mode));

    auto constexpr static dim         = _dim;
    auto constexpr static layout_mode = _layout_mode;
    auto constexpr static alloc_mode  = _alloc_mode;
};


template<typename Param>
struct FaradayTileTest : public ::testing::Test
{
    struct TileSetPatch
    {
        auto constexpr static dim         = Param::dim;
        auto constexpr static interp      = 1;
        auto constexpr static layout_mode = Param::layout_mode;
        auto constexpr static alloc_mode  = Param::alloc_mode;
        auto constexpr static opts        = SimOpts{.dimension    = dim,
                                                    .interp_order = interp,
                                                    .layout_mode  = layout_mode,
                                                    .alloc_mode   = alloc_mode};

        using PhareTypes                 = PHARE_Types<opts>;
        using GridLayout_t               = PhareTypes::Hybrid::GridLayout_t;
        using Field_t                    = PhareTypes::Hybrid::Grid_t;
        using Hybrid_t                   = PhareTypes::Hybrid;
        auto constexpr static field_opts = TensorFieldOptions<Hybrid_t>{};
        using UsableElectromag_t         = UsableElectromag<field_opts>;
        using UsableVecField_t           = UsableVecField<field_opts>;


        TileSetPatch(GridLayout_t const& layout_)
            : layout{layout_}
        {
        }

        GridLayout_t layout;
        UsableElectromag_t em{layout}, emNew{layout};
    };

    struct ContiguousPatch
    {
        // reference (untiled) is always run on CPU, compared against the tiled/possibly-GPU patch
        auto constexpr static dim    = Param::dim;
        auto constexpr static interp = 1;
        auto constexpr static opts   = SimOpts{.dimension = dim, .interp_order = interp};

        using PhareTypes                 = PHARE_Types<opts>;
        using GridLayout_t               = PhareTypes::Hybrid::GridLayout_t;
        using Hybrid_t                   = PhareTypes::Hybrid;
        auto constexpr static field_opts = TensorFieldOptions<Hybrid_t>{};
        using UsableElectromag_t         = UsableElectromag<field_opts>;
        using Electromag_t               = UsableElectromag_t::Super;
        using Grid_t                     = PhareTypes::Hybrid::Grid_t;
        using UsableVecField_t           = UsableVecField<field_opts>;

        ContiguousPatch(GridLayout_t const& layout_)
            : layout{layout_}
        {
        }

        GridLayout_t layout;
        UsableElectromag_t em{layout}, emNew{layout};
    };

    using RefPatch = ContiguousPatch;
    using CmpPatch = TileSetPatch;

    FaradayTileTest() {}

    void init_check() const
    {
        bool static constexpr init_value_check = 0;
        if constexpr (init_value_check)
        {
            // auto ref = make_ions<RefParticleArray_t>(layout);
            // auto cmp = from_ions<CmpParticleArray_t>(layout, *ref);
            // compare(*this->layout, *ref, *cmp);
        }
    }

    void do_compare() const { compare<CmpPatch::alloc_mode>(*cmp_layout, ref_patch, cmp_patch); }

    TestGridLayout<typename RefPatch::GridLayout_t> const ref_layout{cells};
    TestGridLayout<typename CmpPatch::GridLayout_t> const cmp_layout{cells};
    RefPatch ref_patch{ref_layout};
    CmpPatch cmp_patch{cmp_layout};
};

// clang-format off
using Permutations_t = testing::Types<
     TestParam<3, LayoutMode::AoSPCTS, AllocatorMode::CPU>
>;
// clang-format on

TYPED_TEST_SUITE(FaradayTileTest, Permutations_t, );

template<typename FaradayTileTest_t>
auto run(FaradayTileTest_t& self)
{
    ref_do(self.ref_patch);
    cmp_do(self.cmp_patch);
    self.do_compare();
}

TYPED_TEST(FaradayTileTest, dispatch)
{
    run(*this);
}

} // namespace PHARE::core

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
