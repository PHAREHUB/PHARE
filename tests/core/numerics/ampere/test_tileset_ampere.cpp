
#include "core/utilities/types.hpp"
#include "core/utilities/thread_pool.hpp"
#include "core/numerics/ampere/ampere.hpp"
#include "core/models/quantities/hybrid_quantities.hpp"

#include "simulator/simulator_def.hpp"

#include "tests/core/data/gridlayout/test_gridlayout.hpp"
#include "tests/core/data/vecfield/test_vecfield_fixtures.hpp"
#include "tests/core/data/electromag/test_electromag_fixtures.hpp"
#include "tests/core/data/tensorfield/test_tensorfield_fixtures.hpp"

#include "gtest/gtest.h"

#include <cstddef>

namespace PHARE::core
{
auto static const cells = get_env_as("PHARE_CELLS", std::uint32_t{4});

bool static const premain = []() {
    ThreadPool::threads_per_pool = 4;
    return true;
}();

auto& pool = *ThreadPool::INSTANCE().thread_pools[0];


template<typename Patch>
void ref_do(Patch& patch)
{
    using GridLayout_t = Patch::GridLayout_t;
    Ampere<GridLayout_t>{patch.layout}(patch.em.B, patch.J);
}

template<typename Patch>
void cmp_do(Patch& patch)
{
    using GridLayout_t = Patch::GridLayout_t;

    auto const n_tiles = patch.J[0]().size();
    for (std::uint16_t ti = 0; ti < n_tiles; ++ti)
    {
        auto const B = tile_at(patch.em.B, ti);
        auto J       = tile_at(patch.J, ti);
        pool.detach_task([=, layout = patch.em.B[0][ti].layout()]() mutable {
            Ampere<GridLayout_t>{layout}(B, J);
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

    double diff           = 1e-15;
    auto const& patch_box = layout.AMRBox();
    for (std::size_t c = 0; c < n_components; ++c)
    {
        {
            auto const eq = compare_fields(ref.em.B[c], cmp.em.B[c]);
            PHARE_LOG_LINE_SS(eq.why());
            EXPECT_TRUE(eq) << "Failure for B New: " << eq.why();
        }
        {
            auto const eq = compare_fields(ref.em.E[c], cmp.em.E[c]);
            PHARE_LOG_LINE_SS(eq.why());
            EXPECT_TRUE(eq) << "Failure for E New: " << eq.why();
        }
        auto const eq = compare_fields(ref.J[c], cmp.J[c]);
        PHARE_LOG_LINE_SS(eq.why());
        EXPECT_TRUE(eq) << "Failure for J: " << eq.why();
    }
    PHARE_LOG_LINE_SS(*ref.J);
    PHARE_LOG_LINE_SS(*cmp.J);
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
struct AmpereTileTest : public ::testing::Test
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
            PHARE_LOG_LINE_SS(*em.B);
            PHARE_LOG_LINE_SS(*em.E);
        }

        GridLayout_t layout;
        UsableElectromag_t em{layout};
        UsableVecField_t J{"J", layout, HybridQuantity::Vector::J, 0};
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
        UsableElectromag_t em{layout};
        UsableVecField_t J{"J", layout, HybridQuantity::Vector::J, 0};
    };

    using RefPatch = ContiguousPatch;
    using CmpPatch = TileSetPatch;

    AmpereTileTest() {}

    void do_compare() const { compare<CmpPatch::alloc_mode>(*cmp_layout, ref_patch, cmp_patch); }

    TestGridLayout<typename RefPatch::GridLayout_t> const ref_layout{cells};
    TestGridLayout<typename CmpPatch::GridLayout_t> const cmp_layout{cells};
    RefPatch ref_patch{ref_layout};
    CmpPatch cmp_patch{cmp_layout};
};

// clang-format off
using Permutations_t = testing::Types< // ! notice commas !
     TestParam<3, LayoutMode::AoSPCTS, AllocatorMode::CPU>
>;
// clang-format on

TYPED_TEST_SUITE(AmpereTileTest, Permutations_t, );

template<typename AmpereTileTest_t>
auto run(AmpereTileTest_t& self)
{
    ref_do(self.ref_patch);
    cmp_do(self.cmp_patch);
    self.do_compare();
}

TYPED_TEST(AmpereTileTest, dispatch)
{
    run(*this);
}

} // namespace PHARE::core

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
