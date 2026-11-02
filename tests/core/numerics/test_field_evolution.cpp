//  tests/core/numerics/test_field_evolution.cpp
//
//

#include "core/def/phare_config.hpp"

#include "core/utilities/types.hpp"
#include "core/numerics/ohm/ohm.hpp"
#include "core/data/grid/grid_tiles.hpp"
#include "core/numerics/ampere/ampere.hpp"
#include "core/numerics/faraday/faraday.hpp"
#include "core/models/quantities/hybrid_quantities.hpp"


#include "simulator/simulator_def.hpp"

#include "tests/core/data/gridlayout/test_gridlayout.hpp"
#include "tests/core/data/vecfield/test_vecfield_fixtures.hpp"
#include "tests/core/data/electromag/test_electromag_fixtures.hpp"

#include "gtest/gtest.h"

#include <cstddef>

#include "BS_thread_pool.hpp"

namespace PHARE::core
{
// COMPILE TIME SWITCHES
double static constexpr eta = 1;
double static constexpr nu  = 1;


// RUNTIME ENV VAR OVERRIDES
auto static const cells     = get_env_as("PHARE_CELLS", std::uint32_t{22});
auto static const n_patches = get_env_as("PHARE_PATCHES", std::size_t{1});
auto static const dt        = get_env_as("PHARE_TIMESTEP", double{.001});
auto static const do_cmp    = get_env_as("PHARE_COMPARE", std::size_t{1});
static ::BS::thread_pool pool{4};


template<typename Patch>
void ref_do(Patch& patch)
{
    FaradaySingleTransformer{}(patch.layout, *patch.em.B, *patch.em.E, *patch.emNew.B, dt);

    AmpereSingleTransformer{}(patch.layout, *patch.emNew.B, *patch.J);

    OhmSingleTransformer{{eta, nu}}(patch.layout, *patch.n, *patch.V, *patch.P, *patch.emNew.B,
                                    *patch.J, *patch.emNew.E);
}

template<typename Patch>
void cmp_do(Patch& patch)
{
    FaradaySingleTransformer{}(patch.layout, *patch.em.B, *patch.em.E, *patch.emNew.B, dt);

    AmpereSingleTransformer{}(patch.layout, *patch.emNew.B, *patch.J);

    OhmSingleTransformer{{eta, nu}}(patch.layout, *patch.n, *patch.V, *patch.P, *patch.emNew.B,
                                    *patch.J, *patch.emNew.E);
}



template<auto alloc_mode, typename GridLayout_t, typename R, typename C>
void compare(GridLayout_t const& layout, R& ref, C& cmp)
{
    using enum LayoutMode;
    using enum AllocatorMode;
    auto constexpr static n_components = 3;

    double diff = 1e-15;

    // if constexpr (alloc_mode == GPU_UNIFIED)
    //     diff *= 1e3; // atomic no order guaranteed

    auto const& patch_box = layout.AMRBox();

    for (std::size_t c = 0; c < n_components; ++c)

        for (std::size_t c = 0; c < n_components; ++c)
        {
            auto const eq = compare_fields(ref.em.B[c], cmp.em.B[c]);
            PHARE_LOG_LINE_SS(eq.why());
            EXPECT_TRUE(eq) << "Failure for B New: " << eq.why();
        }
    for (std::size_t c = 0; c < n_components; ++c)
    {
        auto const eq = compare_fields(ref.em.E[c], cmp.em.E[c]);
        PHARE_LOG_LINE_SS(eq.why());
        EXPECT_TRUE(eq) << "Failure for E New: " << eq.why();
    }
    for (std::size_t c = 0; c < n_components; ++c)
    {
        auto const eq = compare_fields(ref.J[c], cmp.J[c]);
        PHARE_LOG_LINE_SS(eq.why());
        EXPECT_TRUE(eq) << "Failure for J: " << eq.why();
    }
    // PHARE_LOG_LINE_SS(*ref.J);
    // PHARE_LOG_LINE_SS(*cmp.J);
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

        using PhareTypes   = PHARE_Types<opts>;
        using GridLayout_t = PhareTypes::Hybrid::GridLayout_t;
        using Field_t      = PhareTypes::Hybrid::Grid_t;
        static_assert(is_field_tile_set_v<Field_t>);
        using Hybrid_t                   = PhareTypes::Hybrid;
        auto constexpr static field_opts = TensorFieldOptions<Hybrid_t>{};
        using UsableElectromag_t         = UsableElectromag<field_opts>;
        using UsableVecField_t           = UsableVecField<field_opts>;


        TileSetPatch(GridLayout_t const& layout_)
            : layout{layout_}
        {
            V.fill(.001);
            n.fill(.001);
            P.fill(.001);
        }

        GridLayout_t layout;
        UsableElectromag_t em{layout}, emNew{layout};
        UsableVecField_t J{"J", layout, HybridQuantity::Vector::J};
        UsableVecField_t V{"V", layout, HybridQuantity::Vector::V};
        Field_t n{"rho", layout, HybridQuantity::Scalar::rho};
        Field_t P{"P", layout, HybridQuantity::Scalar::P};
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
            V.fill(.001);
            n.fill(.001);
            P.fill(.001);
        }

        GridLayout_t layout;
        UsableElectromag_t em{layout}, emNew{layout};
        UsableVecField_t J{"J", layout, HybridQuantity::Vector::J};
        UsableVecField_t V{"V", layout, HybridQuantity::Vector::V};
        Grid_t n{"rho", layout, HybridQuantity::Scalar::rho};
        Grid_t P{"P", layout, HybridQuantity::Scalar::P};
    };

    using RefPatch = ContiguousPatch;
    using CmpPatch = TileSetPatch;


    FaradayTileTest() {}


    void do_compare() const { compare<CmpPatch::alloc_mode>(*cmp_layout, ref_patch, cmp_patch); }

    TestGridLayout<typename RefPatch::GridLayout_t> const ref_layout{cells};
    TestGridLayout<typename CmpPatch::GridLayout_t> const cmp_layout{cells};
    RefPatch ref_patch{ref_layout};
    CmpPatch cmp_patch{cmp_layout};
};

// clang-format off
using Permutations_t = testing::Types< // ! notice commas !

PHARE_WITH_MKN_GPU(
    TestParam<1, LayoutMode::AoSTS, AllocatorMode::CPU>
   ,TestParam<2, LayoutMode::AoSTS, AllocatorMode::CPU>
   // ,TestParam<3, LayoutMode::AoSTS, AllocatorMode::CPU>

)

PHARE_WITH_GPU(
   ,TestParam<1, LayoutMode::AoSTS, AllocatorMode::GPU_UNIFIED>
   ,TestParam<2, LayoutMode::AoSTS, AllocatorMode::GPU_UNIFIED>
   // ,TestParam<3, LayoutMode::AoSTS, AllocatorMode::GPU_UNIFIED>
)

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
    // assert(phlop::ScopeTimerMan::INSTANCE().active);
    ::testing::InitGoogleTest(&argc, argv);
    auto r = RUN_ALL_TESTS();
    PHARE_WITH_PHLOP(phlop::threaded::ScopeTimerMan::reset());
    return r;
}
