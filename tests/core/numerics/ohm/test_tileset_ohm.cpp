//  tests/core/numerics/ohm/test_tileset_ohm.cpp
//
// #define PHARE_UNDEF_ASSERT
#define PHARE_SKIP_MPI_IN_CORE

#include "core/logger.hpp"    // scope timing
#include "core/def/phlop.hpp" // scope timing


#include "core/def/phare_config.hpp"

#include "core/utilities/types.hpp"
#include "core/numerics/ohm/ohm.hpp"
#include "core/data/grid/grid_tiles.hpp"
#include "core/hybrid/hybrid_quantities.hpp"
#include "core/data/particles/particle_array_def.hpp"

#include "tests/core/data/field/test_field_fixtures.hpp"
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
auto static const cells = get_env_as("PHARE_CELLS", std::uint32_t{22});

static ::BS::thread_pool pool{4};


template<typename Patch>
void ref_do(Patch& patch)
{
    PHARE_LOG_LINE_SS("ref_do");
    using GridLayout_t = Patch::GridLayout_t;

    Ohm<GridLayout_t>{{eta, nu}, patch.layout}( //
        patch.n, patch.V, patch.P, patch.em.B, patch.J, patch.emNew.E);
}



template<typename Patch>
void cmp_do(Patch& patch)
{
    PHARE_LOG_LINE_SS("cmp_do");
    using GridLayout_t = Patch::GridLayout_t;
    using Field_vt     = Patch::Field_t::View::value_type;
    using VecField_vt  = basic::TensorField<Field_vt, 1>;

    auto constexpr tile_accessor = [](auto& ts, auto const ti) { return ts[ti]; };

    OhmSingleTransformer{{eta, nu}}(patch.layout, *patch.n, *patch.V, *patch.P, *patch.em.B,
                                    *patch.J, *patch.emNew.E);

    // auto const n_tiles = patch.em.B[0]().size();
    // for (std::uint16_t ti = 0; ti < n_tiles; ++ti)
    // {
    //     auto const n = *patch.n[ti];
    //     auto const V = patch.V.template as<VecField_vt>(tile_accessor, ti);
    //     auto const P = *patch.P[ti];
    //     auto const B = patch.em.B.template as<VecField_vt>(tile_accessor, ti);
    //     auto const J = patch.J.template as<VecField_vt>(tile_accessor, ti);
    //     auto ENew    = patch.emNew.E.template as<VecField_vt>(tile_accessor, ti);
    //     ohm(patch.em.E[0][ti].layout(), n, V, P, B, J, ENew);
    // }

    // pool.detach_task([&]() {
    //     auto const n_tiles = patch.em.B[0]().size();
    //     for (std::uint16_t ti = 0; ti < n_tiles; ++ti)
    //     {
    //         auto const n = *patch.n[ti];
    //         auto const P = *patch.P[ti];
    //         auto const V = patch.V.template as<VecField_vt>(tile_accessor, ti);
    //         auto const J = patch.J.template as<VecField_vt>(tile_accessor, ti);
    //         auto const B = patch.em.B.template as<VecField_vt>(tile_accessor, ti);
    //         auto ENew    = patch.emNew.E.template as<VecField_vt>(tile_accessor, ti);

    //         pool.detach_task([=, layout = patch.em.E[0][ti].layout()]() mutable {
    //             Ohm<GridLayout_t>{layout, eta, nu}(n, V, P, B, J, ENew);
    //         });
    //     }
    // });

    // pool.wait();
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
    {
        auto const eq = compare_fields(ref.emNew.E[c], cmp.emNew.E[c]);
        PHARE_LOG_LINE_SS(eq.why());
        EXPECT_TRUE(eq) << "Failure for E New: " << eq.why();
    }

    for (std::size_t c = 0; c < n_components; ++c)
    {
        auto const sum = sum_field(cmp.emNew.E[c]);
        PHARE_LOG_LINE_SS(sum);
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
struct OhmTileTest : public ::testing::Test
{
    auto constexpr static c_order     = true;
    auto constexpr static dim         = Param::dim;
    auto constexpr static interp      = 1;
    auto constexpr static layout_mode = Param::layout_mode;
    auto constexpr static alloc_mode  = Param::alloc_mode;

    auto constexpr static sim_opts = SimOpts{.dimension = dim, .interp_order = interp,
                                             .layout_mode = layout_mode, .alloc_mode = alloc_mode};

    template<typename T, auto am>
    using nd_array_t       = NdArrayVector<dim, T, c_order, am>;
    using PhareTypes       = PHARE_Types<sim_opts>;
    using GridLayout_t     = PhareTypes::GridLayout_t;
    using TestGridLayout_t = TestGridLayout<GridLayout_t>;

    template<auto am>
    struct TileSetPatch
    {
        using GridLayout_t       = OhmTileTest<Param>::GridLayout_t;
        using Field_t            = PhareTypes::Grid_t;
        using UsableElectromag_t = UsableElectromag<GridLayout_t, am, LayoutMode::AoSTS>;
        using UsableVecField_t   = UsableVecField<GridLayout_t, am, LayoutMode::AoSTS>;


        TileSetPatch(GridLayout_t const& layout_)
            : layout{layout_}
        {
            J.fill(.001);
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

    template<auto am = AllocatorMode::CPU>
    struct ContiguousPatch
    {
        using GridLayout_t       = OhmTileTest<Param>::GridLayout_t;
        using UsableElectromag_t = UsableElectromag<GridLayout_t, am>;
        using Electromag_t       = UsableElectromag_t::Super;
        using Grid_t             = Grid<nd_array_t<double, am>, HybridQuantity::Scalar>;
        using UsableVecField_t   = UsableVecField<GridLayout_t, am>;

        ContiguousPatch(GridLayout_t const& layout_)
            : layout{layout_}
        {
            J.fill(.001);
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



    OhmTileTest() {}

    void do_compare() const { compare<alloc_mode>(*layout, ref_patch, cmp_patch); }

    TestGridLayout_t const layout{cells};
    ContiguousPatch<> ref_patch{layout};
    TileSetPatch<alloc_mode> cmp_patch{layout};
};

// clang-format off
using Permutations_t = testing::Types< // ! notice commas !

     TestParam<3, LayoutMode::AoSTS, AllocatorMode::CPU>
    // ,TestParam<3, LayoutMode::SoATS, AllocatorMode::CPU, 2>

PHARE_WITH_GPU(
    ,TestParam<3, LayoutMode::AoSTS, AllocatorMode::GPU_UNIFIED>
    // ,TestParam<3, LayoutMode::SoATS, AllocatorMode::GPU_UNIFIED, 2>
)

>;
// clang-format on

TYPED_TEST_SUITE(OhmTileTest, Permutations_t, );

template<typename OhmTileTest_t>
auto run(OhmTileTest_t& self)
{
    ref_do(self.ref_patch);
    cmp_do(self.cmp_patch);
    self.do_compare();
}


TYPED_TEST(OhmTileTest, dispatch)
{
    run(*this);
}



} // namespace PHARE::core


int main(int argc, char** argv)
{
    PHARE_WITH_PHLOP(phlop::threaded::ScopeTimerMan::INSTANCE().init();)
    ::testing::InitGoogleTest(&argc, argv);
    auto const ret = RUN_ALL_TESTS();
    PHARE_WITH_PHLOP(phlop::threaded::ScopeTimerMan::reset();)
    return ret;
}
