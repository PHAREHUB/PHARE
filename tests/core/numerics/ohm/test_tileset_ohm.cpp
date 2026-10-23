//  tests/core/numerics/ohm/test_tileset_ohm.cpp
//
// #define PHARE_UNDEF_ASSERT

#include "core/logger.hpp"    // scope timing
#include "core/def/phlop.hpp" // scope timing

#include "core/def/phare_config.hpp"

#include "core/utilities/types.hpp"
#include "core/numerics/ohm/ohm.hpp"
#include "core/data/grid/grid_tiles.hpp"
#include "core/data/particles/particle_array_def.hpp"
#include "core/models/quantities/hybrid_quantities.hpp"

#include "amr/solvers/solver_hybrid_field_evolvers.hpp"

#include "simulator/simulator_def.hpp"

#include "tests/core/data/field/test_field_fixtures.hpp"
#include "tests/core/data/gridlayout/test_gridlayout.hpp"
#include "tests/core/data/vecfield/test_vecfield_fixtures.hpp"
#include "tests/core/data/electromag/test_electromag_fixtures.hpp"

#include "gtest/gtest.h"

#include <cstddef>


namespace PHARE::core
{
double static constexpr eta = 1;
double static constexpr nu  = 1;
auto static const cells     = get_env_as("PHARE_CELLS", std::uint32_t{22});

template<typename Patch>
void ref_do(Patch& patch)
{
    using GridLayout_t = Patch::GridLayout_t;
    Ohm<GridLayout_t>{{eta, nu}, patch.layout}( //
        patch.n, patch.V, patch.P, patch.em.B, patch.J, patch.emNew.E);
}

template<typename Patch>
void cmp_do(Patch& patch)
{
    using GridLayout_t = Patch::GridLayout_t;
    using Field_vt     = Patch::Field_t::View::value_type;
    using VecField_vt  = basic::TensorField<Field_vt, 1>;

    auto constexpr tile_accessor = [](auto& ts, auto const ti) { return ts[ti]; };
    solver::OhmSingleTransformer{{eta, nu}}(patch.layout, *patch.n, *patch.V, *patch.P,
                                            *patch.em.B, *patch.J, *patch.emNew.E);

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



    using RefPatch = ContiguousPatch;
    using CmpPatch = TileSetPatch;

    OhmTileTest() {}

    void do_compare() const { compare<CmpPatch::alloc_mode>(*cmp_layout, ref_patch, cmp_patch); }

    TestGridLayout<typename RefPatch::GridLayout_t> const ref_layout{cells};
    TestGridLayout<typename CmpPatch::GridLayout_t> const cmp_layout{cells};
    RefPatch ref_patch{ref_layout};
    CmpPatch cmp_patch{cmp_layout};
};

// clang-format off
using Permutations_t = testing::Types< // ! notice commas !

     TestParam<3, LayoutMode::AoSPCTS, AllocatorMode::CPU>
PHARE_WITH_GPU(
    ,TestParam<3, LayoutMode::AoSPCTS, AllocatorMode::GPU_UNIFIED>
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
