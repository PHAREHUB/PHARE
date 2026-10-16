
#include "phare_core.hpp"
#include "core/utilities/types.hpp"
#include "core/utilities/box/box.hpp"
#include "core/data/particles/particle_array_def.hpp"
#include "core/models/quantities/hybrid_quantities.hpp"
#include "core/data/tiles/tile_set_traversal.hpp"
#include "core/numerics/ion_updater/ion_updater_def.hpp"

#include "simulator/simulator_def.hpp"

#include "tests/core/data/gridlayout/test_gridlayout.hpp"

#include "gtest/gtest.h"

#include <unordered_map>

using namespace PHARE;
using namespace PHARE::core;

std::size_t constexpr static interp = 1;
std::size_t constexpr static cells  = 40;


template<std::size_t _dim, auto lm, auto am>
struct TestParam
{
    auto constexpr static dim = _dim;

    constexpr static auto opts
        = SimOpts{.dimension = dim, .interp_order = interp, .layout_mode = lm, .alloc_mode = am};

    using PhareTypes = core::PHARE_Types<opts>;

    using Box_t            = PHARE::core::Box<int, dim>;
    using GridLayout_t     = PhareTypes::Hybrid::GridLayout_t;
    using TestGridLayout_t = TestGridLayout<GridLayout_t>;
};



template<typename TestParam>
struct Patch
{
    auto constexpr static dim = TestParam::dim;

    using PhareTypes       = TestParam::PhareTypes;
    using GridLayout_t     = TestParam::GridLayout_t;
    using TestGridLayout_t = TestParam::TestGridLayout_t;
    using Box_t            = TestParam::Box_t;

    using Field_t = PhareTypes::Hybrid::Field_t;
    using Grid_t  = PhareTypes::Hybrid::Grid_t;

    Patch(Box_t const& box)
        : layout{box}
    {
    }
    Patch(GridLayout_t const& _layout)
        : layout{_layout}
    {
    }

    TestGridLayout_t const layout;
    Grid_t rho{"rho", layout, HybridQuantity::Scalar::rho};
    Field_t& rho_v = *rho;
    Grid_t Ex{"Ex", layout, HybridQuantity::Scalar::Ex}; // staggered - different ghost/field box
    Field_t& Ex_v = *Ex;
};

namespace PHARE::core
{

// quantity-independent tile adjacency for one patch: every field sharing this patch's
// GridLayout splits the same domain box via the same deterministic Tiler, so tile i is the
// same box across fields even though its ghost_box/field_box differ per physical quantity
// (Yee staggering). Built once (via the links system), reused by every field's overlap cache
// below so only the first field pays for the O(N*26) neighbour discovery.
struct TileNeighbourTopology
{
    template<typename FieldTileSet_t>
    TileNeighbourTopology(FieldTileSet_t& field)
        : neighbours(field().size())
    {
        std::unordered_map<void const*, std::size_t> index_of;
        for (std::size_t i = 0; i < field().size(); ++i)
            index_of[&field()[i]] = i;

        visit_tile_neighbours(field, [&](auto& t0, auto& t1) {
            neighbours[index_of.at(&t0)].push_back(index_of.at(&t1));
        });
    }

    std::vector<std::vector<std::size_t>> neighbours;
};


// owns the precomputed inner-ghost overlap cache for a single field tile set. The pair
// discovery itself is not repeated here - it reuses a TileNeighbourTopology built once for
// the patch - only this field's own ghost_box()/field_box() (from its own tiles' GridLayout +
// physical quantity) are (re)computed. sync_inner_ghosts() then just replays the cache.
// Subclasses UpdaterSelectionBoxing purely to reuse its per-patch layout/box bookkeeping.
template<typename FieldTileSet_t, typename GridLayout_t>
struct UpdaterTileOverlapSelectionBoxing : public UpdaterSelectionBoxing<GridLayout_t>
{
    using Super                     = UpdaterSelectionBoxing<GridLayout_t>;
    using Field_t                   = FieldTileSet_t::field_type; // per-tile field
    auto constexpr static dimension = GridLayout_t::dimension;

    struct Overlap
    {
        Overlap(Field_t* s, Box<std::uint32_t, dimension> const& d,
                Box<std::uint32_t, dimension> const& sb)
            : src{s}
            , dst_lcl{d}
            , src_lcl{sb}
        {
        }

        Field_t* src;
        Box<std::uint32_t, dimension> dst_lcl, src_lcl;
    };

    UpdaterTileOverlapSelectionBoxing(GridLayout_t const& layout_, FieldTileSet_t& field,
                                      TileNeighbourTopology const& topology)
        : Super{layout_, {}}
        , overlaps(field().size())
    {
        for (std::size_t i = 0; i < field().size(); ++i)
        {
            auto& t0 = field()[i];
            for (auto const j : topology.neighbours[i])
            {
                auto& t1 = field()[j];
                if (auto const overlap = t0.ghost_box() * t1.field_box())
                    overlaps[i].emplace_back(
                        &t1(), as_unsigned(shift(*overlap, t0.ghost_box().lower * -1)),
                        as_unsigned(shift(*overlap, t1.ghost_box().lower * -1)));
            }
        }
    }

    void sync_inner_ghosts(FieldTileSet_t& field) const
    {
        for (std::size_t i = 0; i < field().size(); ++i)
            for (auto const& ov : overlaps[i])
                for (auto const& lix : ov.dst_lcl)
                    field()[i]()(lix) = (*ov.src)(lix - ov.dst_lcl.lower + ov.src_lcl.lower);
    }

    std::vector<std::vector<Overlap>> overlaps;
};


// builds one shared topology + one overlap cache per field, per patch - done once here
// (fixture construction) rather than inside a test body, so it's excluded from what a test's
// timing is meant to measure (repeated sync_inner_ghosts() calls).
template<typename TestParam>
void build_precomputed_overlaps(
    std::vector<Patch<TestParam>>& patches, std::vector<TileNeighbourTopology>& topologies,
    std::vector<UpdaterTileOverlapSelectionBoxing<typename Patch<TestParam>::Field_t,
                                                  typename TestParam::GridLayout_t>>& rho_boxings,
    std::vector<UpdaterTileOverlapSelectionBoxing<typename Patch<TestParam>::Field_t,
                                                  typename TestParam::GridLayout_t>>& ex_boxings)
{
    topologies.reserve(patches.size());
    rho_boxings.reserve(patches.size());
    ex_boxings.reserve(patches.size());
    for (auto& patch : patches)
    {
        topologies.emplace_back(patch.rho_v);
        rho_boxings.emplace_back(patch.layout, patch.rho_v, topologies.back());
        ex_boxings.emplace_back(patch.layout, patch.Ex_v, topologies.back());
    }
}

} // namespace PHARE::core

template<std::size_t dim, typename TestParam>
struct AGridFieldTest;

template<typename TestParam>
struct AGridFieldTest<1, TestParam>
{
    using Box_t   = TestParam::Box_t;
    using Field_t = typename Patch<TestParam>::Field_t;
    using Boxing_t
        = PHARE::core::UpdaterTileOverlapSelectionBoxing<Field_t, typename TestParam::GridLayout_t>;

    AGridFieldTest()
    {
        patches.reserve(3);
        auto const off = cells - 1;
        for (std::uint8_t i = 0; i < 3; ++i)
        {
            auto const cellx = i * cells;
            patches.emplace_back(Box_t{Point{cellx}, Point{cellx + off}});
        }
        PHARE::core::build_precomputed_overlaps(patches, topologies, rho_boxings, ex_boxings);
    }

    std::vector<Patch<TestParam>> patches;
    std::vector<PHARE::core::TileNeighbourTopology> topologies;
    std::vector<Boxing_t> rho_boxings, ex_boxings;
    Patch<TestParam> L1{Box_t{Point{3}, Point{12}}};
};


template<typename TestParam>
struct AGridFieldTest<2, TestParam>
{
    using Box_t   = TestParam::Box_t;
    using Field_t = typename Patch<TestParam>::Field_t;
    using Boxing_t
        = PHARE::core::UpdaterTileOverlapSelectionBoxing<Field_t, typename TestParam::GridLayout_t>;

    AGridFieldTest()
    {
        patches.reserve(3 * 3);
        auto const off = cells - 1;
        for (std::uint8_t i = 0; i < 3; ++i)
            for (std::uint8_t j = 0; j < 3; ++j)
            {
                auto const cellx = i * cells;
                auto const celly = j * cells;
                patches.emplace_back(Box_t{Point{cellx, celly}, Point{cellx + off, celly + off}});
            }
        PHARE::core::build_precomputed_overlaps(patches, topologies, rho_boxings, ex_boxings);
    }

    std::vector<Patch<TestParam>> patches;
    std::vector<PHARE::core::TileNeighbourTopology> topologies;
    std::vector<Boxing_t> rho_boxings, ex_boxings;
    Patch<TestParam> L1{Box_t{Point{3, 3}, Point{12, 12}}};
};


template<typename TestParam>
struct AGridFieldTest<3, TestParam>
{
    using Box_t   = TestParam::Box_t;
    using Field_t = typename Patch<TestParam>::Field_t;
    using Boxing_t
        = PHARE::core::UpdaterTileOverlapSelectionBoxing<Field_t, typename TestParam::GridLayout_t>;

    AGridFieldTest()
    {
        patches.reserve(3 * 3 * 3);
        auto const off = cells - 1;
        for (std::uint8_t i = 0; i < 3; ++i)
            for (std::uint8_t j = 0; j < 3; ++j)
                for (std::uint8_t k = 0; k < 3; ++k)
                {
                    auto const cellx = i * cells;
                    auto const celly = j * cells;
                    auto const cellz = k * cells;
                    patches.emplace_back(Box_t{Point{cellx, celly, cellz},
                                               Point{cellx + off, celly + off, cellz + off}});
                }
        PHARE::core::build_precomputed_overlaps(patches, topologies, rho_boxings, ex_boxings);
    }

    std::vector<Patch<TestParam>> patches;
    std::vector<PHARE::core::TileNeighbourTopology> topologies;
    std::vector<Boxing_t> rho_boxings, ex_boxings;
    Patch<TestParam> L1{Box_t{Point{3, 3, 3}, Point{12, 12, 12}}};
};


template<typename TestParam>
struct GridFieldTest : public ::testing::Test, public AGridFieldTest<TestParam::dim, TestParam>
{
    using Super        = AGridFieldTest<TestParam::dim, TestParam>;
    using GridLayout_t = TestParam::GridLayout_t;

    using Super::ex_boxings;
    using Super::L1;
    using Super::patches;
    using Super::rho_boxings;
    using Super::topologies;

    GridFieldTest() {}
};


// clang-format off
using ParticlesDatas = testing::Types< //

   TestParam<1, LayoutMode::AoSPCTS, AllocatorMode::CPU>,
   TestParam<2, LayoutMode::AoSPCTS, AllocatorMode::CPU>,
   TestParam<3, LayoutMode::AoSPCTS, AllocatorMode::CPU>

>;
// clang-format on

TYPED_TEST_SUITE(GridFieldTest, ParticlesDatas);


namespace PHARE::core
{

TYPED_TEST(GridFieldTest, test_compiles)
{
    auto& L1 = this->L1;

    EXPECT_EQ(L1.rho().data(), L1.rho_v().data());

    for (auto& tile : L1.rho())
    {
        Point const lower{tile.lower};
        EXPECT_EQ((*L1.rho).at(lower), L1.rho_v.at(lower));
    }
}


TYPED_TEST(GridFieldTest, test_ghost_sync_links)
{
    using GridLayout_t = TestFixture::GridLayout_t;

    auto& patches = this->patches;

    auto const verify_synced = [](auto& patch, auto& field, auto const qty) {
        auto const patch_gb = patch.layout.AMRGhostBoxFor(qty);
        auto const patch_pb = shrink(patch_gb, GridLayout_t::options.field_ghost_width);
        auto const g_layer  = patch_gb.remove(patch_pb);
        for (auto& tile : field())
        {
            std::size_t expected = tile.ghost_box().size();

            for (auto const& gl : g_layer)
                if (auto const glo = tile.ghost_box() * gl)
                    expected -= glo->size();

            EXPECT_EQ(sum(tile()), expected);
        }
    };

    for (auto& patch : patches)
    {
        patch.rho.fill(0);
        patch.Ex.fill(0);
    }

    for (auto& patch : patches)
        for (auto& field : {&patch.rho_v, &patch.Ex_v})
            for (auto& tile : (*field)())
                for (auto const& lix : tile.layout().domainBoxFor(tile))
                    tile()(lix) = 1;

    // same two fields as test_ghost_sync_precomputed_overlaps, so the two tests are comparable.
    for (auto& patch : patches)
    {
        patch.rho.sync_inner_ghosts();
        patch.Ex.sync_inner_ghosts();

        verify_synced(patch, patch.rho_v, HybridQuantity::Scalar::rho);
        verify_synced(patch, patch.Ex_v, HybridQuantity::Scalar::Ex);
    }
}


TYPED_TEST(GridFieldTest, test_ghost_sync_precomputed_overlaps)
{
    using GridLayout_t = TestFixture::GridLayout_t;

    auto& patches     = this->patches;
    auto& rho_boxings = this->rho_boxings;
    auto& ex_boxings  = this->ex_boxings;

    auto const verify_synced = [](auto& patch, auto& field, auto const qty) {
        auto const patch_gb = patch.layout.AMRGhostBoxFor(qty);
        auto const patch_pb = shrink(patch_gb, GridLayout_t::options.field_ghost_width);
        auto const g_layer  = patch_gb.remove(patch_pb);
        for (auto& tile : field())
        {
            std::size_t expected = tile.ghost_box().size();

            for (auto const& gl : g_layer)
                if (auto const glo = tile.ghost_box() * gl)
                    expected -= glo->size();

            EXPECT_EQ(sum(tile()), expected);
        }
    };

    for (auto& patch : patches)
    {
        patch.rho.fill(0);
        patch.Ex.fill(0);
    }

    for (auto& patch : patches)
        for (auto& field : {&patch.rho_v, &patch.Ex_v})
            for (auto& tile : (*field)())
                for (auto const& lix : tile.layout().domainBoxFor(tile))
                    tile()(lix) = 1;

    // topology + overlap caches were already built once in the fixture constructor, not here.
    for (std::size_t i = 0; i < patches.size(); ++i)
    {
        auto& patch = patches[i];

        rho_boxings[i].sync_inner_ghosts(patch.rho_v);
        ex_boxings[i].sync_inner_ghosts(patch.Ex_v);

        verify_synced(patch, patch.rho_v, HybridQuantity::Scalar::rho);
        verify_synced(patch, patch.Ex_v, HybridQuantity::Scalar::Ex);
    }
}

} // namespace PHARE::core


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
