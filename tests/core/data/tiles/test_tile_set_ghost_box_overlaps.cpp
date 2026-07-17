// deal with ghost box of tiles to find where tile ghost boxes overlap


#include "phare_core.hpp"
#include "core/def/phare_config.hpp"

#include "core/utilities/types.hpp"
#include "core/utilities/box/box.hpp"
#include "core/hybrid/hybrid_quantities.hpp"
#include "core/data/tiles/tile_set_overlaps.hpp"

#include "tests/core/data/particles/test_particles.hpp"
#include "tests/core/data/gridlayout/test_gridlayout.hpp"
#include "tests/core/data/vecfield/test_vecfield_fixtures.hpp"

#include "gtest/gtest.h"



namespace PHARE::core
{

auto static const cells = get_env_as("PHARE_CELLS", std::uint32_t{22});
auto static const ppc   = get_env_as("PHARE_PPC", std::size_t{10});


template<std::size_t _dim, std::size_t _interp, auto lm, auto am>
struct TestParam
{
    static_assert(all_are<LayoutMode>(lm));
    static_assert(all_are<AllocatorMode>(am));

    auto constexpr static dim         = _dim;
    auto constexpr static interp      = _interp;
    auto constexpr static layout_mode = lm;
    auto constexpr static alloc_mode  = am;
    auto static constexpr opts        = SimOpts{.dimension = dim, .interp_order = interp,
                                          .layout_mode = layout_mode, .alloc_mode = alloc_mode};

    using PhareTypes   = PHARE::core::PHARE_Types<opts>;
    using Box_t        = Box<int, dim>;
    using GridLayout_t = TestGridLayout<typename PhareTypes::GridLayout_t>;
    using ParticleArray_t
        = ParticleArray<ParticleArrayOptions{dim, layout_mode, StorageMode::VECTOR, alloc_mode}>;
};


template<typename TestParam>
struct Patch
{
    auto constexpr static dim = TestParam::dim;
    using GridLayout_t        = TestParam::GridLayout_t;
    using ParticleArray_t     = TestParam::ParticleArray_t;

    Patch(GridLayout_t const& _layout)
        : layout{_layout}
    {
    }


    GridLayout_t layout;
    ParticleArray_t domain = make_particles<ParticleArray_t>(*layout);
    ParticleArray_t ghost  = make_particles<ParticleArray_t>(*layout);
};


template<typename TileSetTest_t>
void setup(TileSetTest_t& test)
{
    auto constexpr static dim = TileSetTest_t::dim;
    using Box_t               = Box<int, dim>;

    auto& patches  = test.patches;
    auto const off = cells - 1;
    if constexpr (dim == 1)
    {
        patches.reserve(3);
        for (std::uint8_t i = 0; i < 3; ++i)
        {
            auto const cellx = i * cells;
            patches.emplace_back(Box_t{Point{cellx}, Point{cellx + off}});
        }
    }
    if constexpr (dim == 2)
    {
        patches.reserve(3 * 3);
        for (std::uint8_t i = 0; i < 3; ++i)
            for (std::uint8_t j = 0; j < 3; ++j)
            {
                auto const cellx = i * cells;
                auto const celly = j * cells;
                patches.emplace_back(Box_t{Point{cellx, celly}, Point{cellx + off, celly + off}});
            }
    }
    if constexpr (dim == 3)
    {
        patches.reserve(3 * 3 * 3);
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
    }
}


template<typename TestParam>
struct TileSetOverlapsTest : public ::testing::Test
{
    auto constexpr static dim = TestParam::dim;
    using Patch_t             = Patch<TestParam>;
    using ParticleArray_t     = TestParam::ParticleArray_t;
    using GridLayout_t        = TestParam::GridLayout_t;

    TileSetOverlapsTest()
    {
        setup(*this);
        assert(!any_overlaps_in(patches, [](auto const& patch) { return patch.layout.AMRBox(); }));
    }

    auto static overlap(Patch_t const& src, Patch_t& dst) { return Patch_t::overlap(src, dst); }
    auto static overlap(Patch_t const& src, Patch_t& dst, auto const shift)
    {
        return Patch_t::overlap(src, dst, shift);
    }

    auto mid_pid() const { return int(((patches.size() + 1) / 2) - 1); }

    auto static lcl_cell(auto const& box, auto const& icell)
    {
        return (icell - box.lower).as_unsigned();
    }

    std::vector<Patch_t> patches;
};


template<typename TileSetTest_t>
auto make_nd_span_set_from(TileSetTest_t& self, int const pid)
{
    auto& patch  = self.patches[pid];
    auto spanset = make_particle_nd_span_set_from(patch.ghost(), patch.layout);

    check_tile_span_set(spanset);

    return spanset;
}


template<typename TileSetTest_t>
auto patch_ghost_copy(TileSetTest_t& self)
{
    using ParticleArray_t           = TileSetTest_t::ParticleArray_t;
    using Tile_t                    = ParticleArray_t::VecTile;
    using GridLayout_t              = TileSetTest_t::GridLayout_t;
    auto constexpr static nbrGhosts = GridLayout_t::nbrParticleGhosts();

    auto const pid  = 0;
    auto& p0        = self.patches[pid]; // origin patch
    auto const p0b  = p0.layout.AMRBox();
    auto const p0gb = grow(p0b, nbrGhosts);
    auto ss0        = make_nd_span_set_from(self, pid);

    // add some ghosts - or particles that left domain
    for (auto const& box : p0gb.remove(p0b))
        for (auto const& bix : box) // pick first tile for ghost cell
            add_particles_in((**unique_tiles_for(ss0, self.lcl_cell(p0gb, bix)).begin())(),
                             asBox(bix), ppc);

    auto& pM           = self.patches[self.mid_pid()]; // middle patch
    auto const pMb     = pM.layout.AMRBox();
    auto const overlap = *(p0gb * pMb);
    onBox(ss0, asBox(overlap.lower.as_unsigned()), [&](auto& tile) {
        for (auto const& p : tile())
            if (isIn(p, pMb))
                pM.domain.emplace_back(p);
    });

    EXPECT_EQ(pM.domain.size(), ppc * overlap.size()); // 1 cell overlap in interp 1!
    std::size_t count = 0;
    for (auto const& tile : pM.domain())
        for (auto const& p : tile())
        {
            ++count;
            EXPECT_TRUE(isIn(p, pMb));
        }

    EXPECT_EQ(count, ppc * overlap.size());
}



template<typename TileSetTest_t>
auto patch_tile_field(TileSetTest_t& self)
{
    using GridLayout_t     = TileSetTest_t::GridLayout_t;
    using UsableVecField_t = UsableVecField<GridLayout_t, AllocatorMode::CPU, LayoutMode::AoSTS>;
    auto constexpr static nbrGhosts = GridLayout_t::nbrGhosts();
    auto constexpr static quantity  = HybridQuantity::Scalar::Vx;

    auto const pid = 0;
    auto& p0       = self.patches[pid]; // origin patch
    UsableVecField_t V{"V", p0.layout, HybridQuantity::Vector::V};
    auto const ss0 = make_qty_nd_span_set_from<quantity>(V[0], p0.layout);
    EXPECT_EQ(ss0.box().size(), p0.layout.AMRGhostBoxFor(quantity).size());
    EXPECT_EQ(ss0.box().size(),
              static_cast<std::size_t>(std::pow(cells + (nbrGhosts * 2) + 1, self.dim)));
}



template<typename TileSetTest_t>
auto from_layout_for_quantity(TileSetTest_t& self)
{
    using GridLayout_t     = TileSetTest_t::GridLayout_t;
    using UsableVecField_t = UsableVecField<GridLayout_t, AllocatorMode::CPU, LayoutMode::AoSTS>;
    auto constexpr static nbrGhosts = GridLayout_t::nbrGhosts();
    auto constexpr static quanitity = HybridQuantity::Scalar::Vx;

    auto const pid = 0;
    auto& p0       = self.patches[pid]; // origin patch
    auto const ss0 = make_nd_span_set_for_qty(p0.layout, quanitity);
    EXPECT_EQ(ss0.box().size(), p0.layout.AMRGhostBoxFor(quanitity).size());
    EXPECT_EQ(ss0.box().size(),
              static_cast<std::size_t>(std::pow(cells + (nbrGhosts * 2) + 1, self.dim)));
}


// clang-format off
using Permutations_t = testing::Types<
    TestParam<1, 1, LayoutMode::AoSTS, AllocatorMode::CPU>
   ,TestParam<2, 2, LayoutMode::AoSTS, AllocatorMode::CPU>
   ,TestParam<3, 3, LayoutMode::AoSTS, AllocatorMode::CPU>
PHARE_WITH_GPU(
   ,TestParam<3, 3, LayoutMode::AoSTS, AllocatorMode::GPU_UNIFIED>
)
>;
// clang-format on

TYPED_TEST_SUITE(TileSetOverlapsTest, Permutations_t, );

TYPED_TEST(TileSetOverlapsTest, patch_ghost_copy)
{
    patch_ghost_copy(*this);
}

TYPED_TEST(TileSetOverlapsTest, patch_tile_field)
{
    patch_tile_field(*this);
}


TYPED_TEST(TileSetOverlapsTest, from_layout_for_quantity)
{
    from_layout_for_quantity(*this);
}



} // namespace PHARE::core


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
