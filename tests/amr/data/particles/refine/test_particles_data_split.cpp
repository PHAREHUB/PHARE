
#include "core/utilities/box/box.hpp"
#include "core/utilities/types.hpp"

#include "core/data/particles/particle_array_def.hpp"

#include "amr/data/particles/particles_data.hpp"
#include "amr/data/particles/refine/particles_data_split.hpp"

#include "simulator/simulator_def.hpp"

#include "tests/core/data/particles/test_particles.hpp"
#include "tests/core/data/gridlayout/test_gridlayout.hpp"

#include "tests/amr/amr.hpp"

#include "gtest/gtest.h"
#include <SAMRAI/pdat/CellGeometry.h>
#include <stdexcept>

using namespace PHARE;
using namespace PHARE::core;
using namespace PHARE::amr;


std::size_t constexpr static interp = 1;
std::size_t constexpr static ppc    = 100;
std::size_t constexpr static cells  = 30;

auto constexpr nb_split_parts(auto const dim)
{
    if (dim == 1)
        return 2;
    if (dim == 2)
        return 4;
    if (dim == 3)
        return 6;
    throw std::runtime_error("no impl for dim");
}


// walks whatever the actual storage looks like (flat vector, tiles of flat vectors,
// tiles of per-cell vectors) and accumulates size/capacity from the same leaf
// containers, so the two are always comparable regardless of ParticleArray_t's layout.
template<typename ParticleArray_t>
auto total_size_and_capacity(ParticleArray_t const& particles)
{
    std::size_t size = 0, capacity = 0;

    if constexpr (ParticleArray_t::layout_mode == LayoutMode::AoSPCTS)
    {
        for (auto const& tile : particles())
        {
            auto const& pc = tile();
            for (auto const& bix : pc.ghost_box())
            {
                auto const& cell = pc(pc.local_cell(bix));
                size += cell.size();
                capacity += cell.capacity();
            }
        }
    }
    else if constexpr (is_tiled(ParticleArray_t::layout_mode))
    {
        for (auto const& tile : particles())
        {
            size += tile().size();
            capacity += tile().capacity();
        }
    }
    else
    {
        size     = particles.size();
        capacity = particles.capacity();
    }

    return std::make_pair(size, capacity);
}



template<std::size_t _dim, auto lm, auto am>
struct TestParam
{
    static_assert(all_are<LayoutMode>(lm));
    static_assert(all_are<AllocatorMode>(am));

    auto constexpr static dim         = _dim;
    auto constexpr static layout_mode = lm;
    auto constexpr static alloc_mode  = am;

    using Box_t = PHARE::core::Box<int, dim>;
    using GridLayout_t
        = TestGridLayout<typename core::PHARE_Types<SimOpts{dim, interp}>::Hybrid::GridLayout_t>;
    using ParticleArray_t
        = ParticleArray<ParticleArrayOptions{dim, layout_mode, StorageMode::VECTOR, alloc_mode}>;

    // stands in for a "HybridTypes" (as PHARE::amr::ParticlesRefining expects) that uses
    // this test's specific ParticleArray_t layout/alloc combo instead of the default one -
    // everything else (hybrid_options, particle_ghost_width, GridLayout_t) is unaffected by
    // particle layout/alloc mode, so it's reused as-is.
    struct HybridTypes_t : core::PHARE_Types<SimOpts{dim, interp}>::Hybrid
    {
        using ParticleArray_t = TestParam::ParticleArray_t;
    };
};



template<typename TestParam>
struct Patch
{
    auto constexpr static dim = TestParam::dim;

    using ParticleArray_t = TestParam::ParticleArray_t;
    using ParticlesData_t = ParticlesData<ParticleArray_t>;
    using GridLayout_t    = TestParam::GridLayout_t;
    using Box_t           = TestParam::Box_t;

    SAMRAI::tbox::Dimension static inline const dimension{dim};
    SAMRAI::hier::BlockId static inline const blockId{0};
    SAMRAI::hier::IntVector static inline const ghostVec{
        dimension, GridLayout_t::options.particle_ghost_width};

    Patch(Box_t const& box)
        : layout{box}
    {
    }
    Patch(GridLayout_t const& _layout)
        : layout{_layout}
    {
    }


    GridLayout_t const layout;
    SAMRAI::hier::Box const domain{samrai_box_from(layout.AMRBox())};
    std::shared_ptr<SAMRAI::hier::BoxGeometry> const geom{
        std::make_shared<SAMRAI::pdat::CellGeometry>(domain, ghostVec)};

    std::shared_ptr<ParticlesData<ParticleArray_t>> data{
        std::make_shared<ParticlesData<ParticleArray_t>>(domain, ghostVec, "name")};

    SAMRAI::hier::Box const mask{data->getGhostBox()};
};

template<std::size_t dim, typename TestParam>
struct AParticlesDataTest;

template<typename TestParam>
struct AParticlesDataTest<1, TestParam>
{
    using Box_t = TestParam::Box_t;

    AParticlesDataTest()
    {
        patches.reserve(3);
        auto const off = cells - 1;
        for (std::uint8_t i = 0; i < 3; ++i)
        {
            auto const cellx = i * cells;
            patches.emplace_back(Box_t{Point{cellx}, Point{cellx + off}});
        }
    }

    std::vector<Patch<TestParam>> patches;
    Patch<TestParam> L1{Box_t{Point{4}, Point{9}}};
};


template<typename TestParam>
struct AParticlesDataTest<2, TestParam>
{
    using Box_t = TestParam::Box_t;

    AParticlesDataTest()
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
    }

    std::vector<Patch<TestParam>> patches;
    Patch<TestParam> L1{Box_t{Point{3, 3}, Point{12, 12}}};
};


template<typename TestParam>
struct AParticlesDataTest<3, TestParam>
{
    using Box_t = TestParam::Box_t;

    AParticlesDataTest()
    {
        patches.reserve(3 * 3);
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
    }

    std::vector<Patch<TestParam>> patches;
    Patch<TestParam> L1{Box_t{Point{3, 3, 3}, Point{12, 12, 12}}};
};


template<typename TestParam>
struct ParticlesDataTest : public ::testing::Test,
                           public AParticlesDataTest<TestParam::dim, TestParam>
{
    using Super           = AParticlesDataTest<TestParam::dim, TestParam>;
    using ParticleArray_t = TestParam::ParticleArray_t;
    using HybridTypes_t   = TestParam::HybridTypes_t;

    using GridLayout_t = TestParam::GridLayout_t;

    using Splitter
        = PHARE::amr::Splitter<PHARE::core::DimConst<TestParam::dim>,
                               PHARE::core::InterpConst<interp>,
                               PHARE::core::RefinedParticlesConst<nb_split_parts(TestParam::dim)>>;

    auto static constexpr nghosts = GridLayout_t::options.particle_ghost_width;

    using Super::L1;
    using Super::patches;

    ParticlesDataTest()
    {
        assert(!any_overlaps_in(patches, [](auto const& patch) { return patch.layout.AMRBox(); }));
        for (auto& patch : patches)
            add_particles(patch.data->domainParticles, patch.layout.AMRBox(), ppc);
    }


    auto refineDomain(auto& src, auto& dst)
    {
        using Refiner
            = ParticlesRefining<HybridTypes_t, ParticlesDataSplitType::interior, Splitter>;
        SAMRAI::hier::BoxContainer boxes;
        boxes.pushBack(samrai_box_from(dst.layout.AMRBox()));
        Refiner{*src.data, *dst.data}.forBoxes(boxes);
    }

    auto refineLevelGhost(auto& src, auto& dst)
    {
        using Refiner
            = ParticlesRefining<HybridTypes_t, ParticlesDataSplitType::coarseBoundary, Splitter>;
        auto const sGhostBox  = samrai_box_from(grow(dst.layout.AMRBox(), nghosts));
        auto const sDomainBox = samrai_box_from(dst.layout.AMRBox());
        SAMRAI::hier::BoxContainer boxes;
        boxes.removeIntersections(sGhostBox, sDomainBox);
        Refiner{*src.data, *dst.data}.forBoxes(boxes);
    }
};


// clang-format off
using ParticlesDatas = testing::Types< //

    TestParam<1, LayoutMode::AoSMapped, AllocatorMode::CPU>
   ,TestParam<1, LayoutMode::AoSPCTS, AllocatorMode::CPU>
   ,TestParam<3, LayoutMode::AoSPCTS, AllocatorMode::CPU>

>;
// clang-format on

TYPED_TEST_SUITE(ParticlesDataTest, ParticlesDatas, );

namespace PHARE::amr
{


TYPED_TEST(ParticlesDataTest, splitWorksForDomain)
{
    auto& dst = this->L1;

    for (auto const& src : this->patches)
        this->refineDomain(src, dst);

    // nb_split_parts(dim) refined particles come from one coarse particle, spread across
    // the 2**dim fine cells it refines into - only equal to 2**dim (i.e. density/count
    // preserving) for the 1D/2D patterns used here; the 3D pattern (6) instead conserves
    // total weight via SplitPattern::weight_, not raw per-cell particle count, so the
    // expected count must be scaled by the same ratio the split itself uses.
    std::size_t const expected = this->L1.layout.AMRBox().size() * ppc
                                 * nb_split_parts(TestFixture::ParticleArray_t::dimension)
                                 / (1u << TestFixture::ParticleArray_t::dimension);
    EXPECT_EQ(expected, dst.data->domainParticles.size());

    // reserves are expected to closely track the final size, not grossly over-allocate.
    // AoSPCTS per-cell storage grows each cell via its own vector's strategy while
    // streaming, then shrinks every touched cell to its exact final size once at the end
    // (see aos_refining.hpp), so this holds regardless of how unevenly a split pattern
    // distributes particles across a coarse cell's children.
    auto const [total_size, total_capacity] = total_size_and_capacity(dst.data->domainParticles);
    EXPECT_EQ(expected, total_size);
    EXPECT_LE(total_capacity, static_cast<std::size_t>(total_size * 2.0)); // todo

    dst.data->domainParticles.check();

    per_particle(dst.data->domainParticles, [&](auto const& p) {
        EXPECT_TRUE(isIn(p, dst.layout.AMRBox()));
        EXPECT_NE(p.v()[0], 0);
        EXPECT_NE(p.v()[1], 0);
        EXPECT_NE(p.v()[2], 0);
    });

    // for (auto const& bix : dst.layout.AMRBox())
    // {
    //     EXPECT_TRUE(sum_from(dst.data->domainParticles,
    //                          [&](auto const& p) { return array_equals(p.iCell(), *bix) ? 1 : 0;
    //                          })
    //                 == ppc)
    //         << "failed for " << bix;
    // }
}


TYPED_TEST(ParticlesDataTest, splitWorksForLevelGhost)
{
    auto& dst = this->L1;

    for (auto const& src : this->patches)
        this->refineLevelGhost(src, dst);

    auto const& L1domainBox = this->L1.layout.AMRBox();
    auto const ghost_box    = grow(L1domainBox, TestFixture::nghosts);

    // see splitWorksForDomain for why this scales by nb_split_parts(dim) / 2**dim
    std::size_t const expected = (ghost_box.size() - L1domainBox.size()) * ppc
                                 * nb_split_parts(TestFixture::ParticleArray_t::dimension)
                                 / (1u << TestFixture::ParticleArray_t::dimension);
    EXPECT_EQ(expected, dst.data->levelGhostParticles.size());

    // reserves are expected to closely track the final size, not grossly over-allocate
    // (see splitWorksForDomain for why this holds for AoSPCTS too)
    auto const [total_size, total_capacity]
        = total_size_and_capacity(dst.data->levelGhostParticles);
    EXPECT_EQ(expected, total_size);
    EXPECT_LE(total_capacity, static_cast<std::size_t>(total_size * 2.0));

    dst.data->levelGhostParticles.check();

    per_particle(dst.data->levelGhostParticles, [&](auto const& p) {
        EXPECT_NE(p.v()[0], 0);
        EXPECT_NE(p.v()[1], 0);
        EXPECT_NE(p.v()[2], 0);
    });

    if constexpr (TestFixture::ParticleArray_t::layout_mode == LayoutMode::AoSPCTS)
    {
        for (auto const& tile : dst.data->levelGhostParticles())
        {
            auto const& pc = tile();
            for (auto const& bix : pc.ghost_box())
            {
                auto const local = pc.local_cell(bix);
                for (auto const& p : pc(local))
                    EXPECT_TRUE(array_equals(pc.local_cell(p.iCell()), local))
                        << "levelGhost particle in wrong cell bucket: iCell=" << Point{p.iCell()}
                        << " bucket=" << Point{local} << " tile="
                        << static_cast<
                               core::Box<int, TestFixture::ParticleArray_t::dimension> const&>(
                               tile);
            }
        }
    }
}


} // namespace PHARE::amr


int main(int argc, char** argv)
{
    PHARE::test::amr::SamraiLifeCycle samsam{argc, argv};
    ::testing::InitGoogleTest(&argc, argv);
    auto const r = RUN_ALL_TESTS();
    return r;
}
