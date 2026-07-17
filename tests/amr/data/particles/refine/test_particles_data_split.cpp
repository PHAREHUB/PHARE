
#include "core/data/particles/particle_array_def.hpp"
#include "core/utilities/box/box.hpp"
#include "core/utilities/types.hpp"

#include "amr/data/particles/refine/particles_data_split.hpp"
#include "amr/data/particles/particles_data.hpp"

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
std::size_t constexpr static ppc    = 10;
std::size_t constexpr static cells  = 14;

auto constexpr nb_split_parts(auto const dim)
{
    if (dim == 1)
        return 2;
    if (dim == 2)
        return 4;
    throw std::runtime_error("no impl for dim");
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
        = TestGridLayout<typename core::PHARE_Types<SimOpts{dim, interp}>::GridLayout_t>;
    using ParticleArray_t
        = ParticleArray<ParticleArrayOptions{dim, layout_mode, StorageMode::VECTOR, alloc_mode}>;
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
    SAMRAI::hier::IntVector static inline const ghostVec{dimension,
                                                         GridLayout_t::nbrParticleGhosts()};

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

    using GridLayout_t = TestParam::GridLayout_t;

    using Splitter
        = PHARE::amr::Splitter<PHARE::core::DimConst<TestParam::dim>,
                               PHARE::core::InterpConst<interp>,
                               PHARE::core::RefinedParticlesConst<nb_split_parts(TestParam::dim)>>;

    auto static constexpr nghosts = GridLayout_t::nbrParticleGhosts();

    using Super::L1;
    using Super::patches;

    ParticlesDataTest()
    {
        assert(!any_overlaps_in(patches, [](auto const& patch) { return patch.layout.AMRBox(); }));
        for (auto& patch : patches)
            add_particles_in(patch.data->domainParticles, patch.layout.AMRBox(), ppc);
    }


    auto refineDomain(auto& src, auto& dst)
    {
        using Refiner
            = ParticlesRefining<ParticleArray_t, ParticlesDataSplitType::interior, Splitter>;
        std::array const boxes{dst.layout.AMRBox()};
        Refiner{*src.data, *dst.data}.forBoxes(boxes);
    }

    auto refineLevelGhost(auto& src, auto& dst)
    {
        using Refiner
            = ParticlesRefining<ParticleArray_t, ParticlesDataSplitType::coarseBoundary, Splitter>;
        auto const boxes = grow(dst.layout.AMRBox(), nghosts).remove(dst.layout.AMRBox());
        Refiner{*src.data, *dst.data}.forBoxes(boxes);
    }
};


// clang-format off
using ParticlesDatas = testing::Types< //

    TestParam<1, LayoutMode::AoSMapped, AllocatorMode::CPU>
//    ,TestParam<1, LayoutMode::AoS, AllocatorMode::CPU>

PHARE_WITH_MKN_GPU(
   ,TestParam<1, LayoutMode::AoSTS, AllocatorMode::CPU>
   ,TestParam<1, LayoutMode::AoSPCTS, AllocatorMode::CPU>
)

//    ,TestParam<2, LayoutMode::AoSMapped, AllocatorMode::CPU>
//    ,TestParam<2, LayoutMode::AoS, AllocatorMode::CPU>

// PHARE_WITH_THRUST(
//    // ,TestParam<2, LayoutMode::AoSTS, AllocatorMode::CPU>
// )

// PHARE_WITH_GPU(
//    ,TestParam<2, LayoutMode::AoS, AllocatorMode::GPU_UNIFIED>

// PHARE_WITH_THRUST(
//    // ,TestParam<2, LayoutMode::AoSTS, AllocatorMode::GPU_UNIFIED>
// )

// )


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

    std::size_t const expected = this->L1.layout.AMRBox().size() * ppc;
    EXPECT_EQ(expected, dst.data->domainParticles.size());

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

    std::size_t const expected = (ghost_box.size() - L1domainBox.size()) * ppc;
    EXPECT_EQ(expected, dst.data->levelGhostParticles.size());

    dst.data->levelGhostParticles.check();

    per_particle(dst.data->levelGhostParticles, [&](auto const& p) {
        EXPECT_NE(p.v()[0], 0);
        EXPECT_NE(p.v()[1], 0);
        EXPECT_NE(p.v()[2], 0);
    });
}


} // namespace PHARE::amr


int main(int argc, char** argv)
{
    PHARE::test::amr::SamraiLifeCycle samsam{argc, argv};
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
