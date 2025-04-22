
#include "phare_core.hpp"
#include "phare/phare.hpp"
#include "core/utilities/types.hpp"
#include "amr/data/particles/particles_data.hpp"

#include "tests/core/data/gridlayout/test_gridlayout.hpp"
#include "tests/core/data/particles/test_particles_fixtures.hpp"

#include "gtest/gtest.h"

#include <SAMRAI/pdat/CellGeometry.h>


std::size_t constexpr static interp      = 1;
std::size_t constexpr static cells       = 5;
std::size_t constexpr static ppc         = 10;
std::size_t constexpr static ghost_cells = 1;


template<typename AParticlesData>
struct Patch
{
    auto constexpr static dim          = AParticlesData::dim;
    using GridLayout_t                 = AParticlesData::GridLayout_t;
    using ParticleArray_t              = AParticlesData::ParticleArray_t;
    using ParticlesData_t              = PHARE::amr::ParticlesData<ParticleArray_t>;
    auto constexpr static nbPartGhosts = GridLayout_t::nbrParticleGhosts();

    Patch(auto const& box)
        : layout{box}
    {
    }

    auto static overlap(Patch const& src, Patch& dst)
    {
        return overlap(src, dst, PHARE::core::ConstArray<int, dim>());
    }

    auto static overlap(Patch const& src, Patch& dst, auto const shift)
    {
        bool constexpr static overwriteInterior{true};

        SAMRAI::hier::IntVector shiftVec{dimension};
        for (std::size_t i = 0; i < dim; ++i)
            shiftVec[i] = -shift[i];

        SAMRAI::hier::Box const srcMask{src.data->getGhostBox()};
        SAMRAI::hier::Box const fillBox{dst.data->getGhostBox()};
        SAMRAI::hier::Transformation const transformation{shiftVec};
        return std::shared_ptr<SAMRAI::pdat::CellOverlap>{
            std::dynamic_pointer_cast<SAMRAI::pdat::CellOverlap>(dst.geom->calculateOverlap(
                *src.geom, srcMask, fillBox, overwriteInterior, transformation))};
    }

    SAMRAI::tbox::Dimension static inline const dimension{dim};
    SAMRAI::hier::BlockId static inline const blockId{0};
    SAMRAI::hier::IntVector static inline const ghostVec{dimension, ghost_cells};

    GridLayout_t const layout;
    SAMRAI::hier::Box const domain{PHARE::amr::samrai_box_from(layout.AMRBox())};
    std::shared_ptr<SAMRAI::hier::BoxGeometry> const geom{
        std::make_shared<SAMRAI::pdat::CellGeometry>(domain, ghostVec)};
    std::unique_ptr<ParticlesData_t> data{
        std::make_unique<ParticlesData_t>(domain, ghostVec, "name")};
    SAMRAI::hier::Box const mask{data->getGhostBox()};
};

template<typename T, std::size_t D>
auto static make_shift_for(PHARE::core::Box<T, D> const& box)
{
    int const span = cells * 3;
    int const mid  = cells * 3 / 2;

    if constexpr (D == 1)
    {
        auto shifts = PHARE::core::for_N<1, PHARE::core::for_N_R_mode::make_array>(
            [&](auto i) { return PHARE::core::Point<int, D>{0}; });
        PHARE::core::for_N<D>([&](auto i) {
            int const shift = box.upper[i] < mid ? 1 : -1;
            shifts[i][i]    = span * shift;
        });

        return shifts;
    }
    if constexpr (D == 2)
    {
        auto shifts = PHARE::core::for_N<3, PHARE::core::for_N_R_mode::make_array>(
            [&](auto i) { return PHARE::core::Point<int, D>{0, 0}; });
        PHARE::core::for_N<D>([&](auto i) {
            int const shift = box.upper[i] < mid ? 1 : -1;
            shifts[i][i]    = span * shift;
        });

        shifts[2] = {shifts[0][0], shifts[1][1]};

        return shifts;
    }
    if constexpr (D == 3)
    {
        auto shifts = PHARE::core::for_N<7, PHARE::core::for_N_R_mode::make_array>(
            [&](auto i) { return PHARE::core::Point<int, D>{0, 0, 0}; });
        PHARE::core::for_N<D>([&](auto i) {
            int const shift = box.upper[i] < mid ? 1 : -1;
            shifts[i][i]    = span * shift;
        });

        shifts[3] = {shifts[0][0], shifts[1][1], 0};
        shifts[4] = {0, shifts[1][1], shifts[2][2]};
        shifts[5] = {shifts[0][0], 0, shifts[2][2]};
        shifts[6] = {shifts[0][0], shifts[1][1], shifts[2][2]};

        return shifts;
    }
}

template<std::size_t _dim>
struct TestParam
{
    auto constexpr static dim = _dim;
    using PhareTypes          = PHARE::core::PHARE_Types<dim, interp>;
    using GridLayout_t        = TestGridLayout<typename PhareTypes::GridLayout_t>;
    using Box_t               = PHARE::core::Box<int, dim>;
    using ParticleArray_t     = PHARE::core::ParticleArray<dim>;
};



template<std::size_t dim, typename TestParam>
struct AParticlesDataTest;

template<typename TestParam>
struct AParticlesDataTest<1, TestParam>
{
    using Box_t               = TestParam::Box_t;
    using Patch_t             = Patch<TestParam>;
    auto constexpr static dim = TestParam::dim;

    AParticlesDataTest()
    {
        patches.reserve(3);
        auto const off = cells - 1;
        for (std::uint8_t i = 0; i < 3; ++i)
        {
            auto const cellx = i * cells;
            patches.emplace_back(Box_t{PHARE::core::Point{cellx}, PHARE::core::Point{cellx + off}});
        }
    }

    std::vector<Patch_t> patches;
};


template<typename TestParam>
struct AParticlesDataTest<2, TestParam>
{
    using Box_t               = TestParam::Box_t;
    using Patch_t             = Patch<TestParam>;
    auto constexpr static dim = TestParam::dim;

    AParticlesDataTest()
    {
        patches.reserve(3 * 3);
        auto const off = cells - 1;
        for (std::uint8_t i = 0; i < 3; ++i)
            for (std::uint8_t j = 0; j < 3; ++j)
            {
                auto const cellx = i * cells;
                auto const celly = j * cells;
                patches.emplace_back(Box_t{PHARE::core::Point{cellx, celly},
                                           PHARE::core::Point{cellx + off, celly + off}});
            }
    }

    std::vector<Patch_t> patches;
};


template<typename TestParam>
struct AParticlesDataTest<3, TestParam>
{
    using Box_t               = TestParam::Box_t;
    auto constexpr static dim = TestParam::dim;
    using Patch_t             = Patch<TestParam>;

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
                    patches.emplace_back(
                        Box_t{PHARE::core::Point{cellx, celly, cellz},
                              PHARE::core::Point{cellx + off, celly + off, cellz + off}});
                }
    }

    std::vector<Patch_t> patches;
};


template<typename TestParam>
struct ParticlesDataTest : public ::testing::Test,
                           public AParticlesDataTest<TestParam::dim, TestParam>
{
    using Super               = AParticlesDataTest<TestParam::dim, TestParam>;
    using Patch_t             = Super::Patch_t;
    using ParticleArray_t     = TestParam::ParticleArray_t;
    using GridLayout_t        = TestParam::GridLayout_t;
    auto constexpr static dim = TestParam::dim;
    using Super::patches;

    ParticlesDataTest()
    {
        assert(!PHARE::core::any_overlaps_in(
            patches, [](auto const& patch) { return patch.layout.AMRBox(); }));

        for (auto& patch : patches)
            PHARE::core::add_particles_in(patch.data->domainParticles, patch.layout.AMRBox(), ppc);
    }

    auto static overlap(Patch_t const& src, Patch_t& dst) { return Patch_t::overlap(src, dst); }
    auto static overlap(Patch_t const& src, Patch_t& dst, auto const shift)
    {
        return Patch_t::overlap(src, dst, shift);
    }

    auto mid_pid() const { return int(((patches.size() + 1) / 2) - 1); }
};


using ParticlesDatas = testing::Types<TestParam<1>, TestParam<2> /*,TestParam<3>*/>;


TYPED_TEST_SUITE(ParticlesDataTest, ParticlesDatas);

namespace PHARE::core
{
TYPED_TEST(ParticlesDataTest, copyWorks)
{
    auto constexpr static dim = TestFixture::dim;
    std::size_t const pid     = this->mid_pid();
    auto& dst                 = this->patches[pid];

    for (std::size_t i = 0; i < pid; ++i)
        dst.data->copy(*this->patches[i].data);
    for (std::size_t i = pid + 1; i < this->patches.size(); ++i)
        dst.data->copy(*this->patches[i].data);

    std::size_t const expected
        = (std::pow(cells + (ghost_cells * 2), dim) - std::pow(cells, dim)) * ppc;
    EXPECT_EQ(expected, dst.data->patchGhostParticles.size());
}

TYPED_TEST(ParticlesDataTest, copyAtPeriodicBoundaryWorks)
{
    auto constexpr static dim = TestFixture::dim;
    std::size_t const pid     = 0;
    auto& dst                 = this->patches[pid];
    auto const dst_ghostbox   = grow(dst.layout.AMRBox(), 1);

    // non-periodic neighbours
    for (std::size_t i = pid + 1; i < this->patches.size(); ++i)
        if (auto const overlap = dst_ghostbox * this->patches[i].layout.AMRBox())
            dst.data->copy(*this->patches[i].data);

    // // periodic neighbours
    for (auto const& shifter : make_shift_for(dst.layout.AMRBox()))
    {
        auto const shift_box = shift(dst.layout.AMRBox(), shifter);
        auto const ghostbox  = grow(shift_box, 1);

        for (std::size_t i = pid + 1; i < this->patches.size(); ++i)
            if (auto const overlap = ghostbox * this->patches[i].layout.AMRBox())
            {
                auto const celloverlap = TestFixture::overlap(this->patches[i], dst, shifter);
                dst.data->copy(*this->patches[i].data, *celloverlap);
            }
    }

    std::size_t const expected
        = (std::pow(cells + (ghost_cells * 2), dim) - std::pow(cells, dim)) * ppc;
    EXPECT_EQ(expected, dst.data->patchGhostParticles.size());
}

TYPED_TEST(ParticlesDataTest, packWorks)
{
    auto constexpr static dim = TestFixture::dim;
    std::size_t const pid     = this->mid_pid();
    auto& dst                 = this->patches[pid];

    auto for_neighbour = [&](auto i) {
        auto const cellOverlap = TestFixture::overlap(this->patches[i], dst);
        auto const& src        = this->patches[i];
        SAMRAI::tbox::MessageStream particlesWriteStream;
        src.data->packStream(particlesWriteStream, *cellOverlap);
        SAMRAI::tbox::MessageStream particlesReadStream{particlesWriteStream.getCurrentSize(),
                                                        SAMRAI::tbox::MessageStream::Read,
                                                        particlesWriteStream.getBufferStart()};

        dst.data->unpackStream(particlesReadStream, *cellOverlap);
    };

    for (std::size_t i = 0; i < pid; ++i)
        for_neighbour(i);

    for (std::size_t i = pid + 1; i < this->patches.size(); ++i)
        for_neighbour(i);

    std::size_t const expected
        = (std::pow(cells + (ghost_cells * 2), dim) - std::pow(cells, dim)) * ppc;
    EXPECT_EQ(expected, dst.data->patchGhostParticles.size());
}

TYPED_TEST(ParticlesDataTest, packAtPeriodicBoundaryWorks)
{
    auto constexpr static dim = TestFixture::dim;
    std::size_t const pid     = 0;
    auto& dst                 = this->patches[pid];
    auto const dst_ghostbox   = grow(dst.layout.AMRBox(), 1);

    auto for_neighbour = [&](auto i, auto shift) {
        auto const cellOverlap = TestFixture::overlap(this->patches[i], dst, shift);
        auto const& src        = this->patches[i];
        SAMRAI::tbox::MessageStream particlesWriteStream;
        src.data->packStream(particlesWriteStream, *cellOverlap);
        SAMRAI::tbox::MessageStream particlesReadStream{particlesWriteStream.getCurrentSize(),
                                                        SAMRAI::tbox::MessageStream::Read,
                                                        particlesWriteStream.getBufferStart()};

        dst.data->unpackStream(particlesReadStream, *cellOverlap);
    };

    // non-periodic neighbours
    for (std::size_t i = pid + 1; i < this->patches.size(); ++i)
        if (auto const overlap = dst_ghostbox * this->patches[i].layout.AMRBox())
            for_neighbour(i, PHARE::core::ConstArray<int, dim>());

    // periodic neighbours
    for (auto const& shifter : make_shift_for(dst.layout.AMRBox()))
    {
        auto const shift_box = shift(dst.layout.AMRBox(), shifter);
        auto const ghostbox  = grow(shift_box, 1);

        for (std::size_t i = pid + 1; i < this->patches.size(); ++i)
            if (auto const overlap = ghostbox * this->patches[i].layout.AMRBox())
                for_neighbour(i, shifter);
    }

    std::size_t const expected
        = (std::pow(cells + (ghost_cells * 2), dim) - std::pow(cells, dim)) * ppc;
    EXPECT_EQ(expected, dst.data->patchGhostParticles.size());
}

} // namespace PHARE::core


int main(int argc, char** argv)
{
    PHARE::SamraiLifeCycle samsam{argc, argv};
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
