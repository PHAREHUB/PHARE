
#include <cstddef>
#include <iterator>
#include <string>
#include <vector>
#include <random>

#include "core/utilities/cellmap.h"
#include "core/utilities/box/box.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

using namespace PHARE::core;

struct Obj
{
};

template<std::size_t dim>
struct Particle
{
    std::array<int, dim> iCell;
};

TEST(CellMap, isDefaultConstructible)
{
    auto constexpr dim         = 2u;
    auto constexpr bucket_size = 100u;
    CellMap<dim, Obj, bucket_size> cm;
}




TEST(CellMap, sizeIsNbrOfElements)
{
    auto constexpr dim         = 2u;
    auto constexpr bucket_size = 99u;
    CellMap<dim, Obj, bucket_size> cm;
    cm.addToCell(std::array<int, 2>{14, 27}, Obj{});
    EXPECT_EQ(cm.size({14, 27}), 1);
    cm.addToCell(Point{14, 26}, Obj{});
    EXPECT_EQ(cm.size({14, 26}), 1);
    EXPECT_EQ(cm.size({14, 25}), 0);
    EXPECT_EQ(cm.size({26, 14}), 0);
}




template<std::size_t dim>
auto make_particles_in(PHARE::core::Box<int, dim> box, std::size_t nppc)
{
    std::vector<Particle<dim>> particles;
    particles.reserve(box.size() * nppc);

    for (auto const& cell : box)
    {
        for (auto ip = 0u; ip < nppc; ++ip)
        {
            Particle<dim> p;
            for (auto idim = 0u; idim < dim; ++idim)
            {
                p.iCell[idim] = cell[idim];
            }
            particles.push_back(p);
        }
    }

    return particles;
}



TEST(CellMap, canStoreItemCollection)
{
    auto constexpr dim         = 2u;
    auto constexpr bucket_size = 100u;
    CellMap<dim, Particle<dim>, bucket_size> cm;
    Box<int, dim> patchbox{{10, 11}, {20, 21}};
    auto nppc      = 100u;
    auto particles = make_particles_in(patchbox, nppc);
    EXPECT_EQ(patchbox.size() * nppc, particles.size());
    EXPECT_EQ(cm.size(), 0);
    EXPECT_EQ(cm.capacity(), 0);
    cm.add(particles);
    EXPECT_EQ(cm.size(), particles.size());
    EXPECT_EQ(cm.capacity(), cm.size());
}


TEST(CellMap, canStoreItemCollectionFromIterators)
{
    auto constexpr dim         = 2u;
    auto constexpr bucket_size = 100u;
    CellMap<dim, Particle<dim>, bucket_size> cm;
    Box<int, dim> patchbox{{10, 11}, {20, 21}};
    auto nppc      = 100u;
    auto particles = make_particles_in(patchbox, nppc);
    EXPECT_EQ(patchbox.size() * nppc, particles.size());
    EXPECT_EQ(cm.size(), 0);
    EXPECT_EQ(cm.capacity(), 0);
    cm.add(std::begin(particles), std::end(particles));
    EXPECT_EQ(cm.size(), particles.size());
    EXPECT_EQ(cm.capacity(), cm.size());
}



TEST(CellMap, givesAccessToAllParticlesInACell)
{
    auto constexpr dim         = 2u;
    auto constexpr bucket_size = 100u;
    Box<int, 2> patchbox{{10, 20}, {25, 42}};
    CellMap<dim, Particle<dim>, bucket_size> cm;
    auto nppc      = 100u;
    auto particles = make_particles_in(patchbox, nppc);
    cm.add(particles);

    for (auto const& cell : patchbox)
    {
        auto& blist = cm.list_at(cell);
        EXPECT_EQ(blist.size(), nppc);
        for (auto const& particle : blist)
        {
            EXPECT_EQ(Point{particle->iCell}, cell);
        }
    }
}



TEST(CellMap, emptyLeavesCapacityButZeroSize)
{
    auto constexpr dim         = 2u;
    auto constexpr bucket_size = 100u;
    Box<int, 2> patchbox{{10, 20}, {25, 42}};
    CellMap<dim, Particle<dim>, bucket_size> cm;
    auto nppc      = 100u;
    auto particles = make_particles_in(patchbox, nppc);
    cm.add(particles);

    cm.empty();
    EXPECT_EQ(cm.size(), 0);
    EXPECT_TRUE(cm.is_empty());
    EXPECT_EQ(cm.capacity(), patchbox.size() * nppc);
}




TEST(CellMap, selectParticleInSubsetBox)
{
    auto constexpr dim         = 2u;
    auto constexpr bucket_size = 100u;
    Box<int, 2> patchbox{{10, 20}, {25, 42}};
    CellMap<dim, Particle<dim>, bucket_size> cm;
    auto nppc      = 100u;
    auto particles = make_particles_in(patchbox, nppc);
    cm.add(particles);

    Box<int, 2> selection_box{{12, 22}, {17, 23}};
    auto selected = cm.select(particles, selection_box);
    for (auto const& particle : selected)
    {
        EXPECT_TRUE(isIn(Point{particle.iCell}, selection_box));
    }
    std::size_t isInCounter = 0;
    for (auto const& particle : particles)
    {
        if (isIn(Point{particle.iCell}, selection_box))
        {
            isInCounter++;
        }
    }
    EXPECT_EQ(isInCounter, selected.size());
}


TEST(CellMap, trimMemory)
{
    auto constexpr dim         = 2u;
    auto constexpr bucket_size = 100u;
    Box<int, 2> patchbox{{10, 20}, {25, 42}};
    CellMap<dim, Particle<dim>, bucket_size> cm;
    auto nppc      = 200u;
    auto particles = make_particles_in(patchbox, nppc);
    cm.add(particles);

    EXPECT_EQ(cm.size(), patchbox.size() * nppc);
    EXPECT_EQ(patchbox.size() * nppc, cm.capacity());

    cm.empty();
    EXPECT_EQ(cm.size(), 0);
    EXPECT_EQ(cm.capacity(), patchbox.size() * nppc);
    nppc      = 100;
    particles = make_particles_in(patchbox, nppc);
    cm.add(particles);
    EXPECT_EQ(cm.size(), patchbox.size() * nppc);
    EXPECT_EQ(patchbox.size() * nppc * 2, cm.capacity());
    EXPECT_EQ(0.5, cm.used_mem_ratio());
    cm.trim(1);
    EXPECT_EQ(patchbox.size() * nppc * 2, cm.capacity());
    cm.trim(0);
    EXPECT_EQ(patchbox.size() * nppc * 1, cm.capacity());
}


TEST(CellMap, emptyTrimsWithMaxEmpty1IfCapacityExceeds3TimesSize)
{
    auto constexpr dim         = 2u;
    auto constexpr bucket_size = 100u;
    Box<int, 2> patchbox{{10, 20}, {25, 42}};
    using cellmap_t = CellMap<dim, Particle<dim>, bucket_size>;
    cellmap_t cm;
    auto nppc      = 600u;
    auto particles = make_particles_in(patchbox, nppc);
    cm.add(particles);
    cm.empty();
    nppc      = 100u;
    particles = make_particles_in(patchbox, nppc);
    cm.add(particles);
    EXPECT_EQ(cm.capacity(), 600 * patchbox.size());
    cm.empty();
    EXPECT_EQ(cm.capacity(), particles.size() * 2);
}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}

