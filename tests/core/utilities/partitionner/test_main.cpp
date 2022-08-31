
#include <algorithm>
#include <random>
#include <string>
#include <vector>

#include "core/data/particles/particle_array.hpp"
#include "core/utilities/box/box.hpp"
#include "core/utilities/partitionner/partitionner.hpp"
#include "core/utilities/point/point.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

using PHARE::core::Box;
using PHARE::core::Particle;
using PHARE::core::ParticleArray;
using PHARE::core::Point;




class APartitionner : public ::testing::Test
{
public:
    auto splitLeaving()
    {
        return std::partition(std::begin(particles), std::end(particles),
                              [this](Particle<2> const& part) {
                                  return PHARE::core::isIn(PHARE::core::cellAsPoint(part),
                                                           patchBox);
                              }

        );
    }


    void setupBoxes()
    {
        auto right  = Box<int, 2>{Point<int, 2>{21, 0}, Point<int, 2>{22, 10}};
        auto corner = Box<int, 2>{Point<int, 2>{21, 11}, Point<int, 2>{22, 12}};
        auto up     = Box<int, 2>{Point<int, 2>{0, 11}, Point<int, 2>{20, 11}};

        boundaryBoxes.push_back(std::move(right));
        boundaryBoxes.push_back(std::move(up));
        boundaryBoxes.push_back(std::move(corner));
    }

    void setupParticles()
    {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> disInX(0, 20);
        std::uniform_int_distribution<> disInY(0, 10);
        std::uniform_int_distribution<> disInBox0Y(0, 10);
        std::uniform_int_distribution<> disInBox1X(0, 20);

        for (int i = 0; i < 850; ++i)
        {
            particles[i].iCell[0] = disInX(gen);
            particles[i].iCell[1] = disInY(gen);
        }

        for (int i = 850; i < 900; ++i)
        {
            particles[i].iCell[0] = 21;
            particles[i].iCell[1] = disInBox0Y(gen);
        }

        for (int i = 900; i < 950; ++i)
        {
            particles[i].iCell[0] = disInBox1X(gen);
            particles[i].iCell[1] = 11;
        }

        for (int i = 950; i < 1000; ++i)
        {
            particles[i].iCell[0] = 21;
            particles[i].iCell[1] = 11;
        }

        for (int i = 1000; i < 1025; ++i)
        {
            particles[i].iCell[0] = -1;
            particles[i].iCell[1] = disInBox0Y(gen);
        }
        for (int i = 1025; i < 1050; ++i)
        {
            particles[i].iCell[0] = disInBox1X(gen);
            particles[i].iCell[1] = -1;
        }
    }


    APartitionner()
        : particles(1050)
        , patchBox{Point<int, 2>{0, 0}, Point<int, 2>{20, 10}}
    {
        setupParticles();

        setupBoxes();

        // firstLeaving points to the first particle
        // that is NOT in the patch domain
        // particles from this pointer are either outside in
        // physical boundary boxes or outside the patch but still
        // inside the physical domain.
        firstLeaving = splitLeaving();
    }



protected:
    using ParticleArray_t = ParticleArray<Particle<2>>;
    ParticleArray_t particles;
    Box<int, 2> patchBox;
    std::vector<Box<int, 2>> boundaryBoxes;
    ParticleArray_t::iterator firstLeaving;
};



TEST_F(APartitionner, returnsNbrBoxPlusOneIterators)
{
    auto partitions = partitionner(firstLeaving, std::end(particles), boundaryBoxes);
    EXPECT_EQ(boundaryBoxes.size() + 1, partitions.size());
}




TEST_F(APartitionner, sortsParticlesAccordingToBoxTheyAreIn)
{
    auto partitions = partitionner(firstLeaving, std::end(particles), boundaryBoxes);
    EXPECT_EQ(50, std::distance(partitions[0], partitions[1]));
    EXPECT_TRUE(std::all_of(partitions[0], partitions[1], [this](auto const& part) {
        return isIn(cellAsPoint(part), boundaryBoxes[0]);
    }));
    EXPECT_EQ(50, std::distance(partitions[1], partitions[2]));
    EXPECT_TRUE(std::all_of(partitions[1], partitions[2], [this](auto const& part) {
        return isIn(cellAsPoint(part), boundaryBoxes[1]);
    }));
    EXPECT_EQ(50, std::distance(partitions[2], partitions[3]));
    EXPECT_TRUE(std::all_of(partitions[2], partitions[3], [this](auto const& part) {
        return isIn(cellAsPoint(part), boundaryBoxes[2]);
    }));
}

TEST_F(APartitionner, putsAllLeavingParticlesAtTheEnd)
{
    auto partitions = partitionner(firstLeaving, std::end(particles), boundaryBoxes);
    std::vector<Box<int, 2>> patchAndBoundaries;
    patchAndBoundaries.push_back(patchBox);
    for (auto const& boundaryBox : boundaryBoxes)
        patchAndBoundaries.push_back(boundaryBox);
    EXPECT_TRUE(
        std::none_of(partitions[3], std::end(particles), [patchAndBoundaries](auto const& part) {
            return isIn(cellAsPoint(part), patchAndBoundaries);
        }));
}




int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
