#include "core/utilities/sorting.hpp"

#include "core/data/particles/particle_array.hpp"
#include "core/utilities/box/box.hpp"
#include "core/utilities/sorting.hpp"
#include "core/utilities/types.hpp"

#include "gtest/gtest.h"


#include <type_traits>
#include <array>
#include <random>


using namespace PHARE::core;


template<typename Dim>
class SortingTest : public ::testing::Test
{
public:
    std::size_t constexpr static dim                    = Dim::value;
    std::array<int, dim> constexpr static lastCellIndex = ConstArray<int, dim>(100);
    std::size_t constexpr static nbrParticles           = 2 * product(lastCellIndex);

    SortingTest()
        : box{ConstArray<int, dim>(0), lastCellIndex}
        , particles(box, nbrParticles)
    {
        std::random_device rd;
        std::mt19937_64 gen(rd());
        std::uniform_int_distribution<int> dist(0, lastCellIndex[0]);

        for (std::size_t iPart = 0; iPart < nbrParticles; ++iPart)
        {
            for (std::size_t iDim = 0; iDim < dim; ++iDim)
            {
                particles[iPart].iCell[iDim] = dist(gen);
            }
        }
    }

    Box<int, dim> box;
    ParticleArray<dim> particles;
};

using ConcreteSortingTests = ::testing::Types<std::integral_constant<std::size_t, 1u>,
                                              std::integral_constant<std::size_t, 2u>,
                                              std::integral_constant<std::size_t, 3u>>;

TYPED_TEST_SUITE(SortingTest, ConcreteSortingTests);


TYPED_TEST(SortingTest, CoutingSort)
{
    constexpr auto dim = TypeParam::value;
    CountingSort<ParticleArray<dim>, dim> sorter;

    // when allocation is ok find a way to allocate room enough for
    // multiple sorting
    sorter.setup(this->nbrParticles, this->box);

    // when ready to sort
    auto cellFattener = [&](auto const& particle) {
        if constexpr (dim == 1)
        {
            return particle.iCell[0];
        }
        else if constexpr (dim == 2)
        {
            auto const ny = this->lastCellIndex[1] + 1;
            return particle.iCell[1] + particle.iCell[0] * ny;
        }
        else if constexpr (dim == 3)
        {
            auto const ny = this->lastCellIndex[1] + 1;
            auto const nz = this->lastCellIndex[2] + 1;
            return particle.iCell[2] + particle.iCell[1] * nz + particle.iCell[0] * ny * nz;
        }
    };
    sorter.sort(this->particles, cellFattener);


    // check that the particles are sorted
    for (std::size_t iPart = 0; iPart < this->nbrParticles; ++iPart)
    {
        if (iPart > 0)
        {
            EXPECT_LE(cellFattener(this->particles[iPart - 1]),
                      cellFattener(this->particles[iPart]));
        }
    }
}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
