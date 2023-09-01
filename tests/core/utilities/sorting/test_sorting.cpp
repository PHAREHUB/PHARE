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
    std::size_t constexpr static nbrParticles           = 100 * product(lastCellIndex);

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
    CountingSort<ParticleArray<TypeParam::value>, TypeParam::value> sorter;

    // when allocation is ok find a way to allocate room enough for
    // multiple sorting
    sorter.set_size(this->nbrParticles);

    // when ready to sort
    // auto cellFattener = [](auto const& particle) { return particle.iCell[0]; };
    // sorter.sort(particles, cellFattener);
}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
