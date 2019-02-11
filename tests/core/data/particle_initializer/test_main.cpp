

#include "data/grid/gridlayout.h"
#include "data/grid/gridlayoutimplyee.h"
#include "data/ions/particle_initializers/particle_initializer_factory.h"
#include "data/particles/particle_array.h"
#include "data_provider.h"


#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <type_traits>

using namespace PHARE::core;
using namespace PHARE::initializer;


using GridLayoutT    = GridLayout<GridLayoutImplYee<1, 1>>;
using ParticleArrayT = ParticleArray<1>;



double density(double x)
{
    return 1.;
}


std::array<double, 3> bulkVelocity(double x)
{
    return {{1.0, 0., 0.}};
}



std::array<double, 3> thermalVelocity(double x)
{
    return {{0.2, 0.2, 0.2}};
}




TEST(AParticleIinitializerFactory, takesAPHAREDictToCreateAParticleInitializer)
{
    PHARE::initializer::PHAREDict<1> dict;
    dict["name"]            = std::string{"FluidParticleInitializer"};
    dict["density"]         = static_cast<PHARE::initializer::ScalarFunction<1>>(density);
    dict["bulkVelocity"]    = static_cast<PHARE::initializer::VectorFunction<1>>(bulkVelocity);
    dict["thermalVelocity"] = static_cast<PHARE::initializer::VectorFunction<1>>(thermalVelocity);
    dict["charge"]          = 1.;
    dict["nbrPartPerCell"]  = std::size_t{100};
    dict["basis"]           = std::string{"Cartesian"};

    auto initializer = ParticleInitializerFactory<ParticleArrayT, GridLayoutT>::create(dict);
}



int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
