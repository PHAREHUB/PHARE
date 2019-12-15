
#include "core/data/grid/gridlayout.h"
#include "core/data/grid/gridlayoutimplyee.h"
#include "core/data/ions/particle_initializers/particle_initializer_factory.h"
#include "core/data/particles/particle_array.h"
#include "initializer/data_provider.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <type_traits>

using namespace PHARE::core;
using namespace PHARE::initializer;
using GridLayoutT = GridLayout<GridLayoutImplYee<1, 1>>;

double return_one(double /*x*/)
{
    return 1.;
}

template<typename ParticleArray>
void test_particle_array_types()
{
    PHARE::initializer::PHAREDict dict;
    dict["name"]            = std::string{"MaxwellianParticleInitializer"};
    dict["density"]         = static_cast<PHARE::initializer::ScalarFunction<1>>(return_one);
    dict["bulk_velocity_x"] = static_cast<PHARE::initializer::ScalarFunction<1>>(return_one);
    dict["bulk_velocity_y"] = static_cast<PHARE::initializer::ScalarFunction<1>>(return_one);
    dict["bulk_velocity_z"] = static_cast<PHARE::initializer::ScalarFunction<1>>(return_one);

    dict["thermal_velocity_x"] = static_cast<PHARE::initializer::ScalarFunction<1>>(return_one);
    dict["thermal_velocity_y"] = static_cast<PHARE::initializer::ScalarFunction<1>>(return_one);
    dict["thermal_velocity_z"] = static_cast<PHARE::initializer::ScalarFunction<1>>(return_one);

    dict["charge"]         = 1.;
    dict["nbrPartPerCell"] = int{100};
    dict["basis"]          = std::string{"Cartesian"};

    auto initializer = ParticleInitializerFactory<ParticleArray, GridLayoutT>::create(dict);
}

TEST(AParticleIinitializerFactory, takesAPHAREDictToCreateAParticleInitializer_AoS)
{
    test_particle_array_types<ParticleArray<1, false>>();
}
TEST(AParticleIinitializerFactory, takesAPHAREDictToCreateAParticleInitializer_SoA)
{
    test_particle_array_types<ParticleArray<1, true>>();
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
