

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

double density([[maybe_unused]] double x)
{
    return 1.;
}



double vx(double x)
{
    (void)x;
    return 1.;
}


double vy(double x)
{
    (void)x;
    return 1.;
}


double vz(double x)
{
    (void)x;
    return 1.;
}


double vthx(double x)
{
    (void)x;
    return 1.;
}


double vthy(double x)
{
    (void)x;
    return 1.;
}



double vthz(double x)
{
    (void)x;
    return 1.;
}




TEST(AParticleIinitializerFactory, takesAPHAREDictToCreateAParticleInitializer)
{
    PHARE::initializer::PHAREDict dict;
    dict["name"]            = std::string{"MaxwellianParticleInitializer"};
    dict["density"]         = static_cast<PHARE::initializer::ScalarFunction<1>>(density);
    dict["bulk_velocity_x"] = static_cast<PHARE::initializer::ScalarFunction<1>>(vx);
    dict["bulk_velocity_y"] = static_cast<PHARE::initializer::ScalarFunction<1>>(vx);
    dict["bulk_velocity_z"] = static_cast<PHARE::initializer::ScalarFunction<1>>(vx);

    dict["thermal_velocity_x"] = static_cast<PHARE::initializer::ScalarFunction<1>>(vthx);
    dict["thermal_velocity_y"] = static_cast<PHARE::initializer::ScalarFunction<1>>(vthy);
    dict["thermal_velocity_z"] = static_cast<PHARE::initializer::ScalarFunction<1>>(vthz);

    dict["charge"]         = 1.;
    dict["nbrPartPerCell"] = int{100};
    dict["basis"]          = std::string{"Cartesian"};

    auto initializer = ParticleInitializerFactory<ParticleArrayT, GridLayoutT>::create(dict);
}



int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
