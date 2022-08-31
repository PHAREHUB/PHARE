

#include "core/utilities/types.hpp"
#include "core/utilities/span.hpp"
#include "core/data/grid/gridlayout.hpp"
#include "core/data/grid/gridlayoutimplyee.hpp"
#include "core/data/ions/particle_initializers/particle_initializer_factory.hpp"
#include "core/data/particles/particle.hpp"
#include "core/data/particles/particle_array.hpp"
#include "initializer/data_provider.hpp"


#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <vector>
#include <type_traits>

using namespace PHARE::core;
using namespace PHARE::initializer;

using namespace PHARE::core;
using namespace PHARE::initializer;

#include "tests/initializer/init_functions.hpp"
using namespace PHARE::initializer::test_fn::func_1d; // density/etc are here


using GridLayoutT    = GridLayout<GridLayoutImplYee<1, 1>>;
using ParticleArrayT = ParticleArray<Particle<1>>;


TEST(AParticleIinitializerFactory, takesAPHAREDictToCreateAParticleVectorInitializer)
{
    PHARE::initializer::PHAREDict dict;
    dict["name"]            = std::string{"MaxwellianParticleInitializer"};
    dict["fn_type"]         = std::string{"vector"};
    dict["density"]         = static_cast<PHARE::initializer::InitFunction<1>>(density);
    dict["bulk_velocity_x"] = static_cast<PHARE::initializer::InitFunction<1>>(vx);
    dict["bulk_velocity_y"] = static_cast<PHARE::initializer::InitFunction<1>>(vx);
    dict["bulk_velocity_z"] = static_cast<PHARE::initializer::InitFunction<1>>(vx);

    dict["thermal_velocity_x"] = static_cast<PHARE::initializer::InitFunction<1>>(vthx);
    dict["thermal_velocity_y"] = static_cast<PHARE::initializer::InitFunction<1>>(vthy);
    dict["thermal_velocity_z"] = static_cast<PHARE::initializer::InitFunction<1>>(vthz);

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
