#include <type_traits>



#include "core/data/ions/ion_population/ion_population.hpp"
#include "core/data/ions/ions.hpp"
#include "core/data/ndarray/ndarray_vector.hpp"
#include "core/data/particles/particle_array.hpp"
#include "core/data/vecfield/vecfield.hpp"
#include "core/hybrid/hybrid_quantities.hpp"

#include "core/data/tensorfield//tensorfield.hpp"
#include "core/data/grid/gridlayout.hpp"
#include "core/data/grid/gridlayout_impl.hpp"
#include "core/data/ions/particle_initializers/maxwellian_particle_initializer.hpp"
#include "initializer/data_provider.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "tests/initializer/init_functions.hpp"
using namespace PHARE::initializer::test_fn::func_1d; // density/etc are here

using namespace PHARE::core;

static constexpr std::size_t dim         = 1;
static constexpr std::size_t interpOrder = 1;
using GridImplYee1D                      = GridLayoutImplYee<dim, interpOrder>;
using GridYee1D                          = GridLayout<GridImplYee1D>;
using MaxwellianParticleInitializer1D = MaxwellianParticleInitializer<ParticleArray<1>, GridYee1D>;



class theIons : public ::testing::Test
{
protected:
    using VecField1D       = VecField<NdArrayVector<1>, HybridQuantity>;
    using SymTensorField1D = SymTensorField<NdArrayVector<1>, HybridQuantity>;
    using InitFunctionT    = PHARE::initializer::InitFunction<1>;

    using IonPopulation1D
        = IonPopulation<ParticleArray<1>, VecField1D, SymTensorField1D, GridYee1D>;
    Ions<IonPopulation1D, GridYee1D> ions;

    PHARE::initializer::PHAREDict createIonsDict()
    {
        PHARE::initializer::PHAREDict dict;
        dict["ions"]["nbrPopulations"] = int{2};
        dict["ions"]["pop0"]["name"]   = std::string{"protons"};
        dict["ions"]["pop0"]["mass"]   = 1.;
        dict["ions"]["pop0"]["particle_initializer"]["name"]
            = std::string{"MaxwellianParticleInitializer"};
        dict["ions"]["pop0"]["particle_initializer"]["density"]
            = static_cast<InitFunctionT>(density);

        dict["ions"]["pop0"]["particle_initializer"]["bulk_velocity_x"]
            = static_cast<InitFunctionT>(vx);

        dict["ions"]["pop0"]["particle_initializer"]["bulk_velocity_y"]
            = static_cast<InitFunctionT>(vy);

        dict["ions"]["pop0"]["particle_initializer"]["bulk_velocity_z"]
            = static_cast<InitFunctionT>(vz);


        dict["ions"]["pop0"]["particle_initializer"]["thermal_velocity_x"]
            = static_cast<InitFunctionT>(vthx);

        dict["ions"]["pop0"]["particle_initializer"]["thermal_velocity_y"]
            = static_cast<InitFunctionT>(vthy);

        dict["ions"]["pop0"]["particle_initializer"]["thermal_velocity_z"]
            = static_cast<InitFunctionT>(vthz);


        dict["ions"]["pop0"]["particle_initializer"]["nbrPartPerCell"] = int{100};
        dict["ions"]["pop0"]["particle_initializer"]["charge"]         = -1.;
        dict["ions"]["pop0"]["particle_initializer"]["basis"]          = std::string{"Cartesian"};

        dict["ions"]["pop1"]["name"] = std::string{"alpha"};
        dict["ions"]["pop1"]["mass"] = 1.;
        dict["ions"]["pop1"]["particle_initializer"]["name"]
            = std::string{"MaxwellianParticleInitializer"};
        dict["ions"]["pop1"]["particle_initializer"]["density"]
            = static_cast<InitFunctionT>(density);

        dict["ions"]["pop1"]["particle_initializer"]["bulk_velocity_x"]
            = static_cast<InitFunctionT>(vx);

        dict["ions"]["pop1"]["particle_initializer"]["bulk_velocity_y"]
            = static_cast<InitFunctionT>(vy);

        dict["ions"]["pop1"]["particle_initializer"]["bulk_velocity_z"]
            = static_cast<InitFunctionT>(vz);


        dict["ions"]["pop1"]["particle_initializer"]["thermal_velocity_x"]
            = static_cast<InitFunctionT>(vthx);

        dict["ions"]["pop1"]["particle_initializer"]["thermal_velocity_y"]
            = static_cast<InitFunctionT>(vthy);

        dict["ions"]["pop1"]["particle_initializer"]["thermal_velocity_z"]
            = static_cast<InitFunctionT>(vthz);


        dict["ions"]["pop1"]["particle_initializer"]["nbrPartPerCell"] = int{100};
        dict["ions"]["pop1"]["particle_initializer"]["charge"]         = -1.;
        dict["ions"]["pop1"]["particle_initializer"]["basis"]          = std::string{"Cartesian"};

        return dict;
    }


    theIons()
        : ions{createIonsDict()["ions"]}
    {
    }

public:
    ~theIons();
};

theIons::~theIons() {}


TEST_F(theIons, areAContainerOfIonPopulations)
{
    //
    for (auto& pop : ions)
    {
        (void)pop;
    }
}




TEST_F(theIons, areNotUsableUponConstruction)
{
    EXPECT_FALSE(ions.isUsable());
}




TEST_F(theIons, areSettableUponConstruction)
{
    EXPECT_TRUE(ions.isSettable());
}




TEST_F(theIons, throwIfAccessingDensityWhileNotUsable)
{
    EXPECT_ANY_THROW(auto& n = ions.density());
}



int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
