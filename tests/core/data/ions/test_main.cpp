
#include <type_traits>



#include "data/ions/ion_population/ion_population.h"
#include "data/ions/ions.h"
#include "data/ndarray/ndarray_vector.h"
#include "data/particles/particle_array.h"
#include "data/vecfield/vecfield.h"
#include "hybrid/hybrid_quantities.h"


#include "data/grid/gridlayout.h"
#include "data/grid/gridlayout_impl.h"
#include "data/ions/particle_initializers/maxwellian_particle_initializer.h"
#include "data_provider.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"


using namespace PHARE::core;

static constexpr std::size_t dim         = 1;
static constexpr std::size_t interpOrder = 1;
using GridImplYee1D                      = GridLayoutImplYee<dim, interpOrder>;
using GridYee1D                          = GridLayout<GridImplYee1D>;
using MaxwellianParticleInitializer1D = MaxwellianParticleInitializer<ParticleArray<1>, GridYee1D>;



double density(double x)
{
    return x * x + 2.;
}

std::array<double, 3> bulkVelocity(double x)
{
    return std::array<double, 3>{{1.0, 0.0, 0.0}};
}


std::array<double, 3> thermalVelocity(double x)
{
    return std::array<double, 3>{{0.5, 0.0, 0.0}};
}




class theIons : public ::testing::Test
{
protected:
    using VecField1D     = VecField<NdArrayVector1D<>, HybridQuantity>;
    using ScalarFunction = PHARE::initializer::ScalarFunction<1>;
    using VectorFunction = PHARE::initializer::VectorFunction<1>;

    using IonPopulation1D = IonPopulation<ParticleArray<1>, VecField1D, GridYee1D>;
    Ions<IonPopulation1D, GridYee1D> ions;

    PHARE::initializer::PHAREDict<1> createIonsDict()
    {
        PHARE::initializer::PHAREDict<1> dict;
        dict["name"]                                = std::string{"ions"};
        dict["nbrPopulations"]                      = std::size_t{2};
        dict["pop0"]["name"]                        = std::string{"protons"};
        dict["pop0"]["mass"]                        = 1.;
        dict["pop0"]["ParticleInitializer"]["name"] = std::string{"MaxwellianParticleInitializer"};
        dict["pop0"]["ParticleInitializer"]["density"] = static_cast<ScalarFunction>(density);

        dict["pop0"]["ParticleInitializer"]["bulkVelocity"]
            = static_cast<VectorFunction>(bulkVelocity);

        dict["pop0"]["ParticleInitializer"]["thermalVelocity"]
            = static_cast<VectorFunction>(thermalVelocity);

        dict["pop0"]["ParticleInitializer"]["nbrPartPerCell"] = std::size_t{100};
        dict["pop0"]["ParticleInitializer"]["charge"]         = -1.;
        dict["pop0"]["ParticleInitializer"]["basis"]          = std::string{"Cartesian"};



        dict["pop1"]["name"]                        = std::string{"protons"};
        dict["pop1"]["mass"]                        = 1.;
        dict["pop1"]["ParticleInitializer"]["name"] = std::string{"MaxwellianParticleInitializer"};
        dict["pop1"]["ParticleInitializer"]["density"] = static_cast<ScalarFunction>(density);

        dict["pop1"]["ParticleInitializer"]["bulkVelocity"]
            = static_cast<VectorFunction>(bulkVelocity);

        dict["pop1"]["ParticleInitializer"]["thermalVelocity"]
            = static_cast<VectorFunction>(thermalVelocity);

        dict["pop1"]["ParticleInitializer"]["nbrPartPerCell"] = std::size_t{100};
        dict["pop1"]["ParticleInitializer"]["charge"]         = -1.;
        dict["pop1"]["ParticleInitializer"]["basis"]          = std::string{"Cartesian"};

        return dict;
    }


    theIons()
        : ions{createIonsDict()}
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
    EXPECT_ANY_THROW(ions.density());
}



int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
