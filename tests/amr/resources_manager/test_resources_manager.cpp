#include "test_resources_manager.hpp"
#include "core/data/electromag/electromag.hpp"
#include "core/data/electrons/electrons.hpp"
#include "core/data/grid/grid.hpp"
#include "core/data/grid/gridlayout.hpp"
#include "core/data/ions/particle_initializers/maxwellian_particle_initializer.hpp"
#include "initializer/data_provider.hpp"
#include "core/models/hybrid_state.hpp"
#include "core/data/vecfield/vecfield.hpp"
#include "core/data/tensorfield/tensorfield.hpp"
#include "core/utilities/variants.hpp"

#include <core/hybrid/hybrid_quantities.hpp>
#include <core/utilities/types.hpp>

#include <SAMRAI/tbox/SAMRAIManager.h>
#include <SAMRAI/tbox/SAMRAI_MPI.h>

#include "gtest/gtest.h"

#include <string>
#include <vector>


#include "tests/initializer/init_functions.hpp"
using namespace PHARE::initializer::test_fn::func_1d; // density/etc are here

static constexpr std::size_t dim         = 1;
static constexpr std::size_t interpOrder = 1;
using GridImplYee1D                      = GridLayoutImplYee<dim, interpOrder>;
using GridYee1D                          = GridLayout<GridImplYee1D>;

using GridImplYee1D    = GridLayoutImplYee<dim, interpOrder>;
using GridYee1D        = GridLayout<GridImplYee1D>;
using Field1D          = Field<dim, HybridQuantity::Scalar>;
using Grid1D           = Grid<NdArrayVector<dim>, HybridQuantity::Scalar>;
using VecField1D       = VecField<Field1D, HybridQuantity>;
using SymTensorField1D = SymTensorField<Field1D, HybridQuantity>;
using IonPopulation1D  = IonPopulation<ParticleArray<1>, VecField1D, SymTensorField1D>;
using Ions1D           = Ions<IonPopulation1D, GridYee1D>;
using Electromag1D     = Electromag<VecField1D>;
using Electrons1D      = Electrons<Ions1D>;
using HybridState1D    = HybridState<Electromag1D, Ions1D, Electrons1D>;

using MaxwellianParticleInitializer1D
    = MaxwellianParticleInitializer<ParticleArray<dim>, GridYee1D>;

using InitFunctionT = PHARE::initializer::InitFunction<dim>;

PHARE::initializer::PHAREDict createInitDict()
{
    PHARE::initializer::PHAREDict dict;
    dict["ions"]["nbrPopulations"] = std::size_t{2};
    dict["ions"]["pop0"]["name"]   = std::string{"protons"};
    dict["ions"]["pop0"]["mass"]   = 1.;
    dict["ions"]["pop0"]["particle_initializer"]["name"]
        = std::string{"MaxwellianParticleInitializer"};
    dict["ions"]["pop0"]["particle_initializer"]["density"] = static_cast<InitFunctionT>(density);

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
    dict["ions"]["pop1"]["particle_initializer"]["density"] = static_cast<InitFunctionT>(density);

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

    dict["electromag"]["name"]             = std::string{"EM"};
    dict["electromag"]["electric"]["name"] = std::string{"E"};
    dict["electromag"]["magnetic"]["name"] = std::string{"B"};

    dict["electromag"]["magnetic"]["initializer"]["x_component"] = static_cast<InitFunctionT>(bx);
    dict["electromag"]["magnetic"]["initializer"]["y_component"] = static_cast<InitFunctionT>(by);
    dict["electromag"]["magnetic"]["initializer"]["z_component"] = static_cast<InitFunctionT>(bz);

    dict["electrons"]["pressure_closure"]["name"] = std::string{"isothermal"};
    dict["electrons"]["pressure_closure"]["Te"]   = 0.12;


    return dict;
}
static auto init_dict = createInitDict();

struct IonPopulation1D_P
{
    IonPopulation1D user{init_dict["ions"]["pop0"]};
};


struct VecField1D_P
{
    std::string name = "B";
    HybridQuantity::Vector qty{HybridQuantity::Vector::B};
    VecField1D user{name, qty};
};




struct Ions1D_P
{
    Ions1D user{init_dict["ions"]};
};




struct Electromag1D_P
{
    std::string name = "ElectroTest";
    Electromag1D user;
    Electromag1D_P()
        : user{init_dict["electromag"]}
    {
    }
};



struct HybridState1D_P
{
    HybridState1D user{init_dict};
};




using IonPop1DOnly          = std::tuple<IonPopulation1D_P>;
using VecField1DOnly        = std::tuple<VecField1D_P>;
using Ions1DOnly            = std::tuple<Ions1D_P>;
using VecField1DAndIonPop1D = std::tuple<VecField1D_P, IonPopulation1D_P>;
using Electromag1DOnly      = std::tuple<Electromag1D_P>;
using HybridState1DOnly     = std::tuple<HybridState1D_P>;

TYPED_TEST_SUITE_P(aResourceUserCollection);




TYPED_TEST_P(aResourceUserCollection, hasPointersValidOnlyWithGuard)
{
    TypeParam resourceUserCollection;

    auto check = [this](auto& resourceUserPack) {
        auto& hierarchy_   = this->hierarchy->hierarchy;
        auto& resourceUser = resourceUserPack.user;

        for (int iLevel = 0; iLevel < hierarchy_->getNumberOfLevels(); ++iLevel)
        {
            auto patchLevel = hierarchy_->getPatchLevel(iLevel);
            for (auto const& patch : *patchLevel)
            {
                auto dataOnPatch = this->resourcesManager.setOnPatch(*patch, resourceUser);
                EXPECT_TRUE(resourceUser.isUsable());
                EXPECT_FALSE(resourceUser.isSettable());
            }
            EXPECT_FALSE(resourceUser.isUsable());
            EXPECT_TRUE(resourceUser.isSettable());
        }
    };

    std::apply(check, resourceUserCollection);
}



TYPED_TEST_P(aResourceUserCollection, hasPointersValidWithEnumerate)
{
    TypeParam resourceUserCollection;

    auto check = [this](auto& resourceUserPack) {
        auto& hierarchy_   = this->hierarchy->hierarchy;
        auto& resourceUser = resourceUserPack.user;
        auto& rm           = this->resourcesManager;
        for (int iLevel = 0; iLevel < hierarchy_->getNumberOfLevels(); ++iLevel)
        {
            auto patchLevel = hierarchy_->getPatchLevel(iLevel);
            for (auto const& patch : rm.enumerate(*patchLevel, resourceUser))
            {
                EXPECT_TRUE(resourceUser.isUsable());
                EXPECT_FALSE(resourceUser.isSettable());
            }
            EXPECT_FALSE(resourceUser.isUsable());
            EXPECT_TRUE(resourceUser.isSettable());
        }
    };

    std::apply(check, resourceUserCollection);
}


TEST(usingResourcesManager, toGetTimeOfAResourcesUser)
{
    std::unique_ptr<BasicHierarchy> hierarchy;
    ResourcesManager<GridLayout<GridLayoutImplYee<1, 1>>, Grid1D> resourcesManager;
    IonPopulation1D_P pop;
    static_assert(is_particles_v<ParticlesPack<ParticleArray<1>>>);

    auto s    = inputBase + std::string("/input/input_db_1d");
    hierarchy = std::make_unique<BasicHierarchy>(inputBase + std::string("/input/input_db_1d"));
    hierarchy->init();
    resourcesManager.registerResources(pop.user);
    auto& patchHierarchy = hierarchy->hierarchy;

    double const initDataTime{3.14};

    for (int iLevel = 0; iLevel < patchHierarchy->getNumberOfLevels(); ++iLevel)
    {
        auto patchLevel = patchHierarchy->getPatchLevel(iLevel);
        for (auto& patch : *patchLevel)
        {
            resourcesManager.allocate(pop.user, *patch, initDataTime);
            auto times = resourcesManager.getTimes(pop.user, *patch);

            EXPECT_TRUE(std::equal(std::begin(times) + 1, std::end(times), std::begin(times)));
            EXPECT_DOUBLE_EQ(initDataTime, times[0]);
        }
    }
}




REGISTER_TYPED_TEST_SUITE_P(aResourceUserCollection, hasPointersValidOnlyWithGuard,
                            hasPointersValidWithEnumerate);


typedef ::testing::Types<IonPop1DOnly, VecField1DOnly, Ions1DOnly, Electromag1DOnly,
                         HybridState1DOnly>
    MyTypes;
INSTANTIATE_TYPED_TEST_SUITE_P(testResourcesManager, aResourceUserCollection, MyTypes);




struct FieldResource
{
    NO_DISCARD auto getCompileTimeResourcesViewList() { return std::forward_as_tuple(rho); }

    std::array<std::uint32_t, 1> cells{5};
    Field1D rho{"name", HybridQuantity::Scalar::rho, nullptr, cells};
};
struct VecFieldResource
{
    NO_DISCARD auto getCompileTimeResourcesViewList() { return std::forward_as_tuple(B); }

    VecField1D B{"B", HybridQuantity::Vector::B};
};

struct ResourceUser
{
    using Resources = std::variant<FieldResource>;

    ResourceUser() { resources.resize(2, FieldResource{}); }

    NO_DISCARD std::vector<Resources>& getRunTimeResourcesViewList() { return resources; }

    auto get() { return get_as_tuple_or_throw<FieldResource, FieldResource>(resources); }

    std::vector<Resources> resources;
};


TEST(usingResources, test_variants_helpers)
{
    using Resources = std::variant<FieldResource, VecFieldResource>;

    std::vector<Resources> resources{FieldResource{}, VecFieldResource{}};

    {
        auto const&& [rho, B] = get_as_tuple_or_throw<FieldResource, VecFieldResource>(resources);
    }
    {
        auto& rho = get_as_ref_or_throw<FieldResource>(resources);
        EXPECT_EQ(&rho, &std::get<FieldResource>(resources[0]));
    }
    {
        auto& B = get_as_ref_or_throw<VecFieldResource>(resources);
        EXPECT_EQ(&B, &std::get<VecFieldResource>(resources[1]));
    }
}


TEST(usingResources, test_variants_resource_helpers)
{
    using Resources = std::variant<Field1D, VecField1D>;

    std::array<std::uint32_t, 1> cells{5};
    Field1D moe_{"moe", HybridQuantity::Scalar::rho, nullptr, cells};
    VecField1D B_{"B", HybridQuantity::Vector::B};
    Field1D rho_{"rho", HybridQuantity::Scalar::rho, nullptr, cells};
    VecField1D E_{"E", HybridQuantity::Vector::E};

    std::vector<Resources> resources{rho_, moe_, B_, E_};
    {
        auto const& rho = get_from_variants(resources, rho_);
        EXPECT_EQ(rho.name(), "rho");

        auto const& B = get_from_variants(resources, B_);
        EXPECT_EQ(B.name(), "B");
    }
    auto [rho, B] = get_from_variants(resources, rho_, B_);
    EXPECT_EQ(rho.name(), "rho");
    EXPECT_EQ(B.name(), "B");
}

TEST(usingResourcesManager, test_variants)
{
    ResourceUser resourceUser;

    std::unique_ptr<BasicHierarchy> bhierarchy;
    ResourcesManager<GridLayout<GridLayoutImplYee<1, 1>>, Grid1D> resourcesManager;

    bhierarchy = std::make_unique<BasicHierarchy>(inputBase + std::string("/input/input_db_1d"));
    bhierarchy->init();
    resourcesManager.registerResources(resourceUser);
    auto& hierarchy = bhierarchy->hierarchy;

    for (int iLevel = 0; iLevel < hierarchy->getNumberOfLevels(); ++iLevel)
        for (auto const& patch : *hierarchy->getPatchLevel(iLevel))
            resourcesManager.allocate(resourceUser, *patch, 0);

    std::size_t patches = 0, checks = 0;
    for (int iLevel = 0; iLevel < hierarchy->getNumberOfLevels(); ++iLevel)
        for (auto const& patch : *hierarchy->getPatchLevel(iLevel))
        {
            auto dataOnPatch = resourcesManager.setOnPatch(*patch, resourceUser);
            auto&& [r0, r1]  = resourceUser.get();
            r1.rho.data()[4] = 5;
            ++patches;
        }


    for (int iLevel = 0; iLevel < hierarchy->getNumberOfLevels(); ++iLevel)
        for (auto const& patch : *hierarchy->getPatchLevel(iLevel))
        {
            auto dataOnPatch = resourcesManager.setOnPatch(*patch, resourceUser);
            auto&& [r0, r1]  = resourceUser.get();
            EXPECT_EQ(r1.rho.data()[4], 5);
            ++checks;
        }

    EXPECT_GE(patches, 1);
    EXPECT_EQ(checks, patches);
}




int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    SAMRAI::tbox::SAMRAI_MPI::init(&argc, &argv);
    SAMRAI::tbox::SAMRAIManager::initialize();
    SAMRAI::tbox::SAMRAIManager::startup();


    int testResult = RUN_ALL_TESTS();

    // Finalize
    SAMRAI::tbox::SAMRAIManager::shutdown();
    SAMRAI::tbox::SAMRAIManager::finalize();
    SAMRAI::tbox::SAMRAI_MPI::finalize();

    return testResult;
}
