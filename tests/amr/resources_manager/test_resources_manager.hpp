#ifndef PHARE_TESTS_AMR_TOOLS_RESSOURCE_RESSOURCE_TEST_1D_HPP
#define PHARE_TESTS_AMR_TOOLS_RESSOURCE_RESSOURCE_TEST_1D_HPP



#include "phare_core.hpp"

#include "test_resources_manager_basic_hierarchy.hpp"
#include "core/data/grid/grid.hpp"
#include "core/data/grid/gridlayout.hpp"
#include "core/data/grid/gridlayout_impl.hpp"
#include "core/data/ions/ion_population/ion_population.hpp"
#include "core/data/ions/ions.hpp"
#include "core/data/ndarray/ndarray_vector.hpp"
#include "core/data/particles/particle_array.hpp"
#include "core/data/vecfield/vecfield.hpp"
#include "input_config.h"
#include "amr/resources_manager/resources_manager.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <memory>

using namespace PHARE::core;
using namespace PHARE::amr;




template<typename ResourcesUsers>
class aResourceUserCollection : public ::testing::Test
{
public:
    using core_Types   = PHARE::core::PHARE_Types<1, 1>;
    using Grid_t       = core_Types::Grid_t;
    using GridLayout_t = core_Types::GridLayout_t;

    std::unique_ptr<BasicHierarchy> hierarchy;
    ResourcesManager<GridLayout_t, Grid_t> resourcesManager;

    ResourcesUsers users;

    void SetUp()
    {
        auto s    = inputBase + std::string("/input/input_db_1d");
        hierarchy = std::make_unique<BasicHierarchy>(inputBase + std::string("/input/input_db_1d"));
        hierarchy->init();

        auto registerAndAllocate = [this](auto& resourcesUser) {
            auto& patchHierarchy = hierarchy->hierarchy;

            resourcesManager.registerResources(resourcesUser.user);

            double const initDataTime{0.0};

            for (int iLevel = 0; iLevel < patchHierarchy->getNumberOfLevels(); ++iLevel)
            {
                auto patchLevel = patchHierarchy->getPatchLevel(iLevel);
                for (auto& patch : *patchLevel)
                {
                    resourcesManager.allocate(resourcesUser.user, *patch, initDataTime);
                }
            }
        }; // end lambda

        std::apply(registerAndAllocate, users);
    }
};




#endif
