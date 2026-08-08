#ifndef PHARE_TESTS_AMR_TOOLS_RESSOURCE_RESSOURCE_TEST_1D_HPP
#define PHARE_TESTS_AMR_TOOLS_RESSOURCE_RESSOURCE_TEST_1D_HPP

#include "phare_solver.hpp"

#include "amr/resources_manager/resources_manager.hpp"

#include "test_resources_manager_basic_hierarchy.hpp"

#include "input_config.h"

#include "gtest/gtest.h"

#include <memory>

using namespace PHARE::core;
using namespace PHARE::amr;


template<typename ResourcesUsers>
class aResourceUserCollection : public ::testing::Test
{
public:
    std::size_t constexpr static dimension    = 1;
    std::size_t constexpr static interp_order = 1;
    auto constexpr static opts                = PHARE::SimOpts{dimension, interp_order};

    using PHARETypes         = PHARE::solver::PHARE_Types<opts>;
    using ResourcesManager_t = PHARETypes::ResourcesManager_t;

    std::unique_ptr<BasicHierarchy> hierarchy;
    ResourcesManager_t resourcesManager;

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
