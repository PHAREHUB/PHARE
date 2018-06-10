#ifndef PHARE_TESTS_AMR_TOOLS_RESSOURCE_RESSOURCE_TEST_1D_H
#define PHARE_TESTS_AMR_TOOLS_RESSOURCE_RESSOURCE_TEST_1D_H

#include <memory>


#include "basic_hierarchy.h"
#include "data/field/field.h"
#include "data/grid/gridlayout.h"
#include "data/grid/gridlayout_impl.h"
#include "data/ions/ion_population/ion_population.h"
#include "data/ions/ions.h"
#include "data/ndarray/ndarray_vector.h"
#include "data/particles/particle_array.h"
#include "data/vecfield/vecfield.h"
#include "input_config.h"
#include "tools/resources_manager.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

using namespace PHARE;




template<typename ResourcesUsers>
class aResourceUserCollection : public ::testing::Test
{
public:
    std::unique_ptr<BasicHierarchy> hierarchy;
    ResourcesManager<GridLayout<GridLayoutImplYee<1, 1>>> resourcesManager{
        SAMRAI::tbox::Dimension{1}};

    ResourcesUsers users;

    void SetUp()
    {
        auto s    = inputBase + std::string("/input/input_db_1d");
        hierarchy = std::make_unique<BasicHierarchy>(inputBase + std::string("/input/input_db_1d"));
        hierarchy->init();

        auto registerAndAllocate = [this](auto &resourcesUser) {
            auto &patchHierarchy = hierarchy->hierarchy;

            resourcesManager.registerResources(resourcesUser.user);

            for (int iLevel = 0; iLevel < patchHierarchy->getNumberOfLevels(); ++iLevel)
            {
                auto patchLevel = patchHierarchy->getPatchLevel(iLevel);
                for (auto &patch : *patchLevel)
                {
                    resourcesManager.allocate(resourcesUser.user, *patch);
                }
            }
        }; // end lambda

        std::apply(registerAndAllocate, users);
    }
};




#endif
