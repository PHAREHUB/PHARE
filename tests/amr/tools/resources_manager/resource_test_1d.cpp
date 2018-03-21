#include "resource_test_1d.h"

using PHARE::HybridQuantity;
using PHARE::ResourcesManager;



void ResourcesManagerTest1D::SetUp()
{
    hierarchy = std::make_unique<BasicHierarchy>(inputBase + std::string("/input/input_db_1d"));
    hierarchy->init();

    param = GetParam();


    auto &field = *param.field;

    auto &patchHierarchy = hierarchy->hierarchy;

    resourcesManager.registerResources(field);

    for (int iLevel = 0; iLevel < patchHierarchy->getNumberOfLevels(); ++iLevel)
    {
        auto patchLevel = patchHierarchy->getPatchLevel(iLevel);
        for (auto &patch : *patchLevel)
        {
            resourcesManager.allocate(field, *patch);
        }
    }
}


TEST_P(ResourcesManagerTest1D, FieldPointerCanBeSet)
{
    auto &field          = *param.field;
    auto &patchHierarchy = hierarchy->hierarchy;

    for (int iLevel = 0; iLevel < patchHierarchy->getNumberOfLevels(); ++iLevel)
    {
        auto patchLevel = patchHierarchy->getPatchLevel(iLevel);
        for (auto const &patch : *patchLevel)
        {
            auto guards = resourcesManager.createResourcesGuards(*patch, field);
            EXPECT_TRUE(field.isValid());
        }
    }
    EXPECT_FALSE(field.isValid());
}


ResourcesManagerTest1DParam createResources(std::string const &name, HybridQuantity hq)
{
    ResourcesManagerTest1DParam rc;
    rc.field = std::make_shared<PlaceHolder::CellField>(name, hq);
    return rc;
}

INSTANTIATE_TEST_CASE_P(ResourcesManager, ResourcesManagerTest1D,
                        ::testing::ValuesIn({createResources("Ex", HybridQuantity::Ex),
                                             createResources("Ey", HybridQuantity::Ey),
                                             createResources("Ez", HybridQuantity::Ez),
                                             createResources("Bx", HybridQuantity::Bx),
                                             createResources("By", HybridQuantity::By),
                                             createResources("Bz", HybridQuantity::Bz)

                        }));
