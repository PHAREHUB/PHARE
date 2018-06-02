#include "resource_test_1d.h"



void ResourcesManagerTest1D::SetUp()
{
    auto s    = inputBase + std::string("/input/input_db_1d");
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
            auto guards = resourcesManager.makeResourcesGuard(*patch, field);
            EXPECT_TRUE(field.isValid());
        }
    }
    EXPECT_FALSE(field.isValid());
}


ResourcesManagerTest1DParam createResources(std::string const &name, HybridQuantity::Scalar hq)
{
    ResourcesManagerTest1DParam rc;
    rc.field = std::make_shared<PlaceHolder::CellField>(name, hq);
    return rc;
}

INSTANTIATE_TEST_CASE_P(ResourcesManager, ResourcesManagerTest1D,
                        ::testing::ValuesIn({createResources("Ex", HybridQuantity::Scalar::Ex),
                                             createResources("Ey", HybridQuantity::Scalar::Ey),
                                             createResources("Ez", HybridQuantity::Scalar::Ez),
                                             createResources("Bx", HybridQuantity::Scalar::Bx),
                                             createResources("By", HybridQuantity::Scalar::By),
                                             createResources("Bz", HybridQuantity::Scalar::Bz)

                        }));
