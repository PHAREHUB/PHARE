#ifndef PHARE_TESTS_AMR_TOOLS_RESSOURCE_RESSOURCE_TEST_1D_H
#define PHARE_TESTS_AMR_TOOLS_RESSOURCE_RESSOURCE_TEST_1D_H

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "basic_hierarchy.h"
#include "cell_placeholder.h"
#include "tools/resources_manager.h"

#include <memory>

extern std::string inputBase;

struct ResourcesManagerTest1DParam
{
    std::shared_ptr<PlaceHolder::CellField> field;
};

class ResourcesManagerTest1D : public ::testing::TestWithParam<ResourcesManagerTest1DParam>
{
public:
    ResourcesManagerTest1D() = default;

    virtual void SetUp() override;
    virtual void TearDown() override {}

    std::unique_ptr<BasicHierarchy> hierarchy;
    ResourcesManagerTest1DParam param;

    PHARE::ResourcesManager resourcesManager{"yee", SAMRAI::tbox::Dimension{1}};
};
#endif
