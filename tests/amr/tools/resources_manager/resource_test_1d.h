#ifndef PHARE_TESTS_AMR_TOOLS_RESSOURCE_RESSOURCE_TEST_1D_H
#define PHARE_TESTS_AMR_TOOLS_RESSOURCE_RESSOURCE_TEST_1D_H

#include <memory>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "basic_hierarchy.h"
#include "cell_placeholder.h"
#include "data/grid/gridlayout.h"
#include "data/grid/gridlayout_impl.h"
#include "data/vecfield/vecfield.h"
#include "input_config.h"
#include "tools/resources_manager.h"

using namespace PHARE;

struct ResourcesManagerTest1DParam
{
    std::shared_ptr<VecField<NdArrayVector1D<>, HybridQuantity>> vecfield;
};

class ResourcesManagerTest1D : public ::testing::TestWithParam<ResourcesManagerTest1DParam>
{
public:
    ResourcesManagerTest1D() = default;

    virtual void SetUp() override;
    virtual void TearDown() override {}

    std::unique_ptr<BasicHierarchy> hierarchy;
    ResourcesManagerTest1DParam param;

    ResourcesManager<GridLayout<GridLayoutImplYee<1, 1>>> resourcesManager{
        SAMRAI::tbox::Dimension{1}};
};


#endif
