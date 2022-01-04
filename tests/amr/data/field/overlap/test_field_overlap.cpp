#include <SAMRAI/hier/BoxContainer.h>
#include <SAMRAI/tbox/SAMRAIManager.h>
#include <SAMRAI/tbox/SAMRAI_MPI.h>

#include "gtest/gtest.h"

#include "amr/data/field/field_overlap.hpp"


using namespace PHARE::amr;

TEST(FieldOverlap, withAtLeastABoxIsNotEmpty)
{
    SAMRAI::hier::BoxContainer boxes;
    auto dim   = SAMRAI::tbox::Dimension{1};
    auto lower = SAMRAI::hier::Index(dim, 0);
    auto upper = SAMRAI::hier::Index(dim, 1);
    SAMRAI::hier::Box box{lower, upper, SAMRAI::hier::BlockId{0}};
    boxes.push_back(box);
    FieldOverlap overlap{boxes, SAMRAI::hier::Transformation{SAMRAI::hier::IntVector::getOne(dim)}};

    EXPECT_FALSE(overlap.isOverlapEmpty());
}

TEST(FieldOverlap, returnsTransformationOffset)
{
    SAMRAI::hier::BoxContainer boxes;
    auto dim          = SAMRAI::tbox::Dimension{1};
    auto lower        = SAMRAI::hier::Index(dim, 0);
    auto upper        = SAMRAI::hier::Index(dim, 1);
    auto sourceOffset = SAMRAI::hier::IntVector::getOne(dim);
    SAMRAI::hier::Box box{lower, upper, SAMRAI::hier::BlockId{0}};
    boxes.push_back(box);
    FieldOverlap overlap{boxes, SAMRAI::hier::Transformation{SAMRAI::hier::IntVector::getOne(dim)}};

    EXPECT_EQ(sourceOffset, overlap.getSourceOffset());
}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    SAMRAI::tbox::SAMRAI_MPI::init(&argc, &argv);
    SAMRAI::tbox::SAMRAIManager::initialize();
    SAMRAI::tbox::SAMRAIManager::startup();

    int testResult = RUN_ALL_TESTS();

    SAMRAI::tbox::SAMRAIManager::shutdown();
    SAMRAI::tbox::SAMRAIManager::finalize();
    SAMRAI::tbox::SAMRAI_MPI::finalize();

    return testResult;
}
