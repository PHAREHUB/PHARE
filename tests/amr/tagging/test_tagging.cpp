


#include "amr/tagging/tagger.h"
#include "amr/tagging/tagger_factory.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

using namespace PHARE::amr;


TEST(test_tagger, fromFactory)
{
    auto hybridTagger = TaggerFactory::make("hybridModel", "grad");
    // EXPECT_TRUE(hybridTagger != nullptr);

    auto badTagger = TaggerFactory::make("invalidModel", "grad");
    // EXPECT_TRUE(badTagger == nullptr);
}




int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    int testResult = RUN_ALL_TESTS();
    return testResult;
}
