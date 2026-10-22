#include "phare_mpi.hpp"

#include "amr/messengers/composite_refine_patch_strategy.hpp"

#include <SAMRAI/hier/Box.h>
#include <SAMRAI/hier/BoxContainer.h>
#include <SAMRAI/hier/Connector.h>
#include <SAMRAI/hier/Patch.h>
#include <SAMRAI/hier/PatchDescriptor.h>
#include <SAMRAI/hier/PatchLevel.h>
#include <SAMRAI/tbox/SAMRAIManager.h>
#include <SAMRAI/tbox/SAMRAI_MPI.h>

#include "gtest/gtest.h"

#include <memory>
#include <string>
#include <vector>

using namespace PHARE::amr;

namespace
{
SAMRAI::tbox::Dimension const dim{2};

// records every hook it receives, so that the forwarding and its order can be checked
class RecordingStrategy : public SAMRAI::xfer::RefinePatchStrategy
{
public:
    RecordingStrategy(std::vector<std::string>& log, std::string name, int stencil = 0,
                      bool synchronize = false)
        : log_{log}
        , name_{std::move(name)}
        , stencil_{stencil}
        , synchronize_{synchronize}
    {
    }

    void setPhysicalBoundaryConditions(SAMRAI::hier::Patch&, double const,
                                       SAMRAI::hier::IntVector const&) override
    {
        record_("bc");
    }

    SAMRAI::hier::IntVector getRefineOpStencilWidth(SAMRAI::tbox::Dimension const& d) const override
    {
        return SAMRAI::hier::IntVector{d, stencil_};
    }

    void preprocessRefine(SAMRAI::hier::Patch&, SAMRAI::hier::Patch const&,
                          SAMRAI::hier::Box const&, SAMRAI::hier::IntVector const&) override
    {
        record_("pre");
    }

    void postprocessRefine(SAMRAI::hier::Patch&, SAMRAI::hier::Patch const&,
                           SAMRAI::hier::Box const&, SAMRAI::hier::IntVector const&) override
    {
        record_("post");
    }

    void preprocessRefineBoxes(SAMRAI::hier::Patch&, SAMRAI::hier::Patch const&,
                               SAMRAI::hier::BoxContainer const&,
                               SAMRAI::hier::IntVector const&) override
    {
        record_("preBoxes");
    }

    void postprocessRefineBoxes(SAMRAI::hier::Patch&, SAMRAI::hier::Patch const&,
                                SAMRAI::hier::BoxContainer const&,
                                SAMRAI::hier::IntVector const&) override
    {
        record_("postBoxes");
    }

    void preprocessRefineLevel(
        SAMRAI::hier::PatchLevel&, SAMRAI::hier::PatchLevel const&, SAMRAI::hier::Connector const&,
        SAMRAI::hier::Connector const&,
        std::vector<std::vector<std::shared_ptr<SAMRAI::hier::BoxOverlap>>> const&,
        SAMRAI::xfer::RefineClasses::Data const**) override
    {
        record_("preLevel");
        setNeedRefineSynchronize(synchronize_);
    }

    void postprocessRefineLevel(SAMRAI::hier::PatchLevel&, SAMRAI::hier::PatchLevel const&,
                                SAMRAI::hier::Connector const&,
                                SAMRAI::hier::Connector const&) override
    {
        record_("postLevel");
        setNeedRefineSynchronize(synchronize_);
    }

private:
    void record_(std::string const& hook) const { log_.push_back(name_ + ":" + hook); }

    std::vector<std::string>& log_;
    std::string const name_;
    int const stencil_;
    bool const synchronize_;
};


struct Levels
{
    SAMRAI::hier::PatchLevel fine{dim};
    SAMRAI::hier::PatchLevel coarse{dim};
    SAMRAI::hier::Connector coarseToFine{dim};
    SAMRAI::hier::Connector coarseToUnfilled{dim};
    std::vector<std::vector<std::shared_ptr<SAMRAI::hier::BoxOverlap>>> overlaps;
};

bool preprocessLevelSynchronizes(CompositeRefinePatchStrategy& composite)
{
    Levels levels;
    composite.preprocessRefineLevel(levels.fine, levels.coarse, levels.coarseToFine,
                                    levels.coarseToUnfilled, levels.overlaps, nullptr);
    return composite.needSynchronize();
}

bool postprocessLevelSynchronizes(CompositeRefinePatchStrategy& composite)
{
    Levels levels;
    composite.postprocessRefineLevel(levels.fine, levels.coarse, levels.coarseToFine,
                                     levels.coarseToUnfilled);
    return composite.needSynchronize();
}
} // namespace



TEST(CompositeRefinePatchStrategy, forwardsEveryHookToEveryStrategyInListOrder)
{
    std::vector<std::string> log;
    CompositeRefinePatchStrategy composite{{std::make_shared<RecordingStrategy>(log, "a"),
                                            std::make_shared<RecordingStrategy>(log, "b")}};

    auto const descriptor = std::make_shared<SAMRAI::hier::PatchDescriptor>();
    SAMRAI::hier::Box const box{dim};
    SAMRAI::hier::Patch fine{box, descriptor};
    SAMRAI::hier::Patch const coarse{box, descriptor};
    SAMRAI::hier::BoxContainer const boxes;
    auto const ratio = SAMRAI::hier::IntVector{dim, 2};
    Levels levels;

    composite.setPhysicalBoundaryConditions(fine, 0., SAMRAI::hier::IntVector::getZero(dim));
    composite.preprocessRefine(fine, coarse, box, ratio);
    composite.postprocessRefine(fine, coarse, box, ratio);
    composite.preprocessRefineBoxes(fine, coarse, boxes, ratio);
    composite.postprocessRefineBoxes(fine, coarse, boxes, ratio);
    composite.preprocessRefineLevel(levels.fine, levels.coarse, levels.coarseToFine,
                                    levels.coarseToUnfilled, levels.overlaps, nullptr);
    composite.postprocessRefineLevel(levels.fine, levels.coarse, levels.coarseToFine,
                                     levels.coarseToUnfilled);

    std::vector<std::string> const expected{
        "a:bc",       "b:bc",       "a:pre",       "b:pre",       "a:post",
        "b:post",     "a:preBoxes", "b:preBoxes",  "a:postBoxes", "b:postBoxes",
        "a:preLevel", "b:preLevel", "a:postLevel", "b:postLevel"};
    EXPECT_EQ(expected, log);
}


TEST(CompositeRefinePatchStrategy, stencilWidthIsTheMaxOfTheStrategies)
{
    std::vector<std::string> log;
    CompositeRefinePatchStrategy composite{{std::make_shared<RecordingStrategy>(log, "a", 1),
                                            std::make_shared<RecordingStrategy>(log, "b", 3),
                                            std::make_shared<RecordingStrategy>(log, "c", 2)}};

    EXPECT_EQ(SAMRAI::hier::IntVector(dim, 3), composite.getRefineOpStencilWidth(dim));
}


TEST(CompositeRefinePatchStrategy, synchronizesIfAnyStrategyAsks)
{
    std::vector<std::string> log;

    CompositeRefinePatchStrategy none{{std::make_shared<RecordingStrategy>(log, "a", 0, false),
                                       std::make_shared<RecordingStrategy>(log, "b", 0, false)}};
    EXPECT_FALSE(preprocessLevelSynchronizes(none));
    EXPECT_FALSE(postprocessLevelSynchronizes(none));

    CompositeRefinePatchStrategy last{{std::make_shared<RecordingStrategy>(log, "a", 0, false),
                                       std::make_shared<RecordingStrategy>(log, "b", 0, true)}};
    EXPECT_TRUE(preprocessLevelSynchronizes(last));
    EXPECT_TRUE(postprocessLevelSynchronizes(last));

    CompositeRefinePatchStrategy first{{std::make_shared<RecordingStrategy>(log, "a", 0, true),
                                        std::make_shared<RecordingStrategy>(log, "b", 0, false)}};
    EXPECT_TRUE(preprocessLevelSynchronizes(first));
    EXPECT_TRUE(postprocessLevelSynchronizes(first));
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
