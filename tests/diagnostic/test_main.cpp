
#include <filesystem>

#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5Easy.hpp>

#include "diagnostic/include.h"
#include "diagnostic/samrai_lifecycle.h"

#include "diagnostic/func.h"
#include "diagnostic/defs.h"
using namespace PHARE_test::_1d;

#define PHARE_WITH_HIGHFIVE
#include "diagnostic/samrai_diagnostic.h"

#include "diagnostic/integrator.h"
#include "diagnostic/tag_strat.h"
#include "diagnostic/hierarchy.h"

struct Hi5Diagnostic : public AfullHybridBasicHierarchy
{
    HighFive::File file{"/tmp/new_diagnostic.file.h5", HighFive::File::ReadWrite
                                                           | HighFive::File::Create
                                                           | HighFive::File::Truncate};
};

TEST_F(Hi5Diagnostic, initializesFieldsOnRefinedLevels)
{
    using Samraive       = PHARE::SamraiHighFiveDiagnostic<HybridModelT>;
    constexpr auto level = PHARE::diagnostic::Level::DBG;
    constexpr auto mode  = PHARE::diagnostic::Mode::FULL;
    PHARE::hi5::Diagnostic hi5{this->file, mode};
    Samraive samhighfo{basicHierarchy->getHierarchy(), *hybridModel, hi5};
    PHARE::DiagnosticsManager<Samraive, level> dMan{samhighfo};
    dMan.addDiagDict().addDiagDict().addDiagDict().addDiagDict();
    dMan.dump();
}

int main(int argc, char** argv)
{
    PHARE_test::SamraiLifeCycle samrai_lifecycle(argc, argv);
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
