
#include <iostream>

#include "python_data_provider.h"

#include "diagnostic/include.h"
#include "diagnostic/samrai_lifecycle.h"

#include "diagnostic/func.h"
#include "diagnostic/defs.h"
using namespace PHARE_test::_1d;

#include "diagnostic/detail/samrai_highfive.h"

#include "diagnostic/integrator.h"
#include "diagnostic/tag_strat.h"
#include "diagnostic/hierarchy.h"

namespace PHARE
{
std::unique_ptr<PHARE::initializer::DataProvider> fromCommandLine(int argc, char** argv)
{
    using dataProvider [[maybe_unused]] = std::unique_ptr<PHARE::initializer::DataProvider>;

    switch (argc)
    {
        case 1: return nullptr;
        case 2:
            std::string arg = argv[1];
            if (arg.substr(arg.find_last_of(".") + 1) == "py")
            {
                std::cout << "python input detected, building with python provider...\n";
                auto provider
                    = std::make_unique<PHARE::initializer::PythonDataProvider>(argc, argv[1]);
                return provider;
            }

            break;
    }
    return nullptr;
}

class GlobalLifeCycle
{
    using Writer = PHARE::SamraiHighFiveDiagnostic<HybridModelT>;

public:
    GlobalLifeCycle(std::string path = "/tmp/roflcopter.hi5")
        : hi5{path}
        , samhighfo{fullHybrid.basicHierarchy->getHierarchy(), *fullHybrid.hybridModel, hi5}
        , dMan{new PHARE::DiagnosticsManager<Writer>(samhighfo)}
    {
        PHARE::initializer::PHAREDictHandler::INSTANCE().init<1>();
    }
    ~GlobalLifeCycle() { PHARE::initializer::PHAREDictHandler::INSTANCE().stop<1>(); }

    AfullHybridBasicHierarchy fullHybrid;
    PHARE::hi5::Diagnostic hi5;
    PHARE::SamraiHighFiveDiagnostic<HybridModelT> samhighfo;
    std::unique_ptr<PHARE::ADiagnosticsManager> dMan;
};

// global extern for python addition of diagnostics
PHARE::ADiagnosticsManager* diagnosticManager;
} /*namespace PHARE*/

int main(int argc, char** argv)
{
    std::string const welcome = R"~(
                  _____   _    _            _____   ______
                 |  __ \ | |  | |    /\    |  __ \ |  ____|
                 | |__) || |__| |   /  \   | |__) || |__
                 |  ___/ |  __  |  / /\ \  |  _  / |  __|
                 | |     | |  | | / ____ \ | | \ \ | |____
                 |_|     |_|  |_|/_/    \_\|_|  \_\|______|)~";

    std::cout << welcome << "\n\n";
    PHARE_test::SamraiLifeCycle samrai{argc, argv};
    auto provider = PHARE::fromCommandLine(argc, argv);
    if (!provider)
    {
        std::cerr << "Invalid input detected" << std::endl;
        return 1;
    }
    PHARE::GlobalLifeCycle lifeCycle;
    PHARE::diagnosticManager = lifeCycle.dMan.get();
    provider->read();
}
