

#include "phare/phare.hpp"
#include "simulator/simulator.hpp"
#include "amr/wrappers/hierarchy.hpp"
#include "initializer/python_data_provider.hpp"

#include "core/logger.hpp"

#include <algorithm>

namespace PHARE
{
std::unique_ptr<PHARE::ISimulator> getSimulator(std::shared_ptr<PHARE::amr::Hierarchy>& hierarchy)
{
    PHARE::initializer::PHAREDict const& theDict
        = PHARE::initializer::PHAREDictHandler::INSTANCE().dict();
    auto dim           = theDict["simulation"]["dimension"].template to<int>();
    auto interpOrder   = theDict["simulation"]["interp_order"].template to<int>();
    auto nbRefinedPart = theDict["simulation"]["refined_particle_nbr"].template to<int>();

    return core::makeAtRuntime<SimulatorMaker>(dim, interpOrder, nbRefinedPart,
                                               SimulatorMaker{hierarchy});
}
} /* namespace PHARE */


std::unique_ptr<PHARE::initializer::DataProvider> fromCommandLine(int argc, char** argv)
{
    using dataProvider [[maybe_unused]] = std::unique_ptr<PHARE::initializer::DataProvider>;

    switch (argc)
    {
        case 1: return nullptr;
        case 2:
            std::string arg = argv[1];
            auto moduleName = arg.substr(0, arg.find_last_of("."));
            if (arg.substr(arg.find_last_of(".") + 1) == "py")
            {
                std::replace(moduleName.begin(), moduleName.end(), '/', '.');
                std::cout << "python input detected, building with python provider...\n";
                return std::make_unique<PHARE::initializer::PythonDataProvider>(moduleName);
            }

            break;
    }
    return nullptr;
}

int main(int argc, char** argv)
{
    std::string const welcome = R"~(
                  _____   _    _            _____   ______
                 |  __ \ | |  | |    /\    |  __ \ |  ____|
                 | |__) || |__| |   /  \   | |__) || |__
                 |  ___/ |  __  |  / /\ \  |  _  / |  __|
                 | |     | |  | | / ____ \ | | \ \ | |____
                 |_|     |_|  |_|/_/    \_\|_|  \_\|______|)~";
    std::cout << welcome;
    std::cout << "\n";
    std::cout << "\n";


    PHARE::SamraiLifeCycle slc{argc, argv};

    std::cerr << "creating python data provider\n";
    auto provider = fromCommandLine(argc, argv);

    std::cerr << "reading user inputs...";
    provider->read();
    std::cerr << "done!\n";

    auto& dictHandler = PHARE::initializer::PHAREDictHandler::INSTANCE();

    auto hierarchy = PHARE::amr::Hierarchy::make();

    auto simulator = PHARE::getSimulator(hierarchy);

    std::cout << PHARE::core::to_str(*simulator) << "\n";

    simulator->initialize();

    dictHandler.stop();
    provider.release();

    [[maybe_unused]] auto time = simulator->startTime();

    while (simulator->currentTime() < simulator->endTime())
    {
        simulator->dump(simulator->currentTime(), simulator->timeStep());
        simulator->advance(simulator->timeStep());
        std::cout << simulator->currentTime() << "\n";
        //    time += simulator.timeStep();
    }
}
