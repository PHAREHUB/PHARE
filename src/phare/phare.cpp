

#include "include.h"

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
                return std::make_unique<PHARE::initializer::PythonDataProvider>(argc, argv[1]);
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

    SamraiLifeCycle slc{argc, argv};

    std::cerr << "creating python data provider\n";
    auto provider = std::make_unique<PHARE::initializer::PythonDataProvider>(2, "phare.init");

    std::cerr << "reading user inputs...";
    provider->read();
    std::cerr << "done!\n";

    auto simulator = PHARE::getSimulator();

    std::cout << PHARE::core::to_str(*simulator) << "\n";

    simulator->initialize();
    RuntimeDiagnosticInterface{*simulator}.dump();

    //
    // auto time = simulator.startTime();
    //
    // while (time < simulator.endTime())
    //{
    //    simulator.advance();
    //    time += simulator.timeStep();
    //}
}
