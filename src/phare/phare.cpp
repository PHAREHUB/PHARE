
#include "data/grid/gridlayout.h"
#include "data/grid/gridlayoutimplyee.h"
#include "data/particles/particle_array.h"
#include "models/physical_state.h"
#include "python_data_provider.h"

#include <iostream>



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

    auto provider = fromCommandLine(argc, argv);
    provider->read();
}
