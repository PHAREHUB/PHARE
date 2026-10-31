#ifndef PHARE_PYTHON_CPP_SIMULATOR_OPTS_HPP
#define PHARE_PYTHON_CPP_SIMULATOR_OPTS_HPP

// PHARE_SIM_STR == template string, between <>
// shared by every translation unit of one cpp_<PHARE_SIM_ID> module

#ifndef PHARE_SIM_STR
// mostly for clangformat - errors in cpp file if define is missing
#define PHARE_SIM_STR 1ull, 1ull, 2ull
#endif

#include "simulator/simulator_def.hpp"

namespace PHARE::pydata
{
auto static constexpr resolve_simulator_options()
{
    using namespace PHARE;
    using namespace PHARE::core;
    using namespace PHARE::MHDOpts;
    return SimOpts{PHARE_SIM_STR};
}

auto inline constexpr simulator_options = resolve_simulator_options();

} // namespace PHARE::pydata

#endif /*PHARE_PYTHON_CPP_SIMULATOR_OPTS_HPP*/
