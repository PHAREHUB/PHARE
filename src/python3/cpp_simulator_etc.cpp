// one part of the cpp_<PHARE_SIM_ID> module, see simulator/simulator_decl.hpp
// split to bound compiler memory per translation unit

#if !defined(PHARE_SIM_STR) || !defined(PHARE_SIM_ID)
#error // needs PHARE_SIM_STR! like "1, 1, 1"
#endif // PHARE_SIM_STR

#include "python3/cpp_simulator_etc.hpp"

namespace PHARE::pydata
{
template void declare_etc<simulator_options>(py::module&);
} // namespace PHARE::pydata
