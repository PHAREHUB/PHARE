
#include "python3/cpp_simulator.hpp"

#if !defined(PHARE_SIM_STR) || !defined(PHARE_SIM_ID)
#error // needs PHARE_SIM_STR! like "1, 1, 1"
#endif // PHARE_SIM_STR

#if !defined(PHARE_CPP_MOD_NAME)
#define PHARE_CPP_MOD_NAME PHARE_STR_CAT(cpp_, PHARE_SIM_ID)
#endif

namespace py = pybind11;

namespace PHARE::pydata
{
PYBIND11_MODULE(PHARE_CPP_MOD_NAME, m)
{
    declare_macro_sim(m);
}
} // namespace PHARE::pydata
