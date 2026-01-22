
// PHARE_SIM_STR == template parameter string for SimOpts, e.g., "1, 2, 3"
//                  (dimension, interpolation order, refined particle count)
// PHARE_SIM_ID  == identifier string for module naming, e.g., "1_2_3"

#if !defined(PHARE_SIM_STR) || !defined(PHARE_SIM_ID)
#error "PHARE_SIM_STR and PHARE_SIM_ID must be defined! Example: PHARE_SIM_STR=\"1, 2, 3\" PHARE_SIM_ID=1_2_3"
#endif // PHARE_SIM_STR

#include "python3/cpp_simulator.hpp"

#if !defined(PHARE_CPP_MOD_NAME)
#define PHARE_CPP_MOD_NAME PHARE_STR_CAT(cpp_, PHARE_SIM_ID)
#endif


namespace PHARE::pydata
{

PYBIND11_MODULE(PHARE_CPP_MOD_NAME, m)
{
    declare_macro_sim(m);
}

} // namespace PHARE::pydata
