
#include "python3/cpp_simulator.hpp"
#include "core/utilities/meta/meta_utilities.hpp"

#if !defined(PHARE_CPP_MOD_NAME)
#define PHARE_CPP_MOD_NAME cpp
#endif

namespace PHARE::pydata
{
PYBIND11_MODULE(PHARE_CPP_MOD_NAME, m)
{
    declare_essential(m);

    declareDim<1>(m);
    declareDim<2>(m);
    declareDim<3>(m);

    core::apply(core::possibleSimulators(), [&](auto const& simType) { declare_all(m, simType); });
}
} // namespace PHARE::pydata
