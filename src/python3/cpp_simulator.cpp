
#include "python3/cpp_simulator.hpp"

#if !defined(PHARE_CPP_MOD_NAME)
#define PHARE_CPP_MOD_NAME cpp
#endif

namespace py = pybind11;

namespace PHARE::pydata
{
PYBIND11_MODULE(PHARE_CPP_MOD_NAME, m)
{
    declare_essential(m);

    declareDim<1>(m);
    declareDim<2>(m);
    declareDim<3>(m);

    core::apply(core::possibleSimulators(), [&](auto const& simType) { declare_all(m, simType); });

    declarePatchData<std::vector<double>, 1>(m, "PatchDataVectorDouble_1D");
    declarePatchData<std::vector<double>, 2>(m, "PatchDataVectorDouble_2D");
    declarePatchData<std::vector<double>, 3>(m, "PatchDataVectorDouble_3D");
}
} // namespace PHARE::pydata

#if defined(PHARE_EXEC_PYBIND)
#define PHARE_EXE_NO_MAIN 1
#include "phare/phare.cpp"
#undef PHARE_EXE_NO_MAIN
#endif
