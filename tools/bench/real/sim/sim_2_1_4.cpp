#include "python3/cpp_simulator.hpp"
#include "python3/mhd_defaults/default_mhd_registerer.hpp"

namespace PHARE::pydata
{
PYBIND11_MODULE(cpp_sim_2_1_4, m)
{
    using dim          = std::integral_constant<std::size_t, 2>;
    using interp       = std::integral_constant<std::size_t, 1>;
    using nbRefinePart = std::integral_constant<std::size_t, 4>;

    declare_essential(m);
    DefaultMHDRegisterer<dim, interp, nbRefinePart>::declare_sim(m);
}
} // namespace PHARE::pydata
