
#include "python3/cpp_simulator.h"

namespace PHARE::pydata
{
PYBIND11_MODULE(cpp_sim_2_1_4, m)
{
    using dim          = std::integral_constant<std::size_t, 2>;
    using interp       = std::integral_constant<std::size_t, 1>;
    using nbRefinePart = std::integral_constant<std::size_t, 4>;

    declare_essential(m);
    declare_sim<dim, interp, nbRefinePart>(m);
}
} // namespace PHARE::pydata
