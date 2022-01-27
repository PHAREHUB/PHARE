
#include "python3/cpp_simulator.hpp"

namespace PHARE::pydata
{
PYBIND11_MODULE(cpp_sim_2_3_4, m)
{
    using dim          = std::integral_constant<std::size_t, 2>;
    using interp       = std::integral_constant<std::size_t, 3>;
    using nbRefinePart = std::integral_constant<std::size_t, 4>;

    declare_essential(m);
    declare_sim<dim, interp, nbRefinePart>(m);
}
} // namespace PHARE::pydata
