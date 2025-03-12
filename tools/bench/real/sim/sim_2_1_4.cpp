

#include "python3/cpp_simulator.hpp"

namespace PHARE::pydata
{
template<typename Dimension, typename InterpOrder, typename NbRefinedPart>
constexpr void declare_default_sim(py::module& m)
{
    using DefaultRegisterer
        = Registerer<Dimension, InterpOrder, NbRefinedPart, TimeIntegratorType::TVDRK3,
                     ReconstructionType::WENOZ, SlopeLimiterType::count, RiemannSolverType::Rusanov,
                     false, false, false>;

    std::string type_name = "_" + std::to_string(Dimension{}()) + "_"
                            + std::to_string(InterpOrder{}()) + "_"
                            + std::to_string(NbRefinedPart{}()) + "_" + "tvdrk3_wenoz_rusanov";

    DefaultRegisterer::declare_sim(m, type_name);
}

PYBIND11_MODULE(cpp_sim_2_1_4, m)
{
    using dim          = std::integral_constant<std::size_t, 2>;
    using interp       = std::integral_constant<std::size_t, 1>;
    using nbRefinePart = std::integral_constant<std::size_t, 4>;

    declare_essential(m);
    declare_default_sim<dim, interp, nbRefinePart>(m);
}
} // namespace PHARE::pydata
