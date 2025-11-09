#ifndef PHARE_DEFAULT_MHD_REGISTERER_HPP
#define PHARE_DEFAULT_MHD_REGISTERER_HPP

#include "python3/cpp_mhd_python_registerer.hpp"
#include "python3/mhd_defaults/mhd_defaults.hpp"

namespace PHARE::pydata
{
namespace py = pybind11;

template<typename Dimension, typename InterpOrder, typename NbRefinedPart>
class DefaultMHDRegisterer
{
    using Registerer_t = Registerer<Dimension, InterpOrder, NbRefinedPart, DefaultTimeIntegrator,
                                    DefaultReconstruction, void, DefaultRiemannSolver,
                                    DefaultEquations, false, false, false>;

    static constexpr auto dim           = Dimension{}();
    static constexpr auto interp        = InterpOrder{}();
    static constexpr auto nbRefinedPart = NbRefinedPart{}();

public:
    static inline std::string type_name = "_" + std::to_string(dim) + "_" + std::to_string(interp)
                                          + "_" + std::to_string(nbRefinedPart);

    constexpr static void declare_sim(py::module& m) { Registerer_t::declare_sim(m, type_name); }

    constexpr static void declare_defaults(py::module& m)
    {
        Registerer_t::declare_etc(m, type_name);
        declare_sim(m);
        declare_splitter(m);
    }

private:
    constexpr static void declare_splitter(py::module& m)
    {
        using _Splitter = PHARE::amr::Splitter<Dimension, InterpOrder,
                                               core::RefinedParticlesConst<nbRefinedPart>>;

        std::string name = "Splitter" + type_name;
        py::class_<_Splitter, py::smart_holder>(m, name.c_str())
            .def(py::init<>())
            .def_property_readonly_static("weight", [](py::object) { return _Splitter::weight; })
            .def_property_readonly_static("delta", [](py::object) { return _Splitter::delta; });

        name = "split_pyarray_particles" + type_name;
        m.def(name.c_str(), splitPyArrayParticles<_Splitter>);
    }
};
} // namespace PHARE::pydata

#endif // PHARE_DEFAULT_MHD_REGISTERER_HPP
