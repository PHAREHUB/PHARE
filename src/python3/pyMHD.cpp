#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "core/utilities/types.hpp"
#include "tests/core/numerics/mock_mhd_simulator/test_mhd_solver.hpp"

namespace py = pybind11;


template<std::size_t Constant>
using DimConst = PHARE::core::DimConst<Constant>;

template<std::size_t Constant>
using InterpConst = PHARE::core::InterpConst<Constant>;

template<typename DimConstant, typename InterpConstant>
using MHDMockSimulatorOption = std::tuple<DimConstant, InterpConstant>;

constexpr decltype(auto) possibleMHDMockSimulators()
{
    return std::tuple<MHDMockSimulatorOption<DimConst<1>, InterpConst<1>>,
                      MHDMockSimulatorOption<DimConst<2>, InterpConst<1>>,
                      MHDMockSimulatorOption<DimConst<3>, InterpConst<1>>>{};
}

template<std::size_t dim, std::size_t ord>
std::unique_ptr<MHDMockSimulator<dim, ord>> makeMHDMockSimulator()
{
    return std::make_unique<MHDMockSimulator<dim, ord>>(
        PHARE::initializer::PHAREDictHandler::INSTANCE().dict());
}

template<typename _dim, typename _ord>
void declare_mhd_mock_sim(py::module& m, std::tuple<_dim, _ord> const&)
{
    constexpr auto dim = _dim{}();
    constexpr auto ord = _ord{}();

    std::string type_string = "_" + std::to_string(dim) + "_" + std::to_string(ord);

    using Sim        = MHDMockSimulator<dim, ord>;
    std::string name = "MHDMockSimulator" + type_string;
    declareMHDMockSimulator<Sim>(
        py::class_<Sim, std::shared_ptr<Sim>>(m, name.c_str())
            .def_property_readonly_static("dims", [](py::object) { return Sim::dimension; })
            .def_property_readonly_static("order", [](py::object) { return Sim::order; }));

    name = "make_mhd_mock_simulator" + type_string;
    m.def(name.c_str(),
          []() { return std::shared_ptr<Sim>{std::move(makeMHDMockSimulator<dim, ord>())}; });
}

template<typename Simulator, typename PyClass>
void declareMHDMockSimulator(PyClass&& sim)
{
    sim.def("advance", &Simulator::advance, py::arg("filename"), py::arg("dumpfrequency"));
}


PYBIND11_MODULE(pyMHD, m)
{
    PHARE::core::apply(possibleMHDMockSimulators(),
                       [&](auto const& simType) { declare_mhd_mock_sim(m, simType); });
}
