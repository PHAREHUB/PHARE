
#include "pybind11/stl.h"
#include "pybind11/stl_bind.h"
#include "python3/pybind_def.hpp"
#include "simulator/simulator.hpp"

#include "core/env.hpp"
#include "core/def/phare_config.hpp"

#include "amr/wrappers/hierarchy.hpp" // for HierarchyRestarter::getRestartFileFullPath

namespace py = pybind11;

PYBIND11_MAKE_OPAQUE(std::unordered_map<std::string, PHARE::env::Var*>);

namespace PHARE::pydata
{
auto pybind_version()
{
    std::stringstream ss;
    ss << PYBIND11_VERSION_MAJOR << ".";
    ss << PYBIND11_VERSION_MINOR << ".";
    ss << PHARE_TO_STR(PYBIND11_VERSION_PATCH);
    return ss.str();
}

auto samrai_version()
{
    std::stringstream ss;
    ss << SAMRAI_VERSION_MAJOR << ".";
    ss << SAMRAI_VERSION_MINOR << ".";
    ss << SAMRAI_VERSION_PATCHLEVEL;
    return ss.str();
}

PYBIND11_MODULE(cpp_etc, m)
{
    py::class_<core::Span<double>, std::shared_ptr<core::Span<double>>>(m, "Span");
    py::class_<PyArrayWrapper<double>, std::shared_ptr<PyArrayWrapper<double>>, core::Span<double>>(
        m, "PyWrapper");

    m.def("makePyArrayWrapper", makePyArrayWrapper<double>);

    m.def("phare_deps", []() {
        std::unordered_map<std::string, std::string> versions{{"pybind", pybind_version()},
                                                              {"samrai", samrai_version()}};
        _PHARE_WITH_HIGHFIVE(versions["highfive"] = HIGHFIVE_VERSION_STRING);
        return versions;
    });

    m.def("samrai_restart_file", [](std::string path) {
        return PHARE::amr::HierarchyRestarter::getRestartFileFullPath(path);
    });

    m.def("restart_path_for_time", [](std::string path, double timestamp) {
        return PHARE::amr::Hierarchy::restartFilePathForTime(path, timestamp);
    });

    m.def("phare_build_config", []() { return PHARE::build_config(); });

    m.def("phare_env_exists",
          [](std::string const& s) { return Env::INSTANCE().vars.count(s) > 0; });
    m.def("phare_env_val", [](std::string const& s) { return Env::INSTANCE()(s)(); });
    py::class_<env::Var>(m, "phare_env_var")
        .def_readonly("id", &env::Var::id)
        .def_readonly("desc", &env::Var::desc)
        .def_readonly("options", &env::Var::options)
        .def_readonly("default", &env::Var::_default)
        .def_readonly("results", &env::Var::results);

    py::bind_map<std::unordered_map<std::string, env::Var*>>(m, "EnvVarMap");

    m.def(
        "phare_env_vars", []() -> auto& { return Env::INSTANCE().vars; },
        py::return_value_policy::reference);
}
} // namespace PHARE::pydata
