
#include "python3/pybind_def.hpp"
#include "simulator/simulator.hpp"

#include "core/def/phare_config.hpp"


#include "amr/wrappers/hierarchy.hpp" // for HierarchyRestarter::getRestartFileFullPath



namespace py = pybind11;

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
}
} // namespace PHARE::pydata
