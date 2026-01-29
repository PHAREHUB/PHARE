// This file is for the python module for everything besides C++ Simulators.


#include "core/def.hpp"
#include "core/def/phare_config.hpp"
#include "core/data/particles/particle_array.hpp"

#include "amr/samrai.hpp"             // SamraiLifeCycle without simulators
#include "amr/wrappers/hierarchy.hpp" // for HierarchyRestarter::getRestartFileFullPath

#include "python3/pybind_def.hpp"
#include "python3/patch_data.hpp"

#include "hdf5/phare_hdf5.hpp"

#if PHARE_HAS_HIGHFIVE
#include "hdf5/detail/h5/h5_file.hpp"
#endif

namespace py = pybind11;

namespace PHARE::pydata
{

template<typename Type, std::size_t dimension>
void declarePatchData(py::module& m, std::string key)
{
    using PatchDataType = PatchData<Type, dimension>;
    py::class_<PatchDataType>(m, key.c_str())
        .def_readonly("patchID", &PatchDataType::patchID)
        .def_readonly("origin", &PatchDataType::origin)
        .def_readonly("lower", &PatchDataType::lower)
        .def_readonly("upper", &PatchDataType::upper)
        .def_readonly("nGhosts", &PatchDataType::nGhosts)
        .def_readonly("data", &PatchDataType::data);
}

template<std::size_t dim>
void declareDim(py::module& m)
{
    using CP         = core::ContiguousParticles<dim>;
    std::string name = "ContiguousParticles_" + std::to_string(dim);
    py::class_<CP, std::shared_ptr<CP>>(m, name.c_str())
        .def(py::init<std::size_t>())
        .def_readwrite("iCell", &CP::iCell)
        .def_readwrite("delta", &CP::delta)
        .def_readwrite("weight", &CP::weight)
        .def_readwrite("charge", &CP::charge)
        .def_readwrite("v", &CP::v)
        .def("size", &CP::size);

    name = "PatchData" + name;
    declarePatchData<CP, dim>(m, name.c_str());
}

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
    auto samrai_restart_file = [](std::string path) {
        return PHARE::amr::HierarchyRestarter::getRestartFileFullPath(path);
    };
    py::class_<core::Span<double>, py::smart_holder>(m, "Span");
    py::class_<PyArrayWrapper<double>, py::smart_holder, core::Span<double>>(m, "PyWrapper");


    m.def("mpi_size", []() { return core::mpi::size(); });
    m.def("mpi_rank", []() { return core::mpi::rank(); });
    m.def("mpi_barrier", []() { core::mpi::barrier(); });

    py::class_<SamraiLifeCycle, std::shared_ptr<SamraiLifeCycle>>(m, "SamraiLifeCycle")
        .def(py::init<>())
        .def("reset", &SamraiLifeCycle::reset);

    py::class_<PHARE::amr::Hierarchy, std::shared_ptr<PHARE::amr::Hierarchy>>(m, "AMRHierarchy");
    m.def("make_hierarchy", []() { return PHARE::amr::Hierarchy::make(); });

    m.def("makePyArrayWrapper", makePyArrayWrapper<double>);

    m.def("phare_deps", []() {
        std::unordered_map<std::string, std::string> versions{{"pybind", pybind_version()},
                                                              {"samrai", samrai_version()}};
        _PHARE_WITH_HIGHFIVE(versions["highfive"] = HIGHFIVE_VERSION_STRING);
        return versions;
    });

    m.def("samrai_restart_file", samrai_restart_file);

    m.def("restart_path_for_time", [](std::string path, double timestamp) {
        return PHARE::amr::Hierarchy::restartFilePathForTime(path, timestamp);
    });

    m.def("phare_build_config", []() { return PHARE::build_config(); });

    m.def("patch_data_ids", [&](std::string const& path) -> std::vector<int> {
        _PHARE_WITH_HIGHFIVE({
            auto const& restart_file = samrai_restart_file(path);
            PHARE::hdf5::h5::HighFiveFile h5File{restart_file, HighFive::File::ReadOnly,
                                                 /*para=*/false};
            return h5File.read_data_set<int>("/phare/patch/ids");
        });

        throw std::runtime_error("PHARE not built with highfive support");
    });
    m.def("serialized_simulation_string", [&](std::string const& path) -> std::string {
        _PHARE_WITH_HIGHFIVE({
            auto const& restart_file = samrai_restart_file(path);
            PHARE::hdf5::h5::HighFiveFile h5File{restart_file, HighFive::File::ReadOnly,
                                                 /*para=*/false};
            return h5File.read_attribute("/phare", "serialized_simulation");
        });

        throw std::runtime_error("PHARE not built with highfive support");
    });


    declareDim<1>(m);
    declareDim<2>(m);
    declareDim<3>(m);

    declarePatchData<std::vector<double>, 1>(m, "PatchDataVectorDouble_1D");
    declarePatchData<std::vector<double>, 2>(m, "PatchDataVectorDouble_2D");
    declarePatchData<std::vector<double>, 3>(m, "PatchDataVectorDouble_3D");
}

} // namespace PHARE::pydata
