// This file is for the python module for everything besides C++ Simulators.


#include "core/def.hpp"
#include "core/def/phare_config.hpp" // IWYU pragma: keep
#include "core/data/particles/particle_array.hpp"
#include "core/data/particles/particle_array_def.hpp"

#include "amr/samrai.hpp"             // SamraiLifeCycle without simulators
#include "amr/wrappers/hierarchy.hpp" // for HierarchyRestarter::getRestartFileFullPath


#include "cpp_simulator.hpp"

#include "python3/pybind_def.hpp"
#include "python3/patch_data.hpp"

#include "hdf5/phare_hdf5.hpp"

#if PHARE_HAS_HIGHFIVE
#include "hdf5/detail/h5/h5_file.hpp"
#endif


#include <pybind11/stl_bind.h>
#include <pybind11/native_enum.h>


namespace py = pybind11;

namespace PHARE::pydata
{



template<std::size_t dim>
void declareDim(py::module& m)
{
    using Particle   = core::Particle<dim>;
    std::string name = "Particle_" + std::to_string(dim);
    py::class_<Particle, std::shared_ptr<Particle>>(m, name.c_str())
        .def_readonly("iCell", &Particle::iCell_)
        .def_readonly("delta", &Particle::delta_)
        .def_readonly("v", &Particle::v_);

    using CP = core::SoAParticleArray<dim>;
    name     = "ParticleArray_SOA_" + std::to_string(dim);
    py::class_<CP, py::smart_holder>(m, name.c_str())
        .def(py::init<std::size_t>())
        .def_readwrite("iCell", &CP::iCell_)
        .def_readwrite("delta", &CP::delta_)
        .def_readwrite("weight", &CP::weight_)
        .def_readwrite("charge", &CP::charge_)
        .def_readwrite("v", &CP::v_)
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

auto supported_layouts()
{
    using enum core::LayoutMode;
    std::vector layouts{AoSMapped};
    PHARE_WITH_MKN_GPU({
        layouts.emplace_back(AoSTS);
        layouts.emplace_back(AoSPCTS);
    })
    return layouts;
}

PYBIND11_MODULE(cpp_etc, m, py::mod_gil_not_used())
{
    auto samrai_restart_file = [](std::string path) {
        return PHARE::amr::HierarchyRestarter::getRestartFileFullPath(path);
    };

    py::class_<core::Span<double>, py::smart_holder>(m, "Span");
    m.def("makeSpan", makePySpan<double>);
    py::class_<PyArrayWrapper<double>, py::smart_holder, core::Span<double>>(m, "PyWrapper");


    m.def("mpi_size", []() { return core::mpi::size(); });
    m.def("mpi_rank", []() { return core::mpi::rank(); });
    m.def("mpi_barrier", []() { core::mpi::barrier(); });
    m.def("mpi_initialized", []() { core::mpi::is_init(); });

    py::class_<SamraiLifeCycle, std::shared_ptr<SamraiLifeCycle>>(m, "SamraiLifeCycle")
        .def(py::init<>())
        .def("reset", &SamraiLifeCycle::reset);

    py::class_<PHARE::amr::Hierarchy, py::smart_holder>(m, "AMRHierarchy");
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

    using enum core::LayoutMode;
    py::native_enum<core::LayoutMode>(m, "LayoutMode", "enum.Enum")
        .value("AoSMapped", AoSMapped)
        .value("AoSTS", AoSTS)
        .value("AoSPCTS", AoSPCTS)
        .export_values()
        .finalize();


    m.def("supported_layouts", supported_layouts);

    declareDim<1>(m);
    declareDim<2>(m);
    declareDim<3>(m);

    // declarePatchData<py_array_t<double>, 1>(m, "PatchPyArrayDouble_1D");
    // declarePatchData<py_array_t<double>, 2>(m, "PatchPyArrayDouble_2D");
    // declarePatchData<py_array_t<double>, 3>(m, "PatchPyArrayDouble_3D");
}

} // namespace PHARE::pydata
