#include "phare/include.h"

#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <pybind11/stl.h>


namespace py = pybind11;

namespace PHARE::pydata
{
PYBIND11_MODULE(cpp, m)
{
    StaticSamraiLifeCycle::INSTANCE(); // init

    py::class_<PHARE::amr::Hierarchy, std::shared_ptr<PHARE::amr::Hierarchy>>(m, "AMRHierarchy");

    py::class_<ISimulator, std::shared_ptr<ISimulator>>(m, "ISimulator")
        .def("initialize", &PHARE::ISimulator::initialize)
        .def("advance", &PHARE::ISimulator::advance)
        .def("startTime", &PHARE::ISimulator::startTime)
        .def("currentTime", &PHARE::ISimulator::currentTime)
        .def("endTime", &PHARE::ISimulator::endTime)
        .def("timeStep", &PHARE::ISimulator::timeStep)
        .def("to_str", &PHARE::ISimulator::to_str);

    py::class_<RuntimeDiagnosticInterface, std::shared_ptr<RuntimeDiagnosticInterface>>(
        m, "IDiagnosticsManager")
        .def("dump", &RuntimeDiagnosticInterface::dump);

    m.def("make_hierarchy", []() { return PHARE::amr::Hierarchy::make(); });
    m.def("make_simulator", [](std::shared_ptr<PHARE::amr::Hierarchy>& hier) {
        auto sim = PHARE::getSimulator(hier);
        auto ptr = sim.get();
        sim.release();
        return std::shared_ptr<ISimulator>{ptr};
    });
    m.def("make_diagnostic_manager", [](std::shared_ptr<ISimulator> const& sim,
                                        std::shared_ptr<PHARE::amr::Hierarchy> const& hier) {
        return std::make_shared<RuntimeDiagnosticInterface>(*sim, *hier);
    });

    m.def("mpi_size", []() {
        int mpi_size;
        MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
        return mpi_size;
    });
    m.def("mpi_rank", []() {
        int mpi_rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
        return mpi_rank;
    });
    m.def("reset", []() {
        py::gil_scoped_release release;
        StaticSamraiLifeCycle::reset();
    });
}

} // namespace PHARE::pydata

#if defined(PHARE_DEV_MODE) && PHARE_DEV_MODE
#include "kul/io.hpp"
static void on_library_load() __attribute__((constructor));
void on_library_load()
{
    kul::Dir build{
#include "build.dir"
    };
    kul::File file{"src/python3/cpp.timestamp", build};
    std::string ts{kul::io::Reader{file}.readLine()};

    uint64_t timestamp = {
#include "cpp.timestamp"
    };
    if (std::to_string(timestamp) != ts)
        std::cerr << __FILE__ << " BUILD POSSIBLY STALE" << std::endl;
}
#endif /*PHARE_DEV_MODE*/
