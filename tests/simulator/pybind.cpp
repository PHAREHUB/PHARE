
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>

#include "phare/include.h"

namespace py = pybind11;

namespace PHARE
{
class StreamAppender : public SAMRAI::tbox::Logger::Appender
{
public:
    StreamAppender(std::ostream* stream) { d_stream = stream; }
    void logMessage(const std::string& message, const std::string& filename, const int line)
    {
        (*d_stream) << "At :" << filename << " line :" << line << " message: " << message
                    << std::endl;
    }

private:
    std::ostream* d_stream;
};

class SamraiLifeCycle
{
public:
    static SamraiLifeCycle& INSTANCE()
    {
        static SamraiLifeCycle i;
        return i;
    }

    SamraiLifeCycle()
    {
        SAMRAI::tbox::SAMRAI_MPI::init(0, nullptr);
        SAMRAI::tbox::SAMRAIManager::initialize();
        SAMRAI::tbox::SAMRAIManager::startup();

        std::shared_ptr<SAMRAI::tbox::Logger::Appender> appender
            = std::make_shared<StreamAppender>(StreamAppender{&std::cout});
        SAMRAI::tbox::Logger::getInstance()->setWarningAppender(appender);
    }
    ~SamraiLifeCycle()
    {
        SAMRAI::tbox::SAMRAIManager::shutdown();
        SAMRAI::tbox::SAMRAIManager::finalize();
        SAMRAI::tbox::SAMRAI_MPI::finalize();
    }

    void reset()
    {
        PHARE::initializer::PHAREDictHandler::INSTANCE().stop();
        if (sim_)
            delete sim_;
        sim_ = nullptr;
        if (rdi_)
            rdi_.reset();
        SAMRAI::tbox::SAMRAIManager::shutdown();
        SAMRAI::tbox::SAMRAIManager::startup();
    }

    std::unique_ptr<ISimulator, py::nodelete> getISimulator()
    {
        auto simulation = PHARE::getSimulator();
        sim_            = simulation.get();
        std::unique_ptr<ISimulator, py::nodelete> ptr{sim_};
        simulation.release();
        return ptr;
    }

    std::unique_ptr<diagnostic::IDiagnosticsManager, py::nodelete> getIDiagnosticsManager()
    {
        if (sim_ == nullptr)
            throw std::runtime_error("getDiagnosticsManager requires simulator be created first");
        rdi_ = std::make_unique<RuntimeDiagnosticInterface>(*sim_);
        std::unique_ptr<diagnostic::IDiagnosticsManager, py::nodelete> ptr{rdi_->dMan.get()};
        return ptr;
    }

private:
    ISimulator* sim_ = nullptr;
    std::unique_ptr<RuntimeDiagnosticInterface> rdi_;
};

PYBIND11_MODULE(test_simulator, m)
{
    SamraiLifeCycle::INSTANCE(); // init
    py::class_<ISimulator, std::unique_ptr<ISimulator, py::nodelete>>(m, "ISimulator")
        .def("initialize", &PHARE::ISimulator::initialize)
        .def("advance", &PHARE::ISimulator::advance)
        .def("startTime", &PHARE::ISimulator::startTime)
        .def("endTime", &PHARE::ISimulator::endTime)
        .def("timeStep", &PHARE::ISimulator::timeStep)
        .def("to_str", &PHARE::ISimulator::to_str);

    py::class_<diagnostic::IDiagnosticsManager,
               std::unique_ptr<diagnostic::IDiagnosticsManager, py::nodelete>>(
        m, "IDiagnosticsManager")
        .def("dump", &PHARE::diagnostic::IDiagnosticsManager::dump);

    m.def("make_simulator", []() { return SamraiLifeCycle::INSTANCE().getISimulator(); });
    m.def("make_diagnostic_manager",
          []() { return SamraiLifeCycle::INSTANCE().getIDiagnosticsManager(); });
    m.def("unmake_simulator", []() {
        py::gil_scoped_release release;
        SamraiLifeCycle::INSTANCE().reset();
    });
}
} // namespace PHARE
