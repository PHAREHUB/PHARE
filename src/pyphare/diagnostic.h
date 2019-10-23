#ifndef PHARE_PYPHARE_DIAGNOSTIC_H
#define PHARE_PYPHARE_DIAGNOSTIC_H

#include "diagnostic_manager.h"

namespace PHARE
{
extern ADiagnosticsManager* diagnosticManager;
namespace pybind
{
    void addDiagnostic(size_t compute_every, size_t write_every, size_t start_iteration,
                       size_t end_iteration, std::string name, std::string species,
                       std::string type)
    {
        PHARE::Diagnostic diagnostic{
            compute_every, write_every, start_iteration, end_iteration, name, species, type};
        diagnosticManager->addDiagnostic(diagnostic);
    }

    template<typename PyBindModule>
    void diagnostic(PyBindModule& m)
    {
        // py::class_<PHARE::Diagnostic>(m, "Diagnostic")
        //     .def_readwrite("name", &PHARE::Diagnostic::name)
        //     .def_readwrite("species", &PHARE::Diagnostic::species)
        //     .def_readwrite("type", &PHARE::Diagnostic::type)
        //     .def_readwrite("compute_every", &PHARE::Diagnostic::compute_every)
        //     .def_readwrite("write_every", &PHARE::Diagnostic::write_every)
        //     .def_readwrite("start_iteration", &PHARE::Diagnostic::start_iteration)
        //     .def_readwrite("end_iteration", &PHARE::Diagnostic::end_iteration);

        m.def("addDiagnostic", addDiagnostic, "diagnostic");
    }

} /* namespace pybind*/
} /* namespace PHARE*/

#endif /* PHARE_PYPHARE_DIAGNOSTIC_H */
