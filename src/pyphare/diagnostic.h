#ifndef PHARE_PYPHARE_DIAGNOSTIC_H
#define PHARE_PYPHARE_DIAGNOSTIC_H

#include "diagnostic_manager.h"

namespace PHARE
{
extern IDiagnosticsManager* diagnosticManager;
namespace pybind
{
    void addDiagnostic(size_t compute_every, size_t write_every, size_t start_iteration,
                       size_t end_iteration, std::string name, std::string species,
                       std::string type)
    {
        PHARE::DiagnosticDAO diagnostic{
            compute_every, write_every, start_iteration, end_iteration, name, species, type};
        diagnosticManager->addDiagnostic(diagnostic);
    }

    template<typename PyBindModule>
    void diagnostic(PyBindModule& m)
    {
        m.def("addDiagnostic", addDiagnostic, "diagnostic");
    }

} /* namespace pybind*/
} /* namespace PHARE*/

#endif /* PHARE_PYPHARE_DIAGNOSTIC_H */
