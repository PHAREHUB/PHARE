#ifndef PHARE_PYPHARE_DIAGNOSTIC_H
#define PHARE_PYPHARE_DIAGNOSTIC_H

#include "diagnostic_manager.h"

namespace PHARE
{
extern ADiagnosticsManager* diagnosticManager;
namespace pybind
{
    void dict(std::string type)
    {
        PHARE::initializer::PHAREDict<1> dict;
        dict["diag"]["name"]            = type;
        dict["diag"]["type"]            = type;
        dict["diag"]["species"]         = std::string{"SPECIES_" + type};
        dict["diag"]["compute_every"]   = std::size_t{1};
        dict["diag"]["write_every"]     = std::size_t{1};
        dict["diag"]["start_iteration"] = std::size_t{0};
        dict["diag"]["end_iteration"]   = std::numeric_limits<std::size_t>::max();
        diagnosticManager->addDiagDict(dict);
    }

    template<typename PyBindModule>
    void diagnostic(PyBindModule& m)
    {
        m.def("diagnostic", dict, "diagnostic");
    }

} /* namespace pybind*/
} /* namespace PHARE*/

#endif /* PHARE_PYPHARE_DIAGNOSTIC_H */
