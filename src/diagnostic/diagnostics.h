
#ifndef DIAGNOSTIC_DIAGNOSTICS_H
#define DIAGNOSTIC_DIAGNOSTICS_H

#if !defined(PHARE_HAS_HIGHFIVE)
#pragma message PHARE_HAS_HIGHFIVE undefined, disabled.
#define PHARE_HAS_HIGHFIVE 0
#endif

#include "diagnostic_manager.h"

#include "cppdict/include/dict.hpp"

#if PHARE_HAS_HIGHFIVE

#include "diagnostic_model_view.h"

#include "diagnostic/detail/h5writer.h"
#include "diagnostic/detail/types/electromag.h"
#include "diagnostic/detail/types/particle.h"
#include "diagnostic/detail/types/fluid.h"

#endif

namespace PHARE::diagnostic
{
struct NullOpDiagnosticsManager : public IDiagnosticsManager
{
    void dump(double timeStamp, double timeStep) override { throw std::runtime_error("NOOP"); }
};

struct DiagnosticsManagerResolver
{
    template<typename Hierarchy, typename Model>
    static decltype(auto) make_shared(Hierarchy& hier, Model& model, initializer::PHAREDict& dict)
    {
#if PHARE_HAS_HIGHFIVE
        using ModelView_t = ModelView<Hierarchy, Model>;
        using Writer_t    = h5::Writer<ModelView_t>;
        return DiagnosticsManager<Writer_t>::make_unique(hier, model, dict);
#else
        return std::make_unique<NullOpDiagnosticsManager>();
#endif
    }
};

} // namespace PHARE::diagnostic

#endif // DIAGNOSTIC_DIAGNOSTICS_H
