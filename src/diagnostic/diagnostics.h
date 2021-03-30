
#ifndef DIAGNOSTIC_DIAGNOSTICS_H
#define DIAGNOSTIC_DIAGNOSTICS_H

#include <memory>

#if !defined(PHARE_HAS_HIGHFIVE)
#error // PHARE_HAS_HIGHFIVE expected to be defined as bool
#endif

#if !defined(PHARE_DIAG_DOUBLES)
#define PHARE_DIAG_DOUBLES false
#endif

#if PHARE_HAS_HIGHFIVE
#include "highfive/H5Version.hpp"
#define _PHARE_WITH_HIGHFIVE(...) __VA_ARGS__
#else
#define _PHARE_WITH_HIGHFIVE(...)
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
    void dump(double /*timeStamp*/, double /*timeStep*/) override
    {
        throw std::runtime_error("NOOP");
    }

    void dump_level(std::size_t /*level*/, double /*timestamp*/) override
    {
        throw std::runtime_error("NOOP");
    }
};

struct DiagnosticsManagerResolver
{
    template<typename Hierarchy, typename Model>
    static std::unique_ptr<IDiagnosticsManager> make_unique(Hierarchy& hier, Model& model,
                                                            initializer::PHAREDict const& dict)
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
