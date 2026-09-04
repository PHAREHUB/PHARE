#ifndef DIAGNOSTIC_DIAGNOSTICS_HPP
#define DIAGNOSTIC_DIAGNOSTICS_HPP

#include "core/def.hpp"

#if !defined(PHARE_HAS_HIGHFIVE)
#error // PHARE_HAS_HIGHFIVE expected to be defined as bool
#endif

#if !defined(PHARE_DIAG_DOUBLES)
#define PHARE_DIAG_DOUBLES false
#endif

#include "diagnostic_manager.hpp"

#if PHARE_HAS_HIGHFIVE
#include "diagnostic/detail/h5writer.hpp"
#include "diagnostic/detail/vtkh5_writer.hpp"
#endif // PHARE_HAS_HIGHFIVE

#include "dict.hpp"


#include <memory>

namespace PHARE::diagnostic
{
struct NullOpDiagnosticsManager : public IDiagnosticsManager
{
    bool dump(double /*timeStamp*/, double /*timeStep*/) override
    {
        throw std::runtime_error("NOOP");
    }

    bool dump(double /*timeStamp*/) override { throw std::runtime_error("NOOP"); }

    void dump_level(std::size_t /*level*/, double /*timestamp*/) override
    {
        throw std::runtime_error("NOOP");
    }
};

struct DiagnosticsManagerResolver
{
    template<typename Hierarchy, typename... Models>
    NO_DISCARD static std::unique_ptr<IDiagnosticsManager>
    make_unique(Hierarchy& hier, initializer::PHAREDict const& dict, Models&... models)
    {
#if PHARE_HAS_HIGHFIVE

        auto const format = cppdict::get_value(dict, "format", std::string{"phareh5"});

        if (format == "phareh5")
            return DiagnosticsManager<
                h5::H5Writer<DiagnosticsModelMapper<Hierarchy, Models...>>>::make_unique(hier, dict,
                                                                                         models...);
        if (format == "pharevtkhdf")
            return DiagnosticsManager<vtkh5::H5Writer<
                DiagnosticsModelMapper<Hierarchy, Models...>>>::make_unique(hier, dict, models...);
        throw std::runtime_error("DiagnosticsManagerResolver - unknown format " + format);
#else
        return std::make_unique<NullOpDiagnosticsManager>();
#endif
    }
};

} // namespace PHARE::diagnostic

#endif // DIAGNOSTIC_DIAGNOSTICS_HPP
