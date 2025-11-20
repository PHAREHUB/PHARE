#ifndef DIAGNOSTIC_DIAGNOSTICS_HPP
#define DIAGNOSTIC_DIAGNOSTICS_HPP

#include "core/def.hpp"
#include <memory>

#if !defined(PHARE_HAS_HIGHFIVE)
#error // PHARE_HAS_HIGHFIVE expected to be defined as bool
#endif

#if !defined(PHARE_DIAG_DOUBLES)
#define PHARE_DIAG_DOUBLES false
#endif


#include "hdf5/phare_hdf5.hpp"


#include "diagnostic_manager.hpp"

#include "cppdict/include/dict.hpp"

#if PHARE_HAS_HIGHFIVE

#include "diagnostic_model_view.hpp"

#include "diagnostic/detail/h5writer.hpp"
#include "diagnostic/detail/vtkh5_writer.hpp"

#endif

namespace PHARE::diagnostic
{
struct NullOpDiagnosticsManager : public IDiagnosticsManager
{
    bool dump(double /*timeStamp*/, double /*timeStep*/) override
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
    NO_DISCARD static std::unique_ptr<IDiagnosticsManager>
    make_unique(Hierarchy& hier, Model& model, initializer::PHAREDict const& dict)
    {
#if PHARE_HAS_HIGHFIVE
        using ModelView_t = ModelView<Hierarchy, Model>;
        auto const format = cppdict::get_value(dict, "format", std::string{"phareh5"});

        if (format == "phareh5")
            return DiagnosticsManager<h5::H5Writer<ModelView_t>>::make_unique(hier, model, dict);
        if (format == "pharevtkhdf")
            return DiagnosticsManager<vtkh5::H5Writer<ModelView_t>>::make_unique(hier, model, dict);
        throw std::runtime_error("DiagnosticsManagerResolver - unknown format " + format);
#else
        return std::make_unique<NullOpDiagnosticsManager>();
#endif
    }
};

} // namespace PHARE::diagnostic

#endif // DIAGNOSTIC_DIAGNOSTICS_HPP
