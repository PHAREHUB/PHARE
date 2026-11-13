#ifndef RESTART_RESTARTS_HPP
#define RESTART_RESTARTS_HPP

#include "core/def.hpp"

#if !defined(PHARE_HAS_HIGHFIVE)
#error // PHARE_HAS_HIGHFIVE expected to be defined as bool
#endif

#include "restarts_manager.hpp"

#if PHARE_HAS_HIGHFIVE
#include "restarts/detail/h5writer.hpp"
#endif

#include <memory>
#include <string>
#include <vector>
#include <iterator>
#include <algorithm>

namespace PHARE::restarts
{
namespace detail
{
    // what `a` has that `b` has not, both sorted
    inline std::string only_in(std::vector<std::string> const& a, std::vector<std::string> const& b)
    {
        std::vector<std::string> diff;
        std::set_difference(a.begin(), a.end(), b.begin(), b.end(), std::back_inserter(diff));
        std::string out;
        for (auto const& name : diff)
            out += "\n    " + name;
        return out.empty() ? "\n    (none)" : out;
    }

    // empty file_names: the restart file predates this check
    inline void check_resource_names(std::vector<std::string> file_names,
                                     std::vector<std::string> const& names)
    {
        std::sort(file_names.begin(), file_names.end());
        if (file_names.empty() or file_names == names)
            return;

        throw std::runtime_error(
            "Restart file is incompatible with this PHARE build, the registered resources differ"
            "\n  only in the restart file:"
            + only_in(file_names, names) + "\n  only in this build:" + only_in(names, file_names));
    }
} // namespace detail


// returns a copy of dict with this build's registered resource names under
// simulation/restarts/resource_names, so they flow through addRestartDict ->
// RestartsProperties::fileAttributes -> the restart file. If dict carries the names read back
// from a restart file being loaded ("file_resource_names"), throws if they differ, naming what
// each side has that the other has not, so an incompatible restart fails clearly instead of
// corrupting/segfaulting later.
NO_DISCARD inline initializer::PHAREDict
inject_simulation_information(initializer::PHAREDict const& dict,
                              std::vector<std::string> const& resource_names)
{
    auto injected = dict;
    if (injected["simulation"].contains("restarts"))
    {
        auto& restarts_dict             = injected["simulation"]["restarts"];
        restarts_dict["resource_names"] = resource_names;

        if (restarts_dict.contains("file_resource_names"))
            detail::check_resource_names(
                restarts_dict["file_resource_names"].to<std::vector<std::string>>(),
                resource_names);
    }
    return injected;
}


struct NullOpRestartsManager : public IRestartsManager
{
    bool dump(double /*timeStamp*/, double /*timeStep*/) override
    {
        throw std::runtime_error("NOOP");
    }
};

struct RestartsManagerResolver
{
    template<typename Hierarchy, typename ResourceManager_t>
    NO_DISCARD static std::unique_ptr<IRestartsManager>
    make_unique(Hierarchy& hier, ResourceManager_t& resman, initializer::PHAREDict const& dict)
    {
#if PHARE_HAS_HIGHFIVE
        using Writer_t = h5::Writer<Hierarchy, ResourceManager_t>;
        return RestartsManager<Writer_t>::make_unique(hier, resman, dict);
#else
        return std::make_unique<NullOpRestartsManager>();
#endif
    }
};

} // namespace PHARE::restarts

#endif // RESTART_RESTARTS_HPP
