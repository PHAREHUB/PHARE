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

namespace PHARE::restarts
{
// returns a copy of dict with this build's current resources hash injected under
// simulation/restarts/resources_hash, so it flows through addRestartDict ->
// RestartsProperties::fileAttributes -> the restart file. If dict already carries a
// "file_resources_hash" (the hash read back from an old restart file being loaded),
// compares it against the current hash and throws if they differ, so an incompatible
// restart fails clearly instead of corrupting/segfaulting later.
NO_DISCARD inline initializer::PHAREDict
inject_simulation_information(initializer::PHAREDict const& dict, std::string const& resources_hash)
{
    auto injected = dict;
    if (injected["simulation"].contains("restarts"))
    {
        auto& restarts_dict             = injected["simulation"]["restarts"];
        restarts_dict["resources_hash"] = resources_hash;

        if (restarts_dict.contains("file_resources_hash"))
        {
            std::string const file_hash = restarts_dict["file_resources_hash"];
            if (!file_hash.empty() && file_hash != resources_hash)
                throw std::runtime_error(
                    "Restart file is incompatible with this PHARE build: the set of "
                    "registered resources differs (a field or particle population was "
                    "likely added, removed, or renamed since the restart file was "
                    "written). file hash=" + file_hash + " current hash=" + resources_hash);
        }
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
