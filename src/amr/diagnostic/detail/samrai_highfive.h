#ifndef PHARE_AMR_DIAGNOSTIC_SAMRAI_HIGHFIVE_H
#define PHARE_AMR_DIAGNOSTIC_SAMRAI_HIGHFIVE_H

#include "phare/diagnostic_manager.h"

#if defined(PHARE_WITH_HIGHFIVE)
#include "phare/hi5/diagnostic.h"
#endif

namespace PHARE
{
template<typename MODEL>
class SamraiHighFiveDiagnostic : public PHARE::amr::SamraiDiagnostic<MODEL>
{
public:
    using H = SAMRAI::hier::PatchHierarchy;

    SamraiHighFiveDiagnostic(H& h, MODEL& r, hi5::Diagnostic& hifive)
        : PHARE::amr::SamraiDiagnostic<MODEL>(h, r)
        , hi5_(hifive)
    {
    }

    template<diagnostic::Level L>
    void dump(std::vector<std::shared_ptr<Diagnostic>> const&);

private:
    hi5::Diagnostic& hi5_;

    SamraiHighFiveDiagnostic(const SamraiHighFiveDiagnostic&)             = delete;
    SamraiHighFiveDiagnostic(const SamraiHighFiveDiagnostic&&)            = delete;
    SamraiHighFiveDiagnostic& operator&(const SamraiHighFiveDiagnostic&)  = delete;
    SamraiHighFiveDiagnostic& operator&(const SamraiHighFiveDiagnostic&&) = delete;
};

template<typename MODEL>
template<diagnostic::Level L>
void SamraiHighFiveDiagnostic<MODEL>::dump(
    std::vector<std::shared_ptr<Diagnostic>> const& diagnostics)
{
    (void)diagnostics; // maybe_unused formates weird // TODO REMOVE
    using PATCH     = diagnostic::SamraiPatch<MODEL>;
    auto& model     = this->model_;
    auto& hierarchy = this->hierarchy_;
    auto& hifile    = this->hi5_.file_;

    constexpr diagnostic::Level LEVEL = L;
    if constexpr (LEVEL == diagnostic::Level::DBG) {}

    hifile.createGroup("/timestep");
    size_t ts_idx = 0;
    [](HighFive::File& file, size_t idx) {
        const auto ts_idx_str = std::to_string(idx);
        std::stringstream ss;
        ss << "/timestep/" << ts_idx_str;
        auto ts_group = file.createGroup(ss.str());
        ts_group.createAttribute<std::string>("index", HighFive::DataSpace::From(ts_idx_str))
            .write(ts_idx_str);
        ss << "/patchlevel/";
        file.createGroup(ss.str());
    }(hifile, ts_idx);

    std::string _0_path("/timestep/" + std::to_string(ts_idx) + "/patchlevel/");
    H5Easy::detail::createGroupsToDataSet(hifile, _0_path);
    for (int iLevel = 0; iLevel < hierarchy.getNumberOfLevels(); ++iLevel)
    {
        std::stringstream levelPath;
        levelPath << _0_path << iLevel;
        hifile.createGroup(levelPath.str());

        auto const& level = hierarchy.getPatchLevel(iLevel);

        size_t patch_idx = 0;
        for (auto& patch : *level)
        {
            std::stringstream patchPath;
            patchPath << levelPath.str() << "/patch/" << patch_idx;
            hifile.createGroup(patchPath.str());

            auto spatch = std::make_unique<PATCH>(*patch, model);

            auto writeFields = [&spatch, &patchPath, &hifile](
                                   auto& f, std::vector<diagnostic::KeyVector>&& kvs, auto eb) {
                for (auto kv : spatch->getBorE(f, kvs))
                {
                    std::stringstream fieldPath;
                    fieldPath << patchPath.str() << "/" << eb << kv.k;
                    HighFive::DataSet dataset = hifile.template createDataSet<double>(
                        fieldPath.str(), HighFive::DataSpace::From(kv.v));
                    dataset.write(kv.v);
                }
            };
            if (1 /*writeB*/)
                writeFields(model.state.electromag.B, {"x", "y", "z"}, "B");
            if (1 /*writeE*/)
                writeFields(model.state.electromag.E, {"x", "y", "z"}, "E");

            if constexpr (LEVEL == diagnostic::Level::DBG)
            {
                // writeGhostNodes
            }

            patch_idx++;
        }
    }
}
} /*namespace PHARE*/

#endif /*PHARE_AMR_DIAGNOSTIC_SAMRAI_HIGHFIVE_H*/
