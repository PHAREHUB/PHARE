#ifndef PHARE_DIAGNOSTIC_DETAIL_TYPES_FLUID_H
#define PHARE_DIAGNOSTIC_DETAIL_TYPES_FLUID_H

#include "diagnostic/detail/highfive.h"

namespace PHARE::diagnostic::h5
{
/*
 * It is assumed that each patch has equal number of populations
 *
 * Possible outputs
 *
 * /t#/pl#/p#/ions/density
 * /t#/pl#/p#/ions/bulkVelocity/(x,y,z)
 * /t#/pl#/p#/ions/pop_(1,2,...)/density
 * /t#/pl#/p#/ions/pop_(1,2,...)/bulkVelocity/(x,y,z)
 */
template<typename HighFiveDiagnostic>
class FluidDiagnosticWriter : public Hi5DiagnosticTypeWriter<HighFiveDiagnostic>
{
public:
    using Hi5DiagnosticTypeWriter<HighFiveDiagnostic>::hi5_;
    using Attributes = typename Hi5DiagnosticTypeWriter<HighFiveDiagnostic>::Attributes;
    FluidDiagnosticWriter(HighFiveDiagnostic& hi5)
        : Hi5DiagnosticTypeWriter<HighFiveDiagnostic>(hi5)
    {
    }
    void write(DiagnosticDAO&, Attributes&, Attributes&) override;
    void compute(DiagnosticDAO&) override {}
    void getDataSetInfo(DiagnosticDAO& diagnostic, size_t iLevel, std::string const& patchID,
                        Attributes& patchAttributes) override;
    void initDataSets(DiagnosticDAO& diagnostic,
                      std::unordered_map<size_t, std::vector<std::string>> const& patchIDs,
                      Attributes& patchAttributes, int maxLevel) override;

private:
    std::unordered_map<std::string, std::unique_ptr<HighFiveFile>> fileData;
};


template<typename HighFiveDiagnostic>
void FluidDiagnosticWriter<HighFiveDiagnostic>::getDataSetInfo(DiagnosticDAO& diagnostic,
                                                               size_t iLevel,
                                                               std::string const& patchID,
                                                               Attributes& patchAttributes)
{
    auto& hi5  = this->hi5_;
    auto& ions = hi5.modelView().getIons();
    std::string lvlPatchID{std::to_string(iLevel) + "_" + patchID};

    auto checkActive = [&](auto& tree, auto var) {
        auto active = tree + var;
        bool b      = diagnostic.type == active;
        if (b && !fileData.count(diagnostic.type))
            fileData.emplace(diagnostic.type, hi5.makeFile(diagnostic));
        return b;
    };

    for (auto& pop : ions)
    {
        std::string popId = "fluid_" + pop.name();
        std::string tree{"/ions/pop/" + pop.name() + "/"};
        auto& popAttr = patchAttributes[lvlPatchID][popId];
        if (checkActive(tree, "density"))
            popAttr["density"] = pop.density().size();
        if (checkActive(tree, "flux"))
            for (auto& [id, type] : core::Components::componentMap)
                popAttr["flux"][id] = pop.flux().getComponent(type).size();
    }

    std::string tree{"/ions/"};
    if (checkActive(tree, "density"))
        patchAttributes[lvlPatchID]["ion"]["density"] = ions.density().size();

    if (checkActive(tree, "bulkVelocity"))
        for (auto& [id, type] : core::Components::componentMap)
            patchAttributes[lvlPatchID]["ion"]["bulkVelocity"][id]
                = ions.velocity().getComponent(type).size();
}


template<typename HighFiveDiagnostic>
void FluidDiagnosticWriter<HighFiveDiagnostic>::initDataSets(
    DiagnosticDAO& diagnostic, std::unordered_map<size_t, std::vector<std::string>> const& patchIDs,
    Attributes& patchAttributes, int maxLevel)
{
    auto& hi5  = this->hi5_;
    auto& ions = hi5.modelView().getIons();

    auto checkActive = [&](auto& tree, auto var) { return diagnostic.type == tree + var; };

    auto initPatch = [&](auto& lvl, auto& attr, std::string patchID = "") {
        bool null = patchID.empty();
        std::string path{hi5.getPatchPath("time", lvl, patchID) + "/ions/"};

        for (auto& pop : ions)
        {
            std::string popId{"fluid_" + pop.name()};
            std::string tree{"/ions/pop/" + pop.name() + "/"};
            std::string popPath(path + "pop/" + pop.name() + "/");
            if (checkActive(tree, "density"))
                hi5.template createDataSet<float>(
                    fileData.at(diagnostic.type)->file(), popPath + "density",
                    null ? 0 : attr[popId]["density"].template to<size_t>());

            if (checkActive(tree, "flux"))
                for (auto& [id, type] : core::Components::componentMap)
                    hi5.template createDataSet<float>(
                        fileData.at(diagnostic.type)->file(), popPath + "flux" + "/" + id,
                        null ? 0 : attr[popId]["flux"][id].template to<size_t>());
        }

        std::string tree{"/ions/"};
        if (checkActive(tree, "density"))
            hi5.template createDataSet<float>(
                fileData.at(diagnostic.type)->file(), path + "density",
                null ? 0 : attr["ion"]["density"].template to<size_t>());

        if (checkActive(tree, "bulkVelocity"))
            for (auto& [id, type] : core::Components::componentMap)
                hi5.template createDataSet<float>(
                    fileData.at(diagnostic.type)->file(), path + "bulkVelocity" + "/" + id,
                    null ? 0 : attr["ion"]["bulkVelocity"][id].template to<size_t>());
    };

    Hi5DiagnosticTypeWriter<HighFiveDiagnostic>::initDataSets_(patchIDs, patchAttributes, maxLevel,
                                                               initPatch);
}


template<typename HighFiveDiagnostic>
void FluidDiagnosticWriter<HighFiveDiagnostic>::write(DiagnosticDAO& diagnostic,
                                                      Attributes& fileAttributes,
                                                      Attributes& patchAttributes)
{
    auto& hi5  = this->hi5_;
    auto& ions = hi5.modelView().getIons();

    auto checkActive = [&](auto& tree, auto var) { return diagnostic.type == tree + var; };

    auto writeAttributes = [&](auto& file) {
        hi5.writeAttributeDict(file, fileAttributes, "/");
        hi5.writeAttributeDict(file, patchAttributes, hi5.patchPath());
    };

    auto writeDS = [&](auto& file, auto path, auto* val) {
        hi5.writeDataSet(file, path, val);
        writeAttributes(file);
    };
    auto writeVF = [&](auto& file, auto path, auto& val) {
        hi5.writeVecFieldAsDataset(file, path, val);
        writeAttributes(file);
    };

    std::string path{hi5.patchPath() + "/"};
    for (auto& pop : ions)
    {
        std::string tree{"/ions/pop/" + pop.name() + "/"};
        std::string popPath{path + tree};
        if (checkActive(tree, "density"))
            writeDS(fileData.at(diagnostic.type)->file(), popPath + "density",
                    pop.density().data());
        if (checkActive(tree, "flux"))
            writeVF(fileData.at(diagnostic.type)->file(), popPath + "flux", pop.flux());
    }

    std::string tree{"/ions/"};
    auto& density = ions.density();
    if (checkActive(tree, "density"))
        writeDS(fileData.at(diagnostic.type)->file(), path + tree + "density", density.data());
    if (checkActive(tree, "bulkVelocity"))
        writeVF(fileData.at(diagnostic.type)->file(), path + tree + "bulkVelocity",
                ions.velocity());
}

} // namespace PHARE::diagnostic::h5

#endif /* PHARE_DIAGNOSTIC_DETAIL_TYPES_FLUID_H */
