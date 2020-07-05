#ifndef PHARE_DIAGNOSTIC_DETAIL_TYPES_FLUID_H
#define PHARE_DIAGNOSTIC_DETAIL_TYPES_FLUID_H

#include "diagnostic/detail/h5typewriter.h"

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
class FluidDiagnosticWriter : public H5TypeWriter<HighFiveDiagnostic>
{
public:
    using Super = H5TypeWriter<HighFiveDiagnostic>;
    using Super::hi5_;
    using Super::initDataSets_;
    using Super::writeAttributes_;
    using Super::writeGhostsAttr_;
    using Super::checkCreateFileFor_;
    using Attributes = typename Super::Attributes;
    using GridLayout = typename HighFiveDiagnostic::GridLayout;

    FluidDiagnosticWriter(HighFiveDiagnostic& hi5)
        : H5TypeWriter<HighFiveDiagnostic>(hi5)
    {
    }
    void write(DiagnosticProperties&) override;
    void compute(DiagnosticProperties&) override {}

    void createFiles(DiagnosticProperties& diagnostic) override;

    void getDataSetInfo(DiagnosticProperties& diagnostic, size_t iLevel, std::string const& patchID,
                        Attributes& patchAttributes) override;

    void initDataSets(DiagnosticProperties& diagnostic,
                      std::unordered_map<size_t, std::vector<std::string>> const& patchIDs,
                      Attributes& patchAttributes, size_t maxLevel) override;

    void
    writeAttributes(DiagnosticProperties&, Attributes&,
                    std::unordered_map<size_t, std::vector<std::pair<std::string, Attributes>>>&,
                    size_t maxLevel) override;

    void finalize(DiagnosticProperties& diagnostic) override;

private:
    std::unordered_map<std::string, std::unique_ptr<HighFiveFile>> fileData;
};



template<typename HighFiveDiagnostic>
void FluidDiagnosticWriter<HighFiveDiagnostic>::createFiles(DiagnosticProperties& diagnostic)
{
    for (auto const& pop : this->hi5_.modelView().getIons())
    {
        std::string tree{"/ions/pop/" + pop.name() + "/"};
        checkCreateFileFor_(diagnostic, fileData, tree, "density", "flux");
    }

    std::string tree{"/ions/"};
    checkCreateFileFor_(diagnostic, fileData, tree, "density", "bulkVelocity");
}

template<typename HighFiveDiagnostic>
void FluidDiagnosticWriter<HighFiveDiagnostic>::getDataSetInfo(DiagnosticProperties& diagnostic,
                                                               size_t iLevel,
                                                               std::string const& patchID,
                                                               Attributes& patchAttributes)
{
    auto& hi5  = this->hi5_;
    auto& ions = hi5.modelView().getIons();
    std::string lvlPatchID{std::to_string(iLevel) + "_" + patchID};

    auto checkActive = [&](auto& tree, auto var) { return diagnostic.quantity == tree + var; };

    auto infoDS = [&](auto& field, std::string name, auto& attr) {
        attr[name]             = field.size();
        attr[name + "_ghosts"] = static_cast<size_t>(
            GridLayout::nbrGhosts(GridLayout::centering(field.physicalQuantity())[0]));
    };

    auto infoVF = [&](auto& vecF, std::string name, auto& attr) {
        for (auto& [id, type] : core::Components::componentMap)
        {
            attr[name][id]             = vecF.getComponent(type).size();
            attr[name][id + "_ghosts"] = static_cast<size_t>(GridLayout::nbrGhosts(
                GridLayout::centering(vecF.getComponent(type).physicalQuantity())[0]));
        }
    };

    for (auto& pop : ions)
    {
        std::string tree{"/ions/pop/" + pop.name() + "/"};
        auto& popAttr = patchAttributes[lvlPatchID]["fluid_" + pop.name()];
        if (checkActive(tree, "density"))
            infoDS(pop.density(), "density", popAttr);
        if (checkActive(tree, "flux"))
            infoVF(pop.flux(), "flux", popAttr);
    }

    std::string tree{"/ions/"};
    if (checkActive(tree, "density"))
        infoDS(ions.density(), "density", patchAttributes[lvlPatchID]["ion"]);
    if (checkActive(tree, "bulkVelocity"))
        infoVF(ions.velocity(), "bulkVelocity", patchAttributes[lvlPatchID]["ion"]);
}


template<typename HighFiveDiagnostic>
void FluidDiagnosticWriter<HighFiveDiagnostic>::initDataSets(
    DiagnosticProperties& diagnostic,
    std::unordered_map<size_t, std::vector<std::string>> const& patchIDs,
    Attributes& patchAttributes, size_t maxLevel)
{
    auto& hi5  = this->hi5_;
    auto& ions = hi5.modelView().getIons();
    auto& file = fileData.at(diagnostic.quantity)->file();

    auto checkActive = [&](auto& tree, auto var) { return diagnostic.quantity == tree + var; };

    auto initDS = [&](auto& path, auto& attr, std::string key, auto null) {
        auto dsPath = path + key;
        hi5.template createDataSet<float>(file, dsPath, null ? 0 : attr[key].template to<size_t>());
        this->writeGhostsAttr_(file, dsPath, null ? 0 : attr[key + "_ghosts"].template to<size_t>(),
                               null);
    };
    auto initVF = [&](auto& path, auto& attr, std::string key, auto null) {
        for (auto& [id, type] : core::Components::componentMap)
        {
            auto vFPath = path + key + "_" + id;
            hi5.template createDataSet<float>(file, vFPath,
                                              null ? 0 : attr[key][id].template to<size_t>());
            this->writeGhostsAttr_(
                file, vFPath, null ? 0 : attr[key][id + "_ghosts"].template to<size_t>(), null);
        }
    };

    auto initPatch = [&](auto& lvl, auto& attr, std::string patchID = "") {
        bool null        = patchID.empty();
        std::string path = hi5.getPatchPathAddTimestamp(lvl, patchID) + "/";

        for (auto& pop : ions)
        {
            std::string popId{"fluid_" + pop.name()};
            std::string tree{"/ions/pop/" + pop.name() + "/"};
            std::string popPath(path + "pop/" + pop.name() + "/");
            if (checkActive(tree, "density"))
                initDS(path, attr[popId], "density", null);
            if (checkActive(tree, "flux"))
                initVF(path, attr[popId], "flux", null);
        }

        std::string tree{"/ions/"};
        if (checkActive(tree, "density"))
            initDS(path, attr["ion"], "density", null);
        if (checkActive(tree, "bulkVelocity"))
            initVF(path, attr["ion"], "bulkVelocity", null);
    };

    initDataSets_(patchIDs, patchAttributes, maxLevel, initPatch);
}


template<typename HighFiveDiagnostic>
void FluidDiagnosticWriter<HighFiveDiagnostic>::write(DiagnosticProperties& diagnostic)
{
    auto& hi5  = this->hi5_;
    auto& ions = hi5.modelView().getIons();
    auto& file = fileData.at(diagnostic.quantity)->file();

    auto checkActive = [&](auto& tree, auto var) { return diagnostic.quantity == tree + var; };
    auto writeDS     = [&](auto path, auto& field) { hi5.writeDataSet(file, path, field.data()); };
    auto writeVF     = [&](auto path, auto& vecF) { hi5.writeVecFieldAsDataset(file, path, vecF); };

    std::string path = hi5.patchPath() + "/";
    for (auto& pop : ions)
    {
        std::string tree{"/ions/pop/" + pop.name() + "/"};
        if (checkActive(tree, "density"))
            writeDS(path + "density", pop.density());
        if (checkActive(tree, "flux"))
            writeVF(path + "flux", pop.flux());
    }

    std::string tree{"/ions/"};
    auto& density = ions.density();
    if (checkActive(tree, "density"))
        writeDS(path + "density", density);
    if (checkActive(tree, "bulkVelocity"))
        writeVF(path + "bulkVelocity", ions.velocity());
}


template<typename HighFiveDiagnostic>
void FluidDiagnosticWriter<HighFiveDiagnostic>::writeAttributes(
    DiagnosticProperties& diagnostic, Attributes& fileAttributes,
    std::unordered_map<size_t, std::vector<std::pair<std::string, Attributes>>>& patchAttributes,
    size_t maxLevel)
{
    writeAttributes_(fileData.at(diagnostic.quantity)->file(), diagnostic, fileAttributes,
                     patchAttributes, maxLevel);
}

template<typename HighFiveDiagnostic>
void FluidDiagnosticWriter<HighFiveDiagnostic>::finalize(DiagnosticProperties& diagnostic)
{
    fileData.erase(diagnostic.quantity);
    assert(fileData.count(diagnostic.quantity) == 0);
}


} // namespace PHARE::diagnostic::h5

#endif /* PHARE_DIAGNOSTIC_DETAIL_TYPES_FLUID_H */
