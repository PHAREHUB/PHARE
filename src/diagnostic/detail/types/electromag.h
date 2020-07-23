#ifndef PHARE_DIAGNOSTIC_DETAIL_TYPES_ELECTROMAG_H
#define PHARE_DIAGNOSTIC_DETAIL_TYPES_ELECTROMAG_H

#include "diagnostic/detail/h5typewriter.h"

namespace PHARE::diagnostic::h5
{
/*
 * Possible outputs
 *
 * /t#/pl#/p#/electromag_(B, E)/(x,y,z)
 */
template<typename HighFiveDiagnostic>
class ElectromagDiagnosticWriter : public H5TypeWriter<HighFiveDiagnostic>
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

    ElectromagDiagnosticWriter(HighFiveDiagnostic& hi5)
        : H5TypeWriter<HighFiveDiagnostic>(hi5)
    {
    }
    void write(DiagnosticProperties&) override;
    void compute(DiagnosticProperties&) override {}

    void createFiles(DiagnosticProperties& diagnostic) override;

    void getDataSetInfo(DiagnosticProperties& diagnostic, std::size_t iLevel,
                        std::string const& patchID, Attributes& patchAttributes) override;

    void initDataSets(DiagnosticProperties& diagnostic,
                      std::unordered_map<std::size_t, std::vector<std::string>> const& patchIDs,
                      Attributes& patchAttributes, std::size_t maxLevel) override;

    void writeAttributes(
        DiagnosticProperties&, Attributes&,
        std::unordered_map<std::size_t, std::vector<std::pair<std::string, Attributes>>>&,
        std::size_t maxLevel) override;

    void finalize(DiagnosticProperties& diagnostic) override;

private:
    std::unordered_map<std::string, std::unique_ptr<HighFiveFile>> fileData;
};


template<typename HighFiveDiagnostic>
void ElectromagDiagnosticWriter<HighFiveDiagnostic>::createFiles(DiagnosticProperties& diagnostic)
{
    for (auto* vecField : this->hi5_.modelView().getElectromagFields())
        checkCreateFileFor_(diagnostic, fileData, "/", vecField->name());
}


template<typename HighFiveDiagnostic>
void ElectromagDiagnosticWriter<HighFiveDiagnostic>::getDataSetInfo(
    DiagnosticProperties& diagnostic, std::size_t iLevel, std::string const& patchID,
    Attributes& patchAttributes)
{
    auto& hi5              = this->hi5_;
    auto vecFields         = hi5.modelView().getElectromagFields();
    std::string lvlPatchID = std::to_string(iLevel) + "_" + patchID;

    auto infoVF = [&](auto& vecF, std::string name, auto& attr) {
        for (auto& [id, type] : core::Components::componentMap)
        {
            attr[name][id]             = vecF.getComponent(type).size();
            attr[name][id + "_ghosts"] = static_cast<std::size_t>(GridLayout::nbrGhosts(
                GridLayout::centering(vecF.getComponent(type).physicalQuantity())[0]));
        }
    };

    for (auto* vecField : vecFields)
    {
        auto& name = vecField->name();
        if (diagnostic.quantity == "/" + name)
            infoVF(*vecField, name, patchAttributes[lvlPatchID]);
    }
}


template<typename HighFiveDiagnostic>
void ElectromagDiagnosticWriter<HighFiveDiagnostic>::initDataSets(
    DiagnosticProperties& diagnostic,
    std::unordered_map<std::size_t, std::vector<std::string>> const& patchIDs,
    Attributes& patchAttributes, std::size_t maxLevel)
{
    auto& hi5      = this->hi5_;
    auto& file     = fileData.at(diagnostic.quantity)->file();
    auto vecFields = hi5.modelView().getElectromagFields();

    auto initVF = [&](auto& path, auto& attr, std::string key, auto null) {
        for (auto& [id, type] : core::Components::componentMap)
        {
            auto vFPath = path + "/" + key + "_" + id;
            hi5.template createDataSet<float>(file, vFPath,
                                              null ? 0 : attr[key][id].template to<std::size_t>());
            this->writeGhostsAttr_(file, vFPath,
                                   null ? 0 : attr[key][id + "_ghosts"].template to<std::size_t>(),
                                   null);
        }
    };

    auto initPatch = [&](auto& level, auto& attr, std::string patchID = "") {
        bool null = patchID.empty();
        std::string path{hi5.getPatchPathAddTimestamp(level, patchID)};
        for (auto* vecField : vecFields)
        {
            auto& name = vecField->name();
            if (diagnostic.quantity == "/" + name)
                initVF(path, attr, name, null);
        }
    };

    initDataSets_(patchIDs, patchAttributes, maxLevel, initPatch);
}



template<typename HighFiveDiagnostic>
void ElectromagDiagnosticWriter<HighFiveDiagnostic>::write(DiagnosticProperties& diagnostic)
{
    auto& hi5 = this->hi5_;

    for (auto* vecField : hi5.modelView().getElectromagFields())
    {
        auto& name = vecField->name();
        if (diagnostic.quantity == "/" + name)
        {
            auto& file = fileData.at(diagnostic.quantity)->file();
            hi5.writeVecFieldAsDataset(file, hi5.patchPath() + "/" + name, *vecField);
        }
    }
}



template<typename HighFiveDiagnostic>
void ElectromagDiagnosticWriter<HighFiveDiagnostic>::writeAttributes(
    DiagnosticProperties& diagnostic, Attributes& fileAttributes,
    std::unordered_map<std::size_t, std::vector<std::pair<std::string, Attributes>>>&
        patchAttributes,
    std::size_t maxLevel)
{
    writeAttributes_(fileData.at(diagnostic.quantity)->file(), fileAttributes, patchAttributes,
                     maxLevel);
}


template<typename HighFiveDiagnostic>
void ElectromagDiagnosticWriter<HighFiveDiagnostic>::finalize(DiagnosticProperties& diagnostic)
{
    fileData.erase(diagnostic.quantity);
    assert(fileData.count(diagnostic.quantity) == 0);
}


} // namespace PHARE::diagnostic::h5

#endif /* PHARE_DIAGNOSTIC_DETAIL_TYPES_ELECTROMAG_H */
