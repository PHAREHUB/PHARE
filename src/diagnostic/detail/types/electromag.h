#ifndef PHARE_DIAGNOSTIC_DETAIL_TYPES_ELECTROMAG_H
#define PHARE_DIAGNOSTIC_DETAIL_TYPES_ELECTROMAG_H

#include "diagnostic/detail/highfive.h"

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
    using Attributes = typename Super::Attributes;
    using GridLayout = typename HighFiveDiagnostic::GridLayout;

    ElectromagDiagnosticWriter(HighFiveDiagnostic& hi5)
        : H5TypeWriter<HighFiveDiagnostic>(hi5)
    {
    }
    void write(DiagnosticProperties&) override;
    void compute(DiagnosticProperties&) override {}
    void getDataSetInfo(DiagnosticProperties& diagnostic, size_t iLevel, std::string const& patchID,
                        Attributes& patchAttributes) override;
    void initDataSets(DiagnosticProperties& diagnostic,
                      std::unordered_map<size_t, std::vector<std::string>> const& patchIDs,
                      Attributes& patchAttributes, size_t maxLevel) override;
    void
    writeAttributes(DiagnosticProperties&, Attributes&,
                    std::unordered_map<size_t, std::vector<std::pair<std::string, Attributes>>>&,
                    size_t maxLevel) override;

private:
    std::unordered_map<std::string, std::unique_ptr<HighFiveFile>> fileData;
};


template<typename HighFiveDiagnostic>
void ElectromagDiagnosticWriter<HighFiveDiagnostic>::getDataSetInfo(
    DiagnosticProperties& diagnostic, size_t iLevel, std::string const& patchID,
    Attributes& patchAttributes)
{
    auto& hi5              = this->hi5_;
    auto vecFields         = hi5.modelView().getElectromagFields();
    std::string lvlPatchID = std::to_string(iLevel) + "_" + patchID;

    auto infoVF = [&](auto& vecF, std::string name, auto& attr) {
        for (auto& [id, type] : core::Components::componentMap)
        {
            attr[name][id]             = vecF.getComponent(type).size();
            attr[name][id + "_ghosts"] = static_cast<size_t>(GridLayout::nbrGhosts(
                GridLayout::centering(vecF.getComponent(type).physicalQuantity())[0]));
        }
    };

    for (auto* vecField : vecFields)
    {
        auto& name = vecField->name();
        if (diagnostic.quantity == "/" + name)
        {
            infoVF(*vecField, name, patchAttributes[lvlPatchID]);
            if (!fileData.count(diagnostic.quantity))
                fileData.emplace(diagnostic.quantity, hi5.makeFile(diagnostic));
        }
    }
}


template<typename HighFiveDiagnostic>
void ElectromagDiagnosticWriter<HighFiveDiagnostic>::initDataSets(
    DiagnosticProperties& diagnostic,
    std::unordered_map<size_t, std::vector<std::string>> const& patchIDs,
    Attributes& patchAttributes, size_t maxLevel)
{
    auto& hi5      = this->hi5_;
    auto& file     = fileData.at(diagnostic.quantity)->file();
    auto vecFields = hi5.modelView().getElectromagFields();

    auto initVF = [&](auto& path, auto& attr, std::string key, auto null) {
        for (auto& [id, type] : core::Components::componentMap)
        {
            auto vFPath = path + "/" + key + "/" + id;
            hi5.template createDataSet<float>(file, vFPath,
                                              null ? 0 : attr[key][id].template to<size_t>());
            this->writeGhostsAttr_(
                file, vFPath, null ? 0 : attr[key][id + "_ghosts"].template to<size_t>(), null);
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
    std::unordered_map<size_t, std::vector<std::pair<std::string, Attributes>>>& patchAttributes,
    size_t maxLevel)
{
    writeAttributes_(fileData.at(diagnostic.quantity)->file(), diagnostic, fileAttributes,
                     patchAttributes, maxLevel);
}

} // namespace PHARE::diagnostic::h5

#endif /* PHARE_DIAGNOSTIC_DETAIL_TYPES_ELECTROMAG_H */
