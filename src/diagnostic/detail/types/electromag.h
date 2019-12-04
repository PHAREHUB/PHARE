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
class ElectromagDiagnosticWriter : public Hi5DiagnosticWriter<HighFiveDiagnostic>
{
public:
    using Hi5DiagnosticWriter<HighFiveDiagnostic>::hi5_;
    using Attributes = typename Hi5DiagnosticWriter<HighFiveDiagnostic>::Attributes;
    ElectromagDiagnosticWriter(HighFiveDiagnostic& hi5)
        : Hi5DiagnosticWriter<HighFiveDiagnostic>(hi5)
    {
    }
    void write(DiagnosticDAO&) override;
    void compute(DiagnosticDAO&) override {}
    void getDataSetInfo(DiagnosticDAO& diagnostic, size_t iLevel, std::string const& patchID,
                        Attributes& patchAttributes) override;
    void initDataSets(DiagnosticDAO& diagnostic,
                      std::unordered_map<size_t, std::vector<std::string>> const& patchIDs,
                      Attributes& patchAttributes, int maxLevel) override;
};


template<typename HighFiveDiagnostic>
void ElectromagDiagnosticWriter<HighFiveDiagnostic>::getDataSetInfo(DiagnosticDAO& diagnostic,
                                                                    size_t iLevel,
                                                                    std::string const& patchID,
                                                                    Attributes& patchAttributes)
{
    auto& hi5      = this->hi5_;
    auto vecFields = hi5.modelView().getElectromagFields();
    std::string lvlPatchID{std::to_string(iLevel) + "_" + patchID};

    for (auto* vecField : vecFields)
    {
        auto& name = vecField->name();
        if (diagnostic.subtype == "/" + name)
            for (auto& [id, type] : core::Components::componentMap)
                patchAttributes[lvlPatchID][name][id] = vecField->getComponent(type).size();
    }
}


template<typename HighFiveDiagnostic>
void ElectromagDiagnosticWriter<HighFiveDiagnostic>::initDataSets(
    DiagnosticDAO& diagnostic, std::unordered_map<size_t, std::vector<std::string>> const& patchIDs,
    Attributes& patchAttributes, int maxLevel)
{
    auto& hi5      = this->hi5_;
    auto vecFields = hi5.modelView().getElectromagFields();

    auto initPatch = [&](auto& level, auto& attributes, std::string patchID = "") {
        bool null = patchID.empty();
        std::string path{hi5.getPatchPath("time", level, patchID)};
        for (auto* vecField : vecFields)
        {
            auto& name = vecField->name();
            if (diagnostic.subtype == "/" + name)
                for (auto& [id, type] : core::Components::componentMap)
                    hi5.template createDataSet<float>(
                        path + "/" + name + "/" + id,
                        null ? 0 : attributes[name][id].template to<size_t>());
        }
    };

    Hi5DiagnosticWriter<HighFiveDiagnostic>::initDataSets_(patchIDs, patchAttributes, maxLevel,
                                                           initPatch);
}



template<typename HighFiveDiagnostic>
void ElectromagDiagnosticWriter<HighFiveDiagnostic>::write([
    [maybe_unused]] DiagnosticDAO& diagnostic)
{
    auto& hi5 = this->hi5_;

    for (auto* vecField : hi5.modelView().getElectromagFields())
    {
        auto& name = vecField->name();
        if (diagnostic.subtype == "/" + name)
            hi5.writeVecFieldAsDataset(hi5.patchPath() + "/" + name, *vecField);
    }
}
} // namespace PHARE::diagnostic::h5

#endif /* PHARE_DIAGNOSTIC_DETAIL_TYPES_ELECTROMAG_H */
