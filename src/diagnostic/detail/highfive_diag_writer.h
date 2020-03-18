#ifndef HIGHFIVEDIAGNOSTICWRITER_H
#define HIGHFIVEDIAGNOSTICWRITER_H

#include <string>
#include <algorithm>
#include <unordered_map>

#include "diagnostic/diagnostic_writer.h"

namespace PHARE::diagnostic::h5
{
template<typename HighFiveDiagnostic>
class Hi5DiagnosticTypeWriter : public PHARE::diagnostic::DiagnosticTypeWriter
{
public:
    using Attributes = typename HighFiveDiagnostic::Attributes;
    Hi5DiagnosticTypeWriter(HighFiveDiagnostic& hi5)
        : hi5_(hi5)
    {
    }

    virtual void getDataSetInfo(DiagnosticDAO& diagnostic, size_t iLevel,
                                std::string const& patchID, Attributes& patchAttributes)
        = 0;
    virtual void initDataSets(DiagnosticDAO& diagnostic,
                              std::unordered_map<size_t, std::vector<std::string>> const& patchIDs,
                              Attributes& patchAttributes, size_t maxLevel)
        = 0;

    virtual void
    writeAttributes(DiagnosticDAO&, Attributes&,
                    std::unordered_map<size_t, std::vector<std::pair<std::string, Attributes>>>&,
                    size_t maxLevel)
        = 0;

protected:
    template<typename InitPatch>
    void initDataSets_(std::unordered_map<size_t, std::vector<std::string>> const& patchIDs,
                       Attributes& patchAttributes, size_t maxLevel, InitPatch&& initPatch)
    {
        for (size_t lvl = hi5_.minLevel; lvl <= maxLevel; lvl++)
        {
            auto& lvlPatches  = patchIDs.at(lvl);
            size_t patchNbr   = lvlPatches.size();
            size_t maxPatches = hi5_.mpiGetMaxOf(patchNbr);
            for (size_t i = 0; i < patchNbr; i++)
                initPatch(lvl, patchAttributes[std::to_string(lvl) + "_" + lvlPatches[i]],
                          lvlPatches[i]);
            for (size_t i = patchNbr; i < maxPatches; i++)
                initPatch(lvl, patchAttributes);
        }
    }

    void
    writeAttributes_(HighFive::File& file, DiagnosticDAO& diagnostic, Attributes& fileAttributes,
                     std::unordered_map<size_t, std::vector<std::pair<std::string, Attributes>>>&
                         patchAttributes,
                     std::size_t maxLevel)
    {
        for (size_t lvl = hi5_.minLevel; lvl <= maxLevel; lvl++)
        {
            auto& lvlPatches  = patchAttributes.at(lvl);
            size_t patchNbr   = lvlPatches.size();
            size_t maxPatches = hi5_.mpiGetMaxOf(patchNbr);
            for (auto const& [patch, attr] : lvlPatches)
                hi5_.writeAttributeDict(file, attr, hi5_.getPatchPathAddTimestamp(lvl, patch));
            for (size_t i = patchNbr; i < maxPatches; i++)
                hi5_.writeAttributeDict(file, hi5_.modelView().getEmptyPatchAttributes(), "");
        }

        hi5_.writeAttributeDict(file, fileAttributes, "/");
    }

    void writeGhostsAttr_(HighFive::File& file, std::string path, size_t ghosts, bool null)
    {
        Attributes dsAttr;
        dsAttr["ghosts"] = ghosts;
        hi5_.writeAttributeDict(file, dsAttr, null ? "" : path);
    }

    HighFiveDiagnostic& hi5_;
};

} // namespace PHARE::diagnostic::h5

#endif // HIGHFIVEDIAGNOSTICWRITER_H
