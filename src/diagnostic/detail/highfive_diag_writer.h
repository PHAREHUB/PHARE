#ifndef HIGHFIVEDIAGNOSTICWRITER_H
#define HIGHFIVEDIAGNOSTICWRITER_H

#include <string>
#include <unordered_map>

#include "diagnostic/diagnostic_writer.h"

namespace PHARE::diagnostic::h5
{
template<typename HighFiveDiagnostic>
class Hi5DiagnosticWriter : public DiagnosticWriter
{
public:
    using Attributes = typename HighFiveDiagnostic::Attributes;
    Hi5DiagnosticWriter(HighFiveDiagnostic& hi5)
        : hi5_(hi5)
    {
    }

    virtual void getDataSetInfo(DiagnosticDAO& diagnostic, size_t iLevel,
                                std::string const& patchID, Attributes& patchAttributes)
        = 0;
    virtual void initDataSets(DiagnosticDAO& diagnostic,
                              std::unordered_map<size_t, std::vector<std::string>> const& patchIDs,
                              Attributes& patchAttributes, int maxLevel)
        = 0;

protected:
    template<typename InitPatch>
    void initDataSets_(std::unordered_map<size_t, std::vector<std::string>> const& patchIDs,
                       Attributes& patchAttributes, std::size_t levelNbr, InitPatch&& initPatch)
    {
        for (size_t lvl = 0; lvl <= levelNbr; lvl++)
        {
            auto& lvlPatches  = patchIDs.at(lvl);
            size_t patchNbr   = lvlPatches.size();
            size_t maxPatches = hi5_.getMaxOfPerMPI(patchNbr);
            for (size_t i = 0; i < patchNbr; i++)
            {
                std::string lvlPatchID = std::to_string(lvl) + "_" + lvlPatches[i];
                initPatch(lvl, patchAttributes[lvlPatchID], lvlPatches[i]);
            }
            for (size_t i = patchNbr; i < maxPatches; i++)
                initPatch(lvl, patchAttributes);
        }
    }



    HighFiveDiagnostic& hi5_;
};

} // namespace PHARE::diagnostic::h5

#endif // HIGHFIVEDIAGNOSTICWRITER_H
