#ifndef HIGHFIVEDIAGNOSTICWRITER_H
#define HIGHFIVEDIAGNOSTICWRITER_H

#include <string>
#include <algorithm>
#include <unordered_map>

#include "diagnostic/diagnostic_writer.h"

#include "core/utilities/mpi_utils.h"

#include "highfive/H5File.hpp"

namespace PHARE::diagnostic::h5
{
template<typename Writer>
class H5TypeWriter : public PHARE::diagnostic::TypeWriter
{
public:
    using Attributes = typename Writer::Attributes;
    H5TypeWriter(Writer& hi5)
        : hi5_(hi5)
    {
    }

    virtual void createFiles(DiagnosticProperties& diagnostic) = 0;

    virtual void getDataSetInfo(DiagnosticProperties& diagnostic, std::size_t iLevel,
                                std::string const& patchID, Attributes& patchAttributes)
        = 0;

    virtual void
    initDataSets(DiagnosticProperties& diagnostic,
                 std::unordered_map<std::size_t, std::vector<std::string>> const& patchIDs,
                 Attributes& patchAttributes, std::size_t maxLevel)
        = 0;

    virtual void writeAttributes(
        DiagnosticProperties&, Attributes&,
        std::unordered_map<std::size_t, std::vector<std::pair<std::string, Attributes>>>&,
        std::size_t maxLevel)
        = 0;

    virtual void finalize(DiagnosticProperties& diagnostic) = 0;


protected:
    template<typename InitPatch>
    void initDataSets_(std::unordered_map<std::size_t, std::vector<std::string>> const& patchIDs,
                       Attributes& patchAttributes, std::size_t maxLevel, InitPatch&& initPatch)
    {
        for (size_t lvl = hi5_.minLevel; lvl <= maxLevel; lvl++)
        {
            auto& lvlPatches       = patchIDs.at(lvl);
            std::size_t patchNbr   = lvlPatches.size();
            std::size_t maxPatches = core::mpi::max(patchNbr);
            for (size_t i = 0; i < patchNbr; i++)
                initPatch(lvl, patchAttributes[std::to_string(lvl) + "_" + lvlPatches[i]],
                          lvlPatches[i]);
            for (size_t i = patchNbr; i < maxPatches; i++)
                initPatch(lvl, patchAttributes);
        }
    }

    void writeAttributes_(
        HighFive::File& file, Attributes& fileAttributes,
        std::unordered_map<std::size_t, std::vector<std::pair<std::string, Attributes>>>&
            patchAttributes,
        std::size_t maxLevel)
    {
        for (size_t lvl = hi5_.minLevel; lvl <= maxLevel; lvl++)
        {
            auto& lvlPatches       = patchAttributes.at(lvl);
            std::size_t patchNbr   = lvlPatches.size();
            std::size_t maxPatches = core::mpi::max(patchNbr);
            for (auto const& [patch, attr] : lvlPatches)
                hi5_.writeAttributeDict(file, attr, hi5_.getPatchPathAddTimestamp(lvl, patch));
            for (size_t i = patchNbr; i < maxPatches; i++)
                hi5_.writeAttributeDict(file, hi5_.modelView().getEmptyPatchProperties(), "");
        }

        hi5_.writeAttributeDict(file, fileAttributes, "/");
    }

    void writeGhostsAttr_(HighFive::File& file, std::string path, std::size_t ghosts, bool null)
    {
        Attributes dsAttr;
        dsAttr["ghosts"] = ghosts;
        hi5_.writeAttributeDict(file, dsAttr, null ? "" : path);
    }

    template<typename FileMap, typename... Quantities>
    void checkCreateFileFor_(DiagnosticProperties const& diagnostic, FileMap& fileData,
                             std::string const tree, Quantities const... vars)
    {
        core::apply(std::forward_as_tuple(vars...), [&](auto const& var) {
            if (diagnostic.quantity == tree + var and !fileData.count(diagnostic.quantity))
                fileData.emplace(diagnostic.quantity, this->hi5_.makeFile(diagnostic));
        });
    }


    Writer& hi5_;
};

} // namespace PHARE::diagnostic::h5

#endif // HIGHFIVEDIAGNOSTICWRITER_H
