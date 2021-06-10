#ifndef HIGHFIVEDIAGNOSTICWRITER_H
#define HIGHFIVEDIAGNOSTICWRITER_H

#include <string>
#include <algorithm>
#include <unordered_map>

#include "diagnostic/diagnostic_writer.h"
#include "diagnostic/detail/h5file.h"

#include "core/utilities/mpi_utils.h"

#include "highfive/H5File.hpp"

namespace PHARE::diagnostic::h5
{
template<typename Writer>
class H5TypeWriter : public PHARE::diagnostic::TypeWriter
{
public:
    using Attributes                = typename Writer::Attributes;
    static constexpr auto dimension = Writer::dimension;

    H5TypeWriter(Writer& h5Writer)
        : h5Writer_{h5Writer}
    {
    }

    //------  defined by each concrete H5TypeWriter---------------------------
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
    //------------------------------------------------------------------------



    //------ valid for all h5type writers -------------------------------------
    void finalize(DiagnosticProperties& diagnostic)
    {
        // we close the file by removing the associated file
        // from the map. This is done only at flush time otherwise
        ++diagnostic.dumpIdx;

        assert(diagnostic.params.contains("flush_every"));

        std::size_t flushEvery = diagnostic.param<std::size_t>("flush_every");

        if (flushEvery != Writer::flush_never and diagnostic.dumpIdx % flushEvery == 0)
        {
            fileData_.erase(diagnostic.quantity);
            assert(fileData_.count(diagnostic.quantity) == 0);
        }
    }
    //------------------------------------------------------------------------


protected:
    template<typename InitPatch>
    void initDataSets_(std::unordered_map<std::size_t, std::vector<std::string>> const& patchIDs,
                       Attributes& patchAttributes, std::size_t maxLevel, InitPatch&& initPatch)
    {
        for (std::size_t lvl = h5Writer_.minLevel; lvl <= maxLevel; lvl++)
        {
            auto& lvlPatches       = patchIDs.at(lvl);
            std::size_t patchNbr   = lvlPatches.size();
            std::size_t maxPatches = core::mpi::max(patchNbr);
            for (std::size_t i = 0; i < patchNbr; ++i)
                initPatch(lvl, patchAttributes[std::to_string(lvl) + "_" + lvlPatches[i]],
                          lvlPatches[i]);
            for (std::size_t i = patchNbr; i < maxPatches; i++)
                initPatch(lvl, patchAttributes);
        }
    }

    void writeAttributes_(
        DiagnosticProperties& diagnostic, HighFive::File& file, Attributes& fileAttributes,
        std::unordered_map<std::size_t, std::vector<std::pair<std::string, Attributes>>>&
            patchAttributes,
        std::size_t maxLevel)
    {
        for (std::size_t lvl = h5Writer_.minLevel; lvl <= maxLevel; lvl++)
        {
            auto& lvlPatches       = patchAttributes.at(lvl);
            std::size_t patchNbr   = lvlPatches.size();
            std::size_t maxPatches = core::mpi::max(patchNbr);
            for (auto const& [patch, attr] : lvlPatches)
                h5Writer_.writeAttributeDict(file, attr,
                                             h5Writer_.getPatchPathAddTimestamp(lvl, patch));
            for (std::size_t i = patchNbr; i < maxPatches; i++)
                h5Writer_.writeAttributeDict(file, h5Writer_.modelView().getEmptyPatchProperties(),
                                             "");
        }

        if (diagnostic.nAttributes > 0)
            h5Writer_.writeAttributeDict(file, diagnostic.fileAttributes, "/py_attrs");
        h5Writer_.writeAttributeDict(file, fileAttributes, "/");
    }

    template<typename ParticlePopulation>
    void writeIonPopAttributes_(HighFive::File& file, ParticlePopulation const& pop)
    {
        auto& h5Writer = this->h5Writer_;

        Attributes popAttributes;
        popAttributes["pop_mass"] = pop.mass();
        h5Writer.writeAttributeDict(file, popAttributes, "/");
    }

    void writeGhostsAttr_(HighFive::File& file, std::string path, std::size_t ghosts, bool null)
    {
        Attributes dsAttr;
        dsAttr["ghosts"] = ghosts;
        h5Writer_.writeAttributeDict(file, dsAttr, null ? "" : path);
    }

    template<typename FileMap, typename... Quantities>
    void checkCreateFileFor_(DiagnosticProperties const& diagnostic, FileMap& fileData,
                             std::string const tree, Quantities const... vars)
    {
        core::apply(std::forward_as_tuple(vars...), [&](auto const& var) {
            if (diagnostic.quantity == tree + var and !fileData.count(diagnostic.quantity))
                fileData.emplace(diagnostic.quantity, this->h5Writer_.makeFile(diagnostic));
        });
    }


    Writer& h5Writer_;
    std::unordered_map<std::string, std::unique_ptr<HighFiveFile>> fileData_;
};

} // namespace PHARE::diagnostic::h5

#endif // HIGHFIVEDIAGNOSTICWRITER_H
