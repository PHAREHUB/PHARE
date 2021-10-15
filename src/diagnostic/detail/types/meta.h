#ifndef PHARE_DIAGNOSTIC_DETAIL_TYPES_META_H
#define PHARE_DIAGNOSTIC_DETAIL_TYPES_META_H

#include "diagnostic/detail/h5typewriter.h"


namespace PHARE::diagnostic::h5
{
/*
 * Possible outputs
 *
 * /t#/pl#/p#/tags
 */
template<typename H5Writer>
class MetaDiagnosticWriter : public H5TypeWriter<H5Writer>
{
public:
    using Super = H5TypeWriter<H5Writer>;
    using Super::checkCreateFileFor_;
    using Super::fileData_;
    using Super::h5Writer_;
    using Super::initDataSets_;
    using Super::writeAttributes_;
    using Attributes = typename Super::Attributes;
    using GridLayout = typename H5Writer::GridLayout;
    using FloatType  = typename H5Writer::FloatType;

    MetaDiagnosticWriter(H5Writer& h5Writer)
        : Super{h5Writer}
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




    std::unordered_map<std::string, int*> tags;
};


template<typename H5Writer>
void MetaDiagnosticWriter<H5Writer>::createFiles(DiagnosticProperties& diagnostic)
{
    std::string tree{"/"};
    checkCreateFileFor_(diagnostic, fileData_, tree, "tags");
}


template<typename H5Writer>
void MetaDiagnosticWriter<H5Writer>::getDataSetInfo(DiagnosticProperties& diagnostic,
                                                    std::size_t iLevel, std::string const& patchID,
                                                    Attributes& patchAttributes)
{
    auto& h5Writer         = this->h5Writer_;
    std::string lvlPatchID = std::to_string(iLevel) + "_" + patchID;
    std::string path{h5Writer.getPatchPathAddTimestamp(iLevel, patchID)};

    if (diagnostic.quantity == "/tags" and h5Writer.modelView().hasTagsVectorFor(iLevel, patchID))
    {
        auto& model_tags = *h5Writer.modelView().getTagsVectorFor(iLevel, patchID);
        auto& shape      = model_tags.shape();
        tags[path]       = model_tags.data();
        patchAttributes[lvlPatchID]["tags"]
            = std::vector<std::size_t>(shape.data(), shape.data() + shape.size());
    }
}


template<typename H5Writer>
void MetaDiagnosticWriter<H5Writer>::initDataSets(
    DiagnosticProperties& diagnostic,
    std::unordered_map<std::size_t, std::vector<std::string>> const& patchIDs,
    Attributes& patchAttributes, std::size_t maxLevel)
{
    auto& h5Writer = this->h5Writer_;

    auto initPatch = [&](auto& iLevel, auto& attr, std::string patchID = "") {
        bool null = patchID.empty();

        std::string path{h5Writer.getPatchPathAddTimestamp(iLevel, patchID)};

        if (diagnostic.quantity == "/tags")
            h5Writer.template createDataSet<int>(
                fileData_.at(diagnostic.quantity)->file(), path + "/tags",
                null or tags.count(path) == 0
                    ? std::vector<std::size_t>(GridLayout::dimension, 0)
                    : attr["tags"].template to<std::vector<std::size_t>>());
    };

    initDataSets_(patchIDs, patchAttributes, maxLevel, initPatch);
}



template<typename H5Writer>
void MetaDiagnosticWriter<H5Writer>::write(DiagnosticProperties& diagnostic)
{
    auto& h5Writer = this->h5Writer_;

    if (diagnostic.quantity == "/tags")
    {
        auto& path = h5Writer.patchPath();

        if (tags.count(path) > 0)
        {
            auto& h5 = fileData_.at(diagnostic.quantity);
            h5->template write_data_set_flat<GridLayout::dimension>(path + "/tags", tags[path]);
            tags.erase(path);
        }
    }
}


template<typename H5Writer>
void MetaDiagnosticWriter<H5Writer>::writeAttributes(
    DiagnosticProperties& diagnostic, Attributes& fileAttributes,
    std::unordered_map<std::size_t, std::vector<std::pair<std::string, Attributes>>>&
        patchAttributes,
    std::size_t maxLevel)
{
    writeAttributes_(diagnostic, fileData_.at(diagnostic.quantity)->file(), fileAttributes,
                     patchAttributes, maxLevel);
}


} // namespace PHARE::diagnostic::h5

#endif /* PHARE_DIAGNOSTIC_DETAIL_TYPES_META_H */
