#ifndef PHARE_DIAGNOSTIC_DETAIL_TYPES_INFO_HPP
#define PHARE_DIAGNOSTIC_DETAIL_TYPES_INFO_HPP

#include "diagnostic/detail/h5typewriter.hpp"


namespace PHARE::diagnostic::h5
{


template<typename H5Writer>
class InfoDiagnosticWriter : public H5TypeWriter<H5Writer>
{
public:
    using Super = H5TypeWriter<H5Writer>;
    using Super::checkCreateFileFor_;
    using Super::fileData_;
    using Super::h5Writer_;
    using Super::initDataSets_;
    using Super::writeAttributes_;
    using Attributes                = Super::Attributes;
    using ModelViewVariant          = Super::ModelViewVariant;
    using FloatType                 = H5Writer::FloatType;
    static constexpr auto dimension = H5Writer::dimension;

    InfoDiagnosticWriter(H5Writer& h5Writer)
        : Super{h5Writer}
    {
    }

    void write(DiagnosticProperties&) override;

    void compute_as_needed(DiagnosticProperties&) override {}

    void createFiles(DiagnosticProperties& diagnostic) override;

    void getDataSetInfo(DiagnosticProperties& diagnostic, std::size_t iLevel,
                        std::string const& patchID, Attributes& patchAttributes,
                        ModelViewVariant& modelView) override;

    void initDataSets(DiagnosticProperties& diagnostic,
                      std::unordered_map<std::size_t, std::vector<std::string>> const& patchIDs,
                      Attributes& patchAttributes, std::size_t maxLevel) override;

    void writeAttributes(
        DiagnosticProperties&, Attributes&,
        std::unordered_map<std::size_t, std::vector<std::pair<std::string, Attributes>>>&,
        std::size_t maxLevel) override;
};


template<typename H5Writer>
void InfoDiagnosticWriter<H5Writer>::createFiles(DiagnosticProperties& diagnostic)
{
    std::string tree{"/"};
    checkCreateFileFor_(diagnostic, fileData_, tree, "particle_count");
}


template<typename H5Writer>
void InfoDiagnosticWriter<H5Writer>::getDataSetInfo(DiagnosticProperties& diagnostic,
                                                    std::size_t iLevel, std::string const& patchID,
                                                    Attributes& patchAttributes,
                                                    ModelViewVariant& /*modelView*/)
{
}


template<typename H5Writer>
void InfoDiagnosticWriter<H5Writer>::initDataSets(
    DiagnosticProperties& diagnostic,
    std::unordered_map<std::size_t, std::vector<std::string>> const& patchIDs,
    Attributes& patchAttributes, std::size_t maxLevel)
{
}



template<typename H5Writer>
void InfoDiagnosticWriter<H5Writer>::write(DiagnosticProperties& diagnostic)
{
}


template<typename H5Writer>
void InfoDiagnosticWriter<H5Writer>::writeAttributes(
    DiagnosticProperties& diagnostic, Attributes& fileAttributes,
    std::unordered_map<std::size_t, std::vector<std::pair<std::string, Attributes>>>&
        patchAttributes,
    std::size_t maxLevel)
{
    auto& h5Writer = this->h5Writer_;

    // every patch on every level needs the attribute (writes are collective per key), so visit
    //  through the mapper: count on hybrid levels, zero on levels owned by models without particles
    std::size_t lvl_idx = -1, p_idx = 0;
    auto gatherParticleCounts = [&](auto& /*gridLayout*/, std::string const& patchID,
                                    std::size_t iLevel, auto& modelViewVariant) {
        if (iLevel != lvl_idx)
        {
            lvl_idx = iLevel;
            p_idx   = 0;
        }

        auto& patches = patchAttributes[iLevel];
        assert(patches[p_idx].first == patchID);
        patches[p_idx].second["particle_count"] = std::visit(
            [](auto& modelView) -> std::size_t {
                using Model_t = std::decay_t<decltype(modelView)>::Model_t;
                if constexpr (solver::is_hybrid_model_v<Model_t>)
                    return sum_from(modelView.getIons(),
                                    [](auto const& pop) { return pop.domainParticles().size(); });
                else
                    return 0;
            },
            modelViewVariant);
        ++p_idx;
    };

    Attributes defaultPatchAttributes;

    if (diagnostic.quantity == "/particle_count")
    {
        h5Writer.mapper().visitHierarchy(gatherParticleCounts, h5Writer_.minLevel, maxLevel);
        defaultPatchAttributes["particle_count"] = std::size_t{0};
    }

    writeAttributes_(diagnostic, Super::h5FileForQuantity(diagnostic), fileAttributes,
                     patchAttributes, maxLevel, defaultPatchAttributes);
}


} // namespace PHARE::diagnostic::h5

#endif /* PHARE_DIAGNOSTIC_DETAIL_TYPES_INFO_H */
