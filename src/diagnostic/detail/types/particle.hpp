#ifndef PHARE_DIAGNOSTIC_DETAIL_TYPES_PARTICLE_HPP
#define PHARE_DIAGNOSTIC_DETAIL_TYPES_PARTICLE_HPP

#include "diagnostic/detail/h5typewriter.hpp"
#include "diagnostic/detail/h5_utils.hpp"

#include "core/data/particles/particle_packer.hpp"

#include "amr/data/particles/particles_data.hpp"

#include <unordered_map>
#include <string>
#include <memory>



namespace PHARE::diagnostic::h5
{
/*
 * It is assumed thateach patch has equal number of populations
 *
 * Possible outputs
 *
 * /t#/pl#/p#/ions/pop_(1,2,...)/domain/(weight, charge, iCell, delta, v)
 * /t#/pl#/p#/ions/pop_(1,2,...)/levelGhost/(weight, charge, iCell, delta, v)
 * /t#/pl#/p#/ions/pop_(1,2,...)/patchGhost/(weight, charge, iCell, delta, v)
 */
template<typename H5Writer>
class ParticlesDiagnosticWriter : public H5TypeWriter<H5Writer>
{
public:
    using Super = H5TypeWriter<H5Writer>;
    using Super::checkCreateFileFor_;
    using Super::fileData_;
    using Super::h5Writer_;
    using Super::initDataSets_;
    using Super::writeAttributes_;
    using Super::writeGhostsAttr_;
    using Super::writeIonPopAttributes_;
    static constexpr auto dimension   = H5Writer::dimension;
    static constexpr auto interpOrder = H5Writer::interpOrder;
    using Attributes                  = typename Super::Attributes;
    using Packer                      = core::ParticlePacker<dimension>;
    using FloatType                   = typename H5Writer::FloatType;

    ParticlesDiagnosticWriter(H5Writer& h5Writer)
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

private:
    // PGI compiler (nvc++ 21.3-0) doesn't like static initializations of arrays
    std::array<std::string, 5> packer_keys_ = core::packer_keys();
};


template<typename H5Writer>
void ParticlesDiagnosticWriter<H5Writer>::createFiles(DiagnosticProperties& diagnostic)
{
    for (auto const& pop : this->h5Writer_.modelView().getIons())
    {
        std::string tree{"/ions/pop/" + pop.name() + "/"};
        checkCreateFileFor_(diagnostic, fileData_, tree, "domain", "levelGhost", "patchGhost");
    }
}

template<typename H5Writer>
void ParticlesDiagnosticWriter<H5Writer>::getDataSetInfo(DiagnosticProperties& diagnostic,
                                                         std::size_t iLevel,
                                                         std::string const& patchID,
                                                         Attributes& patchAttributes)
{
    auto& h5Writer         = this->h5Writer_;
    std::string lvlPatchID = std::to_string(iLevel) + "_" + patchID;

    auto getSize = [&](auto const& value, auto n_particles) {
        using ValueType = std::decay_t<decltype(value)>;
        if (n_particles == 0)
            return std::vector<std::size_t>{0};
        if constexpr (is_array_dataset<ValueType, dimension>)
            return std::vector<std::size_t>{n_particles, value.size()};
        else /* not an array so value one of type ValueType*/
            return std::vector<std::size_t>{n_particles, 1};
    };

    auto particleInfo = [&](auto& attr, auto& particles) {
        std::size_t part_idx = 0;
        auto const& keys     = packer_keys_;
        core::apply(Packer::empty(), [&](auto const& arg) {
            attr[keys[part_idx]] = getSize(arg, particles.size());
            ++part_idx;
        });
    };

    auto checkInfo = [&](auto& tree, auto pType, auto& attr, auto& ps) {
        std::string active{tree + pType};
        if (diagnostic.quantity == active)
            particleInfo(attr[pType], ps);
    };

    for (auto& pop : h5Writer.modelView().getIons())
    {
        std::string tree{"/ions/pop/" + pop.name() + "/"};
        auto& popAttr = patchAttributes[lvlPatchID][pop.name()];
        checkInfo(tree, "domain", popAttr, pop.domainParticles());
        checkInfo(tree, "levelGhost", popAttr, pop.levelGhostParticles());
        checkInfo(tree, "patchGhost", popAttr, pop.patchGhostParticles());
    }
}


template<typename H5Writer>
void ParticlesDiagnosticWriter<H5Writer>::initDataSets(
    DiagnosticProperties& diagnostic,
    std::unordered_map<std::size_t, std::vector<std::string>> const& patchIDs,
    Attributes& patchAttributes, std::size_t maxLevel)
{
    auto& h5Writer = this->h5Writer_;
    auto& h5file   = fileData_.at(diagnostic.quantity)->file();

    auto createDataSet = [&](auto&& path, auto& attr, auto& key, auto& value, auto null) {
        using ValueType = std::decay_t<decltype(value)>;

        auto shape = null ? std::vector<std::size_t>(0)
                          : attr[key].template to<std::vector<std::size_t>>();

        if constexpr (is_array_dataset<ValueType, dimension>)
        {
            return h5Writer.template createDataSet<typename ValueType::value_type>(h5file, path,
                                                                                   shape);
        }
        else
            return h5Writer.template createDataSet<ValueType>(h5file, path, shape);
    };

    auto initDataSet = [&](auto& lvl, auto& patchID, auto& attr) {
        bool null = patchID.empty();
        std::string path{h5Writer_.getPatchPathAddTimestamp(lvl, patchID) + "/"};
        std::size_t part_idx = 0;
        core::apply(Packer::empty(), [&](auto const& arg) {
            auto const& keys = packer_keys_;
            createDataSet(path + keys[part_idx], attr, keys[part_idx], arg, null);
            ++part_idx;
        });
        this->writeGhostsAttr_(h5file, path, core::ghostWidthForParticles<interpOrder>(), null);
    };

    auto initIfActive = [&](auto& lvl, auto& tree, auto& attr, auto& pop, auto& patch, auto var) {
        if (diagnostic.quantity == tree + var)
            initDataSet(lvl, patch, patch.empty() ? attr : attr[pop][var]);
    };

    auto initPatch = [&](auto& lvl, auto& attr, std::string patchID = "") {
        for (auto& pop : h5Writer.modelView().getIons())
        {
            std::string tree{"/ions/pop/" + pop.name() + "/"};
            initIfActive(lvl, tree, attr, pop.name(), patchID, "domain");
            initIfActive(lvl, tree, attr, pop.name(), patchID, "levelGhost");
            initIfActive(lvl, tree, attr, pop.name(), patchID, "patchGhost");
        }
    };

    initDataSets_(patchIDs, patchAttributes, maxLevel, initPatch);
}


template<typename H5Writer>
void ParticlesDiagnosticWriter<H5Writer>::write(DiagnosticProperties& diagnostic)
{
    auto& h5Writer = this->h5Writer_;

    auto writeParticles = [&](auto path, auto& particles) {
        if (particles.size() == 0)
            return;

        auto& h5file = *fileData_.at(diagnostic.quantity);
        Packer packer(particles);
        core::ContiguousParticles<dimension> copy{particles.size()};
        packer.pack(copy);

        auto const& keys = packer_keys_;
        h5file.template write_data_set_flat<2>(path + keys[0], copy.weight.data());
        h5file.template write_data_set_flat<2>(path + keys[1], copy.charge.data());
        h5file.template write_data_set_flat<2>(path + keys[2], copy.iCell.data());
        h5file.template write_data_set_flat<2>(path + keys[3], copy.delta.data());
        h5file.template write_data_set_flat<2>(path + keys[4], copy.v.data());
    };

    auto checkWrite = [&](auto& tree, auto pType, auto& ps) {
        std::string active{tree + pType};
        if (diagnostic.quantity == active)
            writeParticles(h5Writer.patchPath() + "/", ps);
    };

    for (auto& pop : h5Writer.modelView().getIons())
    {
        std::string tree{"/ions/pop/" + pop.name() + "/"};
        checkWrite(tree, "domain", pop.domainParticles());
        checkWrite(tree, "levelGhost", pop.levelGhostParticles());
        checkWrite(tree, "patchGhost", pop.patchGhostParticles());
    }
}


template<typename H5Writer>
void ParticlesDiagnosticWriter<H5Writer>::writeAttributes(
    DiagnosticProperties& diagnostic, Attributes& fileAttributes,
    std::unordered_map<std::size_t, std::vector<std::pair<std::string, Attributes>>>&
        patchAttributes,
    std::size_t maxLevel)
{
    auto& h5Writer = this->h5Writer_;
    auto& h5file   = fileData_.at(diagnostic.quantity)->file();

    auto checkWrite = [&](auto& tree, std::string pType, auto const& pop) {
        if (diagnostic.quantity == tree + pType)
            this->writeIonPopAttributes_(h5file, pop);
    };

    for (auto& pop : h5Writer.modelView().getIons())
    {
        std::string tree = "/ions/pop/" + pop.name() + "/";
        checkWrite(tree, "domain", pop);
        checkWrite(tree, "levelGhost", pop);
        checkWrite(tree, "patchGhost", pop);
    }

    writeAttributes_(diagnostic, h5file, fileAttributes, patchAttributes, maxLevel);
}

} // namespace PHARE::diagnostic::h5


#endif /* PHARE_DIAGNOSTIC_DETAIL_TYPES_PARTICLE_H */
