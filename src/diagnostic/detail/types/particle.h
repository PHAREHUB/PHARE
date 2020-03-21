#ifndef PHARE_DIAGNOSTIC_DETAIL_TYPES_PARTICLE_H
#define PHARE_DIAGNOSTIC_DETAIL_TYPES_PARTICLE_H

#include "diagnostic/detail/highfive.h"
#include "diagnostic/detail/h5_utils.h"

#include "core/data/particles/particle_packer.h"

#include "amr/data/particles/particles_data.h"

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
template<typename HighFiveDiagnostic>
class ParticlesDiagnosticWriter : public H5TypeWriter<HighFiveDiagnostic>
{
public:
    using Super = H5TypeWriter<HighFiveDiagnostic>;
    using Super::hi5_;
    using Super::initDataSets_;
    using Super::writeAttributes_;
    using Super::writeGhostsAttr_;
    static constexpr auto dimension   = HighFiveDiagnostic::dimension;
    static constexpr auto interpOrder = HighFiveDiagnostic::interpOrder;
    using Attributes                  = typename Super::Attributes;
    using Packer                      = core::ParticlePacker<dimension>;

    ParticlesDiagnosticWriter(HighFiveDiagnostic& hi5)
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
    size_t levels_ = 0;
    std::unordered_map<std::string, std::unique_ptr<HighFiveFile>> fileData;
};


template<typename HighFiveDiagnostic>
void ParticlesDiagnosticWriter<HighFiveDiagnostic>::getDataSetInfo(DiagnosticProperties& diagnostic,
                                                                   size_t iLevel,
                                                                   std::string const& patchID,
                                                                   Attributes& patchAttributes)
{
    auto& hi5              = this->hi5_;
    std::string lvlPatchID = std::to_string(iLevel) + "_" + patchID;

    auto getSize = [&](auto const& value) {
        using ValueType = std::decay_t<decltype(value)>;
        if constexpr (is_array_dataset<ValueType, dimension>)
            return value.size();
        else
            return 1; /* not an array so value one of type ValueType*/
    };

    auto particleInfo = [&](auto& attr, auto& particles) {
        size_t part_idx = 0;
        core::apply(Packer::empty(), [&](auto const& arg) {
            attr[Packer::keys()[part_idx]] = getSize(arg) * particles.size();
            part_idx++;
        });
    };

    auto checkInfo = [&](auto& tree, auto pType, auto& attr, auto& ps) {
        std::string active{tree + pType};
        if (diagnostic.quantity == active)
        {
            particleInfo(attr[pType], ps);
            if (!fileData.count(diagnostic.quantity))
                fileData.emplace(diagnostic.quantity, hi5.makeFile(diagnostic));
        }
    };

    for (auto& pop : hi5.modelView().getIons())
    {
        std::string tree{"/ions/pop/" + pop.name() + "/"};
        auto& popAttr = patchAttributes[lvlPatchID][pop.name()];
        checkInfo(tree, "domain", popAttr, pop.domainParticles());
        checkInfo(tree, "levelGhost", popAttr, pop.levelGhostParticles());
        checkInfo(tree, "patchGhost", popAttr, pop.patchGhostParticles());
    }
    levels_ = iLevel;
}


template<typename HighFiveDiagnostic>
void ParticlesDiagnosticWriter<HighFiveDiagnostic>::initDataSets(
    DiagnosticProperties& diagnostic,
    std::unordered_map<size_t, std::vector<std::string>> const& patchIDs,
    Attributes& patchAttributes, size_t maxLevel)
{
    auto& hi5  = this->hi5_;
    auto& file = fileData.at(diagnostic.quantity)->file();

    auto createDataSet = [&](auto&& path, auto size, auto const& value) {
        using ValueType = std::decay_t<decltype(value)>;
        if constexpr (is_array_dataset<ValueType, dimension>)
            return hi5.template createDataSet<typename ValueType::value_type>(file, path, size);
        else
            return hi5.template createDataSet<ValueType>(file, path, size);
    };

    auto initDataSet = [&](auto& lvl, auto& patchID, auto tree, auto& attr) {
        bool null = patchID.empty();
        std::string path{hi5_.getPatchPathAddTimestamp(lvl, patchID) + tree + "/"};
        size_t part_idx = 0;
        core::apply(Packer::empty(), [&](auto const& arg) {
            createDataSet(path + Packer::keys()[part_idx],
                          null ? 0 : attr[Packer::keys()[part_idx]].template to<size_t>(), arg);
            part_idx++;
        });
        this->writeGhostsAttr_(file, path, amr::ghostWidthForParticles<interpOrder>(), null);
    };

    auto initIfActive = [&](auto& lvl, auto& tree, auto& attr, auto& pop, auto& patch, auto var) {
        if (diagnostic.quantity == tree + var)
            initDataSet(lvl, patch, tree + var, patch.empty() ? attr : attr[pop][var]);
    };

    auto initPatch = [&](auto& lvl, auto& attr, std::string patchID = "") {
        for (auto& pop : hi5.modelView().getIons())
        {
            std::string tree{"/ions/pop/" + pop.name() + "/"};
            initIfActive(lvl, tree, attr, pop.name(), patchID, "domain");
            initIfActive(lvl, tree, attr, pop.name(), patchID, "levelGhost");
            initIfActive(lvl, tree, attr, pop.name(), patchID, "patchGhost");
        }
    };

    initDataSets_(patchIDs, patchAttributes, maxLevel, initPatch);
}


template<typename HighFiveDiagnostic>
void ParticlesDiagnosticWriter<HighFiveDiagnostic>::write(DiagnosticProperties& diagnostic)
{
    auto& hi5 = this->hi5_;

    auto writeParticles = [&](auto path, auto& particles) {
        if (particles.size() == 0)
            return;
        auto& hfile = fileData.at(diagnostic.quantity)->file();
        Packer packer(particles);
        core::ContiguousParticles<dimension> copy{particles.size()};
        packer.pack(copy);

        hi5.writeDataSet(hfile, path + packer.keys()[0], copy.weight.data());
        hi5.writeDataSet(hfile, path + packer.keys()[1], copy.charge.data());
        hi5.writeDataSet(hfile, path + packer.keys()[2], copy.iCell.data());
        hi5.writeDataSet(hfile, path + packer.keys()[3], copy.delta.data());
        hi5.writeDataSet(hfile, path + packer.keys()[4], copy.v.data());
    };

    auto checkWrite = [&](auto& tree, auto pType, auto& ps) {
        std::string active{tree + pType};
        if (diagnostic.quantity == active)
            writeParticles(hi5.patchPath() + active + "/", ps);
    };

    for (auto& pop : hi5.modelView().getIons())
    {
        std::string tree{"/ions/pop/" + pop.name() + "/"};
        checkWrite(tree, "domain", pop.domainParticles());
        checkWrite(tree, "levelGhost", pop.levelGhostParticles());
        checkWrite(tree, "patchGhost", pop.patchGhostParticles());
    }
}


template<typename HighFiveDiagnostic>
void ParticlesDiagnosticWriter<HighFiveDiagnostic>::writeAttributes(
    DiagnosticProperties& diagnostic, Attributes& fileAttributes,
    std::unordered_map<size_t, std::vector<std::pair<std::string, Attributes>>>& patchAttributes,
    size_t maxLevel)
{
    writeAttributes_(fileData.at(diagnostic.quantity)->file(), diagnostic, fileAttributes,
                     patchAttributes, maxLevel);
}

} // namespace PHARE::diagnostic::h5

#endif /* PHARE_DIAGNOSTIC_DETAIL_TYPES_PARTICLE_H */
