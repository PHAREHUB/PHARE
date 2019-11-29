#ifndef PHARE_DIAGNOSTIC_DETAIL_TYPES_PARTICLE_H
#define PHARE_DIAGNOSTIC_DETAIL_TYPES_PARTICLE_H

#include "detail/highfive.h"
#include "detail/h5_utils.h"

namespace PHARE::diagnostic::h5
{
/*
 * It is assumed that MPI processes will have at least one patch, and each patch has equal
 *  number of populations
 * It is assumed  that all processes have at least one patch on all levels
 *
 * Possible outputs
 *
 * /#t/pl#/p#/ions/pop_(1,2,...)/domain/(weight, charge, iCell, delta, v)
 * /#t/pl#/p#/ions/pop_(1,2,...)/levelGhost/(weight, charge, iCell, delta, v)
 * /#t/pl#/p#/ions/pop_(1,2,...)/patchGhost/(weight, charge, iCell, delta, v)
 */
template<typename HighFiveDiagnostic>
class ParticlesDiagnosticWriter : public Hi5DiagnosticWriter<HighFiveDiagnostic>
{
public:
    using Hi5DiagnosticWriter<HighFiveDiagnostic>::hi5_;
    using Attributes                = typename Hi5DiagnosticWriter<HighFiveDiagnostic>::Attributes;
    static constexpr auto dimension = HighFiveDiagnostic::dimension;
    using Packer                    = ParticlePacker<dimension>;

    ParticlesDiagnosticWriter(HighFiveDiagnostic& hi5)
        : Hi5DiagnosticWriter<HighFiveDiagnostic>(hi5)
    {
    }
    void write(DiagnosticDAO&) override;
    void compute(DiagnosticDAO&) override {}
    void getDataSetInfo(DiagnosticDAO& diagnostic, size_t iLevel, std::string const& patchID,
                        Attributes& patchAttributes) override;
    void initDataSets(DiagnosticDAO& diagnostic,
                      std::unordered_map<size_t, std::vector<std::string>> const& patchIDs,
                      Attributes& patchAttributes) override;

private:
    size_t levels_ = 0;
};


template<typename HighFiveDiagnostic>
void ParticlesDiagnosticWriter<HighFiveDiagnostic>::getDataSetInfo(DiagnosticDAO& diagnostic,
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
        if (diagnostic.subtype == active)
            particleInfo(attr[pType], ps);
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
    DiagnosticDAO& diagnostic, std::unordered_map<size_t, std::vector<std::string>> const& patchIDs,
    Attributes& patchAttributes)
{
    auto& hi5 = this->hi5_;

    auto createDataSet = [&hi5](auto&& path, auto size, auto const& value) {
        using ValueType = std::decay_t<decltype(value)>;
        if constexpr (is_array_dataset<ValueType, dimension>)
            return hi5.template createDataSet<typename ValueType::value_type>(path, size);
        else
            return hi5.template createDataSet<ValueType>(path, size);
    };

    auto initDataSet = [&](auto& lvl, auto& patchID, auto tree, auto& attr) {
        std::string path{hi5_.getPatchPath("time", lvl, patchID) + tree + "/"};
        size_t part_idx = 0;
        core::apply(Packer::empty(), [&](auto const& arg) {
            createDataSet(
                path + Packer::keys()[part_idx],
                patchID.empty() ? 0 : attr[Packer::keys()[part_idx]].template to<size_t>(), arg);
            part_idx++;
        });
    };

    auto initIfActive = [&](auto& lvl, auto& tree, auto& attr, auto& pop, auto& patch, auto var) {
        if (diagnostic.subtype == tree + var)
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

    Hi5DiagnosticWriter<HighFiveDiagnostic>::initDataSets_(patchIDs, patchAttributes, levels_,
                                                           initPatch);
}


template<typename HighFiveDiagnostic>
void ParticlesDiagnosticWriter<HighFiveDiagnostic>::write([
    [maybe_unused]] DiagnosticDAO& diagnostic)
{
    auto& hi5 = this->hi5_;

    auto writeDataSet = [&hi5](auto path, auto& start, auto const& value) {
        using ValueType = std::decay_t<decltype(value)>;
        if constexpr (is_array_dataset<ValueType, dimension>)
            hi5.writeDataSetPart(path, start * value.size(), value.size(), value);
        else
            hi5.writeDataSetPart(path, start, 1, value); /*not array, write 1 value*/
    };

    auto writeParticles = [&](auto path, auto& particles) {
        if (!particles.size())
            return;

        Packer packer(particles);

        size_t part_idx = 0, idx = 0;

        while (packer.hasNext())
        {
            part_idx = 0;
            core::apply(packer.next(), [&](auto const& arg) {
                writeDataSet(path + packer.keys()[part_idx], idx, arg);
                part_idx++;
            });
            idx++;
        }
    };

    auto checkWrite = [&](auto& tree, auto pType, auto& ps) {
        std::string active{tree + pType};
        if (diagnostic.subtype == active)
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

} // namespace PHARE::diagnostic::h5

#endif /* PHARE_DIAGNOSTIC_DETAIL_TYPES_PARTICLE_H */
