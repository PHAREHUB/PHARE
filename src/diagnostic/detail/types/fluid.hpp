#ifndef PHARE_DIAGNOSTIC_DETAIL_TYPES_FLUID_HPP
#define PHARE_DIAGNOSTIC_DETAIL_TYPES_FLUID_HPP

#include "diagnostic/detail/h5typewriter.hpp"
#include "core/numerics/interpolator/interpolator.hpp"

#include "core/data/vecfield/vecfield_component.hpp"

namespace PHARE::diagnostic::h5
{
/*
 * It is assumed that each patch has equal number of populations
 *
 * Possible outputs
 *
 * /t#/pl#/p#/ions/density
 * /t#/pl#/p#/ions/bulkVelocity/(x,y,z)
 * /t#/pl#/p#/ions/pop_(1,2,...)/density
 * /t#/pl#/p#/ions/pop_(1,2,...)/bulkVelocity/(x,y,z)
 */
template<typename H5Writer>
class FluidDiagnosticWriter : public H5TypeWriter<H5Writer>
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
    using Attributes = typename Super::Attributes;
    using GridLayout = typename H5Writer::GridLayout;
    using FloatType  = typename H5Writer::FloatType;

    static constexpr auto dimension    = GridLayout::dimension;
    static constexpr auto interp_order = GridLayout::interp_order;


    FluidDiagnosticWriter(H5Writer& h5Writer)
        : Super{h5Writer}
    {
    }
    void write(DiagnosticProperties&) override;
    void compute(DiagnosticProperties&) override;

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
    auto isActiveDiag(DiagnosticProperties const& diagnostic, std::string const& tree,
                      std::string var)
    {
        return diagnostic.quantity == tree + var;
    };
};




template<typename H5Writer>
void FluidDiagnosticWriter<H5Writer>::compute(DiagnosticProperties& diagnostic)
{
    core::MomentumTensorInterpolator<dimension, interp_order> interpolator;

    auto& h5Writer  = this->h5Writer_;
    auto& modelView = h5Writer.modelView();
    auto& ions      = modelView.getIons();
    auto minLvl     = this->h5Writer_.minLevel;
    auto maxLvl     = this->h5Writer_.maxLevel;
    // compute the momentum tensor for each population that requires it
    // compute for all ions but that requires the computation of all pop
    std::string tree{"/ions/"};
    if (isActiveDiag(diagnostic, tree, "momentum_tensor"))
    {
        auto computeMomentumTensor
            = [&](GridLayout& layout, std::string patchID, std::size_t iLvel) {
                  for (auto& pop : ions)
                  {
                      std::string tree{"/ions/pop/" + pop.name() + "/"};
                      auto& pop_momentum_tensor = pop.momentumTensor();
                      pop_momentum_tensor.zero();
                      auto domainParts     = core::makeIndexRange(pop.domainParticles());
                      auto patchGhostParts = core::makeIndexRange(pop.patchGhostParticles());

                      // dumps occur after the last substep but before the next first substep
                      // at this time, levelGhostPartsNew is emptied and not yet filled
                      // and the former levelGhostPartsNew has been moved to levelGhostPartsOld
                      auto levelGhostParts = core::makeIndexRange(pop.levelGhostParticlesOld());

                      interpolator(domainParts, pop_momentum_tensor, layout, pop.mass());
                      interpolator(patchGhostParts, pop_momentum_tensor, layout, pop.mass());
                      interpolator(levelGhostParts, pop_momentum_tensor, layout, pop.mass());
                  }
                  ions.computeFullMomentumTensor();
              };
        modelView.visitHierarchy(computeMomentumTensor, minLvl, maxLvl);
    }
    else // if not computing total momentum tensor, user may want to compute it for some pop
    {
        for (auto& pop : ions)
        {
            std::string tree{"/ions/pop/" + pop.name() + "/"};

            auto computePopMomentumTensor
                = [&](GridLayout& layout, std::string patchID, std::size_t iLvel) {
                      auto& pop_momentum_tensor = pop.momentumTensor();
                      pop_momentum_tensor.zero();
                      auto domainParts     = core::makeIndexRange(pop.domainParticles());
                      auto patchGhostParts = core::makeIndexRange(pop.patchGhostParticles());

                      // dumps occur after the last substep but before the next first substep
                      // at this time, levelGhostPartsNew is emptied and not yet filled
                      // and the former levelGhostPartsNew has been moved to levelGhostPartsOld
                      auto levelGhostParts = core::makeIndexRange(pop.levelGhostParticlesOld());

                      interpolator(domainParts, pop_momentum_tensor, layout, pop.mass());
                      interpolator(patchGhostParts, pop_momentum_tensor, layout, pop.mass());
                      interpolator(levelGhostParts, pop_momentum_tensor, layout, pop.mass());
                  };
            if (isActiveDiag(diagnostic, tree, "momentum_tensor"))
            {
                modelView.visitHierarchy(computePopMomentumTensor, minLvl, maxLvl);
            }
        }
    }
}



template<typename H5Writer>
void FluidDiagnosticWriter<H5Writer>::createFiles(DiagnosticProperties& diagnostic)
{
    for (auto const& pop : this->h5Writer_.modelView().getIons())
    {
        std::string tree{"/ions/pop/" + pop.name() + "/"};
        checkCreateFileFor_(diagnostic, fileData_, tree, "density", "flux", "momentum_tensor");
    }

    std::string tree{"/ions/"};
    checkCreateFileFor_(diagnostic, fileData_, tree, "density", "mass_density", "bulkVelocity",
                        "momentum_tensor");
}




template<typename H5Writer>
void FluidDiagnosticWriter<H5Writer>::getDataSetInfo(DiagnosticProperties& diagnostic,
                                                     std::size_t iLevel, std::string const& patchID,
                                                     Attributes& patchAttributes)
{
    auto& h5Writer = this->h5Writer_;
    auto& ions     = h5Writer.modelView().getIons();
    std::string lvlPatchID{std::to_string(iLevel) + "_" + patchID};


    auto setGhostNbr = [](auto const& field, auto& attr, auto const& name) {
        auto ghosts              = GridLayout::nDNbrGhosts(field.physicalQuantity());
        attr[name + "_ghosts_x"] = static_cast<std::size_t>(ghosts[0]);
        if constexpr (GridLayout::dimension > 1)
            attr[name + "_ghosts_y"] = static_cast<std::size_t>(ghosts[1]);
        if constexpr (GridLayout::dimension > 2)
            attr[name + "_ghosts_z"] = static_cast<std::size_t>(ghosts[2]);
    };

    auto infoDS = [&](auto& field, std::string name, auto& attr) {
        // highfive doesn't accept uint32 which ndarray.shape() is
        auto const& shape = field.shape();
        attr[name]        = std::vector<std::size_t>(shape.data(), shape.data() + shape.size());
        setGhostNbr(field, attr, name);
    };

    auto infoVF = [&](auto& vecF, std::string name, auto& attr) {
        for (auto const& [id, type] : core::VectorComponents::map())
            infoDS(vecF.getComponent(type), name + "_" + id, attr);
    };

    auto infoTF = [&](auto& tensorF, std::string name, auto& attr) {
        for (auto const& [id, type] : core::TensorComponents::map())
            infoDS(tensorF.getComponent(type), name + "_" + id, attr);
    };

    for (auto& pop : ions)
    {
        std::string tree{"/ions/pop/" + pop.name() + "/"};
        auto& popAttr = patchAttributes[lvlPatchID]["fluid_" + pop.name()];
        if (isActiveDiag(diagnostic, tree, "density"))
            infoDS(pop.density(), "density", popAttr);
        if (isActiveDiag(diagnostic, tree, "flux"))
            infoVF(pop.flux(), "flux", popAttr);
        if (isActiveDiag(diagnostic, tree, "momentum_tensor"))
            infoTF(pop.momentumTensor(), "momentum_tensor", popAttr);
    }

    std::string tree{"/ions/"};
    if (isActiveDiag(diagnostic, tree, "density"))
        infoDS(ions.density(), "density", patchAttributes[lvlPatchID]["ion"]);
    if (isActiveDiag(diagnostic, tree, "mass_density"))
        infoDS(ions.massDensity(), "mass_density", patchAttributes[lvlPatchID]["ion"]);
    if (isActiveDiag(diagnostic, tree, "bulkVelocity"))
        infoVF(ions.velocity(), "bulkVelocity", patchAttributes[lvlPatchID]["ion"]);
    if (isActiveDiag(diagnostic, tree, "momentum_tensor"))
        infoTF(ions.momentumTensor(), "momentum_tensor", patchAttributes[lvlPatchID]["ion"]);
}




template<typename H5Writer>
void FluidDiagnosticWriter<H5Writer>::initDataSets(
    DiagnosticProperties& diagnostic,
    std::unordered_map<std::size_t, std::vector<std::string>> const& patchIDs,
    Attributes& patchAttributes, std::size_t maxLevel)
{
    auto& h5Writer = this->h5Writer_;
    auto& ions     = h5Writer.modelView().getIons();
    auto& h5file   = *fileData_.at(diagnostic.quantity);

    auto writeGhosts = [&](auto& path, auto& attr, std::string key, auto null) {
        this->writeGhostsAttr_(h5file, path,
                               null ? 0 : attr[key + "_ghosts_x"].template to<std::size_t>(), null);
        if constexpr (GridLayout::dimension > 1)
            this->writeGhostsAttr_(
                h5file, path, null ? 0 : attr[key + "_ghosts_y"].template to<std::size_t>(), null);
        if constexpr (GridLayout::dimension > 2)
            this->writeGhostsAttr_(
                h5file, path, null ? 0 : attr[key + "_ghosts_z"].template to<std::size_t>(), null);
    };

    auto initDS = [&](auto& path, auto& attr, std::string key, auto null) {
        auto dsPath = path + key;
        h5Writer.template createDataSet<FloatType>(
            h5file, dsPath,
            null ? std::vector<std::size_t>(GridLayout::dimension, 0)
                 : attr[key].template to<std::vector<std::size_t>>());
        writeGhosts(dsPath, attr, key, null);
    };
    auto initVF = [&](auto& path, auto& attr, std::string key, auto null) {
        for (auto& [id, type] : core::Components::componentMap())
            initDS(path, attr, key + "_" + id, null);
    };
    auto initTF = [&](auto& path, auto& attr, std::string key, auto null) {
        for (auto& [id, type] : core::Components::componentMap<2>())
            initDS(path, attr, key + "_" + id, null);
    };

    auto initPatch = [&](auto& lvl, auto& attr, std::string patchID = "") {
        bool null        = patchID.empty();
        std::string path = h5Writer.getPatchPathAddTimestamp(lvl, patchID) + "/";

        for (auto& pop : ions)
        {
            std::string popId{"fluid_" + pop.name()};
            std::string tree{"/ions/pop/" + pop.name() + "/"};
            std::string popPath(path + "pop/" + pop.name() + "/");
            if (isActiveDiag(diagnostic, tree, "density"))
                initDS(path, attr[popId], "density", null);
            if (isActiveDiag(diagnostic, tree, "flux"))
                initVF(path, attr[popId], "flux", null);
            if (isActiveDiag(diagnostic, tree, "momentum_tensor"))
                initTF(path, attr[popId], "momentum_tensor", null);
        }

        std::string tree{"/ions/"};
        if (isActiveDiag(diagnostic, tree, "density"))
            initDS(path, attr["ion"], "density", null);
        if (isActiveDiag(diagnostic, tree, "mass_density"))
            initDS(path, attr["ion"], "mass_density", null);
        if (isActiveDiag(diagnostic, tree, "bulkVelocity"))
            initVF(path, attr["ion"], "bulkVelocity", null);
        if (isActiveDiag(diagnostic, tree, "momentum_tensor"))
            initTF(path, attr["ion"], "momentum_tensor", null);
    };

    initDataSets_(patchIDs, patchAttributes, maxLevel, initPatch);
}


template<typename H5Writer>
void FluidDiagnosticWriter<H5Writer>::write(DiagnosticProperties& diagnostic)
{
    auto& h5Writer = this->h5Writer_;
    auto& ions     = h5Writer.modelView().getIons();
    auto& h5file   = *fileData_.at(diagnostic.quantity);

    auto writeDS = [&](auto path, auto& field) {
        h5file.template write_data_set_flat<GridLayout::dimension>(path, &(*field.begin()));
    };
    auto writeTF
        = [&](auto path, auto& vecF) { h5Writer.writeTensorFieldAsDataset(h5file, path, vecF); };

    std::string path = h5Writer.patchPath() + "/";
    for (auto& pop : ions)
    {
        std::string tree{"/ions/pop/" + pop.name() + "/"};
        if (isActiveDiag(diagnostic, tree, "density"))
            writeDS(path + "density", pop.density());
        if (isActiveDiag(diagnostic, tree, "flux"))
            writeTF(path + "flux", pop.flux());
        if (isActiveDiag(diagnostic, tree, "momentum_tensor"))
            writeTF(path + "momentum_tensor", pop.momentumTensor());
    }

    std::string tree{"/ions/"};
    if (isActiveDiag(diagnostic, tree, "density"))
        writeDS(path + "density", ions.density());
    if (isActiveDiag(diagnostic, tree, "mass_density"))
        writeDS(path + "mass_density", ions.massDensity());
    if (isActiveDiag(diagnostic, tree, "bulkVelocity"))
        writeTF(path + "bulkVelocity", ions.velocity());
    if (isActiveDiag(diagnostic, tree, "momentum_tensor"))
        writeTF(path + "momentum_tensor", ions.momentumTensor());
}


template<typename H5Writer>
void FluidDiagnosticWriter<H5Writer>::writeAttributes(
    DiagnosticProperties& diagnostic, Attributes& fileAttributes,
    std::unordered_map<std::size_t, std::vector<std::pair<std::string, Attributes>>>&
        patchAttributes,
    std::size_t maxLevel)
{
    auto& h5Writer = this->h5Writer_;
    auto& h5file   = *fileData_.at(diagnostic.quantity);

    auto checkWrite = [&](auto& tree, std::string qty, auto const& pop) {
        if (diagnostic.quantity == tree + qty)
            this->writeIonPopAttributes_(h5file, pop);
    };

    for (auto& pop : h5Writer.modelView().getIons())
    {
        std::string tree = "/ions/pop/" + pop.name() + "/";
        checkWrite(tree, "density", pop);
        checkWrite(tree, "flux", pop);
        checkWrite(tree, "momentum_tensor", pop);
    }

    writeAttributes_(diagnostic, h5file, fileAttributes, patchAttributes, maxLevel);
}

} // namespace PHARE::diagnostic::h5

#endif /* PHARE_DIAGNOSTIC_DETAIL_TYPES_FLUID_H */
