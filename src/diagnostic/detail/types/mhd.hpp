#ifndef PHARE_DIAGNOSTIC_DETAIL_TYPES_MHD_HPP
#define PHARE_DIAGNOSTIC_DETAIL_TYPES_MHD_HPP

#include "core/numerics/primite_conservative_converter/to_primitive_converter.hpp"
#include "diagnostic/detail/h5typewriter.hpp"

#include "core/data/vecfield/vecfield_component.hpp"

namespace PHARE::diagnostic::h5
{
/* Possible outputs
 * /t#/pl#/p#/mhd/density
 * /t#/pl#/p#/mhd/velocity/(x,y,z)
 * /t#/pl#/p#/mhd/B/(x,y,z)
 * /t#/pl#/p#/mhd/pressure
 * /t#/pl#/p#/mhd/rhoV/(x,y,z)
 * /t#/pl#/p#/mhd/Etot
 */
template<typename H5Writer>
class MHDDiagnosticWriter : public H5TypeWriter<H5Writer>
{
public:
    using Super = H5TypeWriter<H5Writer>;
    using Super::checkCreateFileFor_;
    using Super::fileData_;
    using Super::h5Writer_;
    using Super::initDataSets_;
    using Super::writeAttributes_;
    using Super::writeGhostsAttr_;
    using Attributes = typename Super::Attributes;
    using GridLayout = typename H5Writer::GridLayout;
    using FloatType  = typename H5Writer::FloatType;

    static constexpr auto dimension    = GridLayout::dimension;
    static constexpr auto interp_order = GridLayout::interp_order;

    MHDDiagnosticWriter(H5Writer& h5Writer)
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
    auto isActiveDiag(DiagnosticProperties& diagnostic, std::string const& tree,
                      std::string const& name) const
    {
        return diagnostic.quantity == tree + name;
    }
};

template<typename H5Writer>
void MHDDiagnosticWriter<H5Writer>::createFiles(DiagnosticProperties& diagnostic)
{
    std::string tree{"/mhd/"};
    checkCreateFileFor_(diagnostic, fileData_, tree, "density", "velocity", "B", "pressure", "rhoV",
                        "Etot");
}

template<typename H5Writer>
void MHDDiagnosticWriter<H5Writer>::compute(DiagnosticProperties& diagnostic)
{
    core::ToPrimitiveConverter<GridLayout> toPrim;

    auto& h5Writer  = this->h5Writer_;
    auto& modelView = h5Writer.modelView();
    auto minLvl     = h5Writer.minLevel;
    auto maxLvl     = h5Writer.maxLevel;

    auto& rho  = modelView.getRho();
    auto& V    = modelView.getV();
    auto& B    = modelView.getB();
    auto& P    = modelView.getP();
    auto& rhoV = modelView.getRhoV();
    auto& Etot = modelView.getEtot();

    std::string tree{"/mhd/"};
    if (isActiveDiag(diagnostic, tree, "velocity"))
    {
        auto computeVelocity = [&](GridLayout& layout, std::string& patchID, std::size_t iLevel) {
            auto _sl = core::SetLayout(&layout, toPrim);
            toPrim.rhoVToV(rho, rhoV, V);
        };
        modelView.visitHierarchy(computeVelocity, minLvl, maxLvl);
    }
    if (isActiveDiag(diagnostic, tree, "pressure"))
    {
        auto computePressure = [&](GridLayout& layout, std::string& patchID, std::size_t iLevel) {
            auto _sl = core::SetLayout(&layout, toPrim);
            toPrim.eosEtotToP(rho, rhoV, B, Etot, P);
        };
        modelView.visitHierarchy(computePressure, minLvl, maxLvl);
    }
}

template<typename H5Writer>
void MHDDiagnosticWriter<H5Writer>::getDataSetInfo(DiagnosticProperties& diagnostic,
                                                   std::size_t iLevel, std::string const& patchID,
                                                   Attributes& patchAttributes)
{
    auto& h5Writer         = this->h5Writer_;
    auto& rho              = h5Writer.modelView().getRho();
    auto& V                = h5Writer.modelView().getV();
    auto& B                = h5Writer.modelView().getB();
    auto& P                = h5Writer.modelView().getP();
    auto& rhoV             = h5Writer.modelView().getRhoV();
    auto& Etot             = h5Writer.modelView().getEtot();
    std::string lvlPatchID = std::to_string(iLevel) + "_" + patchID;

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

    std::string tree{"/mhd/"};
    if (isActiveDiag(diagnostic, tree, "density"))
        infoDS(rho, "density", patchAttributes[lvlPatchID]);
    if (isActiveDiag(diagnostic, tree, "velocity"))
        infoVF(V, "velocity", patchAttributes[lvlPatchID]);
    if (isActiveDiag(diagnostic, tree, "B"))
        infoVF(B, "B", patchAttributes[lvlPatchID]);
    if (isActiveDiag(diagnostic, tree, "pressure"))
        infoDS(P, "pressure", patchAttributes[lvlPatchID]);
    if (isActiveDiag(diagnostic, tree, "rhoV"))
        infoVF(rhoV, "rhoV", patchAttributes[lvlPatchID]);
    if (isActiveDiag(diagnostic, tree, "Etot"))
        infoDS(Etot, "Etot", patchAttributes[lvlPatchID]);
}

template<typename H5Writer>
void MHDDiagnosticWriter<H5Writer>::initDataSets(
    DiagnosticProperties& diagnostic,
    std::unordered_map<std::size_t, std::vector<std::string>> const& patchIDs,
    Attributes& patchAttributes, std::size_t maxLevel)
{
    auto& h5Writer = this->h5Writer_;
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

    auto initPatch = [&](auto& lvl, auto& attr, std::string patchID = "") {
        bool null        = patchID.empty();
        std::string path = h5Writer.getPatchPathAddTimestamp(lvl, patchID) + "/";

        std::string tree{"/mhd/"};
        if (isActiveDiag(diagnostic, tree, "density"))
            initDS(path, attr["mhd"], "density", null);
        if (isActiveDiag(diagnostic, tree, "velocity"))
            initVF(path, attr["mhd"], "velocity", null);
        if (isActiveDiag(diagnostic, tree, "B"))
            initVF(path, attr["mhd"], "B", null);
        if (isActiveDiag(diagnostic, tree, "pressure"))
            initDS(path, attr["mhd"], "pressure", null);
        if (isActiveDiag(diagnostic, tree, "rhoV"))
            initVF(path, attr["mhd"], "rhoV", null);
        if (isActiveDiag(diagnostic, tree, "Etot"))
            initDS(path, attr["mhd"], "Etot", null);
    };

    initDataSets_(patchIDs, patchAttributes, maxLevel, initPatch);
}

template<typename H5Writer>
void MHDDiagnosticWriter<H5Writer>::write(DiagnosticProperties& diagnostic)
{
    auto& h5Writer = this->h5Writer_;
    auto& rho      = h5Writer.modelView().getRho();
    auto& V        = h5Writer.modelView().getV();
    auto& B        = h5Writer.modelView().getB();
    auto& P        = h5Writer.modelView().getP();
    auto& rhoV     = h5Writer.modelView().getRhoV();
    auto& Etot     = h5Writer.modelView().getEtot();
    auto& h5file   = *fileData_.at(diagnostic.quantity);

    auto writeDS = [&](auto path, auto& field) {
        h5file.template write_data_set_flat<GridLayout::dimension>(path, field.data());
    };
    auto writeTF
        = [&](auto path, auto& vecF) { h5Writer.writeTensorFieldAsDataset(h5file, path, vecF); };

    std::string path = h5Writer.patchPath() + "/";
    std::string tree{"/mhd/"};
    if (isActiveDiag(diagnostic, tree, "density"))
        writeDS(path + "density", rho);
    if (isActiveDiag(diagnostic, tree, "velocity"))
        writeTF(path + "velocity", V);
    if (isActiveDiag(diagnostic, tree, "B"))
        writeTF(path + "B", B);
    if (isActiveDiag(diagnostic, tree, "pressure"))
        writeDS(path + "pressure", P);
    if (isActiveDiag(diagnostic, tree, "rhoV"))
        writeTF(path + "rhoV", rhoV);
    if (isActiveDiag(diagnostic, tree, "Etot"))
        writeDS(path + "Etot", Etot);
}

template<typename H5Writer>
void MHDDiagnosticWriter<H5Writer>::writeAttributes(
    DiagnosticProperties& diagnostic, Attributes& fileAttributes,
    std::unordered_map<std::size_t, std::vector<std::pair<std::string, Attributes>>>&
        patchAttributes,
    std::size_t maxLevel)
{
    writeAttributes_(diagnostic, *fileData_.at(diagnostic.quantity), fileAttributes,
                     patchAttributes, maxLevel);
}

} // namespace PHARE::diagnostic::h5


#endif
