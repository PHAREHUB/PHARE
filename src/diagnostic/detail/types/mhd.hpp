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
    checkCreateFileFor_(diagnostic, fileData_, tree, "rho", "V", "P", "rhoV", "Etot");
}

template<typename H5Writer>
void MHDDiagnosticWriter<H5Writer>::compute(DiagnosticProperties& diagnostic)
{
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
    if (isActiveDiag(diagnostic, tree, "V"))
    {
        auto computeVelocity = [&](GridLayout& layout, std::string patchID, std::size_t iLevel) {
            core::ToPrimitiveConverter_ref<GridLayout> toPrim{layout};
            toPrim.rhoVToVOnGhostBox(rho, rhoV, V);
        };
        modelView.visitHierarchy(computeVelocity, minLvl, maxLvl);
    }
    if (isActiveDiag(diagnostic, tree, "P"))
    {
        auto computePressure = [&](GridLayout& layout, std::string patchID, std::size_t iLevel) {
            auto const gamma = diagnostic.fileAttributes["heat_capacity_ratio"]
                                   .template to<double>(); // or FloatType if we want to expose that
                                                           // to DiagnosticProperties
            core::ToPrimitiveConverter_ref<GridLayout> toPrim{layout};
            toPrim.eosEtotToPOnGhostBox(gamma, rho, rhoV, B, Etot, P);
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
    if (isActiveDiag(diagnostic, tree, "rho"))
        infoDS(rho, "rho", patchAttributes[lvlPatchID]["mhd"]);
    if (isActiveDiag(diagnostic, tree, "V"))
        infoVF(V, "V", patchAttributes[lvlPatchID]["mhd"]);
    if (isActiveDiag(diagnostic, tree, "P"))
        infoDS(P, "P", patchAttributes[lvlPatchID]["mhd"]);
    if (isActiveDiag(diagnostic, tree, "rhoV"))
        infoVF(rhoV, "rhoV", patchAttributes[lvlPatchID]["mhd"]);
    if (isActiveDiag(diagnostic, tree, "Etot"))
        infoDS(Etot, "Etot", patchAttributes[lvlPatchID]["mhd"]);
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
        if (isActiveDiag(diagnostic, tree, "rho"))
            initDS(path, attr["mhd"], "rho", null);
        if (isActiveDiag(diagnostic, tree, "V"))
            initVF(path, attr["mhd"], "V", null);
        if (isActiveDiag(diagnostic, tree, "P"))
            initDS(path, attr["mhd"], "P", null);
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
    auto& P        = h5Writer.modelView().getP();
    auto& rhoV     = h5Writer.modelView().getRhoV();
    auto& Etot     = h5Writer.modelView().getEtot();
    auto& h5file   = *fileData_.at(diagnostic.quantity);

    auto hasNaN = [](auto const& container) {
        return std::any_of(container.begin(), container.end(),
                           [](auto const& x) { return std::isnan(x); });
    };

    auto checkNaN = [&](std::string const& name, auto const& field) {
        if (hasNaN(field))
        {
            throw std::runtime_error("NaN detected in field '" + name + "'");
        }
    };

    auto writeDS = [&](auto path, auto& field) {
        h5file.template write_data_set_flat<GridLayout::dimension>(path, field.data());
        checkNaN(path, field);
    };

    auto writeTF = [&](auto path, auto& vecF) {
        h5Writer.writeTensorFieldAsDataset(h5file, path, vecF);
        for (std::size_t d = 0; d < vecF.size(); ++d)
            checkNaN(path + "[" + std::to_string(d) + "]", vecF[d]);
    };

    std::string path = h5Writer.patchPath() + "/";
    std::string tree{"/mhd/"};

    if (isActiveDiag(diagnostic, tree, "rho"))
        writeDS(path + "rho", rho);
    if (isActiveDiag(diagnostic, tree, "V"))
        writeTF(path + "V", V);
    if (isActiveDiag(diagnostic, tree, "P"))
        writeDS(path + "P", P);
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
