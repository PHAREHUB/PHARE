#ifndef PHARE_DIAGNOSTIC_DETAIL_TYPES_MHD_HPP
#define PHARE_DIAGNOSTIC_DETAIL_TYPES_MHD_HPP

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
    using Attributes       = typename Super::Attributes;
    using ModelViewVariant = typename Super::ModelViewVariant;
    using MHDModel_t       = typename H5Writer::ModelMapper_t::HasModels_t::MHDModel_t;
    using GridLayout       = typename MHDModel_t::gridlayout_type;
    using FloatType        = typename H5Writer::FloatType;

    static constexpr auto dimension = GridLayout::dimension;

    MHDDiagnosticWriter(H5Writer& h5Writer)
        : Super{h5Writer}
    {
    }
    void write(DiagnosticProperties&) override;
    void compute(DiagnosticProperties&) override;

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
    using ModelView_t = std::decay_t<decltype(this->h5Writer_.mapper().mhdModelView())>;
    using Computers   = typename ModelView_t::MHDComputers;

    auto& h5Writer  = this->h5Writer_;
    auto& modelView = h5Writer.mapper().mhdModelView();
    auto& computers = Computers::getOrCreateFor(modelView);
    // or FloatType if we want to expose that to DiagnosticProperties
    double const gamma = diagnostic.fileAttributes["heat_capacity_ratio"];

    computers.compute(diagnostic.quantity, modelView, gamma, h5Writer.minLevel, h5Writer.maxLevel,
                      h5Writer.timestamp());
}

template<typename H5Writer>
void MHDDiagnosticWriter<H5Writer>::getDataSetInfo(DiagnosticProperties& diagnostic,
                                                   std::size_t iLevel, std::string const& patchID,
                                                   Attributes& patchAttributes,
                                                   ModelViewVariant& /*modelView*/)
{
    using ModelView_t = std::decay_t<decltype(this->h5Writer_.mapper().mhdModelView())>;
    using Accessors   = typename ModelView_t::MHDAccessors;

    auto& h5Writer         = this->h5Writer_;
    auto& modelView         = h5Writer.mapper().mhdModelView();
    auto& model             = modelView.model();
    auto& accessors         = Accessors::getOrCreateFor(modelView);
    auto const& qty         = diagnostic.quantity;
    std::string lvlPatchID = std::to_string(iLevel) + "_" + patchID;

    auto setGhostNbr = [](auto const& field, auto& attr, auto const& name) {
        auto ghosts              = GridLayout::nDNbrGhosts(field.physicalQuantity());
        attr[name + "_ghosts_x"] = static_cast<std::size_t>(ghosts[0]);
        if constexpr (GridLayout::dimension > 1)
            attr[name + "_ghosts_y"] = static_cast<std::size_t>(ghosts[1]);
        if constexpr (GridLayout::dimension > 2)
            attr[name + "_ghosts_z"] = static_cast<std::size_t>(ghosts[2]);
    };

    auto infoDS = [&](auto& field, std::string const& name, auto& attr) {
        // highfive doesn't accept uint32 which ndarray.shape() is
        auto const& shape = field.shape();
        attr[name]        = std::vector<std::size_t>(shape.data(), shape.data() + shape.size());
        setGhostNbr(field, attr, name);
    };

    auto const action = [&](auto& field, std::string const& name, auto& attr) {
        using Quantity = std::decay_t<decltype(field)>;
        if constexpr (std::is_same_v<Quantity, typename Accessors::Field_t>)
            infoDS(field, name, attr);
        else
            for (auto const& [id, type] : core::VectorComponents::map())
                infoDS(field.getComponent(type), name + "_" + id, attr);
    };

    accessors.dispatch(qty, model,
                       [&](auto& field, std::string const& name, std::string const& ownerKey) {
                           action(field, name, patchAttributes[lvlPatchID][ownerKey]);
                       });
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

    using ModelView_t = std::decay_t<decltype(h5Writer.mapper().mhdModelView())>;
    using Accessors   = typename ModelView_t::MHDAccessors;

    auto& modelView = h5Writer.mapper().mhdModelView();
    auto& model     = modelView.model();
    auto& accessors = Accessors::getOrCreateFor(modelView);
    auto const& qty = diagnostic.quantity;

    auto initPatch = [&](auto& lvl, auto& attr, std::string patchID = "") {
        bool null        = patchID.empty();
        std::string path = h5Writer.getPatchPathAddTimestamp(lvl, patchID) + "/";

        auto const action = [&](auto& field, std::string const& name, auto& fieldAttr) {
            using Quantity = std::decay_t<decltype(field)>;
            if constexpr (std::is_same_v<Quantity, typename Accessors::Field_t>)
                initDS(path, fieldAttr, name, null);
            else
                for (auto& [id, type] : core::VectorComponents::map())
                    initDS(path, fieldAttr, name + "_" + id, null);
        };

        accessors.dispatch(qty, model,
                           [&](auto& field, std::string const& name, std::string const& ownerKey) {
                               action(field, name, attr[ownerKey]);
                           });
    };

    initDataSets_(patchIDs, patchAttributes, maxLevel, initPatch);
}

template<typename H5Writer>
void MHDDiagnosticWriter<H5Writer>::write(DiagnosticProperties& diagnostic)
{
    using ModelView_t = std::decay_t<decltype(this->h5Writer_.mapper().mhdModelView())>;
    using Accessors   = typename ModelView_t::MHDAccessors;

    auto& h5Writer   = this->h5Writer_;
    auto& modelView  = h5Writer.mapper().mhdModelView();
    auto& model      = modelView.model();
    auto& accessors  = Accessors::getOrCreateFor(modelView);
    auto& h5file     = *fileData_.at(diagnostic.quantity);
    auto const& qty  = diagnostic.quantity;
    std::string path = h5Writer.patchPath() + "/";

    auto const action = [&](auto& field, std::string const& name) {
        using Quantity = std::decay_t<decltype(field)>;
        if constexpr (std::is_same_v<Quantity, typename Accessors::Field_t>)
            h5file.template write_data_set_flat<GridLayout::dimension>(path + name, field.data());
        else
            h5Writer.writeTensorFieldAsDataset(h5file, path + name, field);
    };

    accessors.dispatch(qty, model, [&](auto& field, std::string const& name, std::string const&) {
        action(field, name);
    });
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
