#ifndef PHARE_DIAGNOSTIC_DETAIL_TYPES_FLUID_HPP
#define PHARE_DIAGNOSTIC_DETAIL_TYPES_FLUID_HPP

#include "core/data/vecfield/vecfield_component.hpp"

#include "diagnostic/detail/h5typewriter.hpp"

#include <stdexcept>

namespace PHARE::diagnostic::h5
{
/*
 * It is assumed that each patch has equal number of populations
 *
 * Possible outputs
 *
 * /t#/pl#/p#/ions/density
 * /t#/pl#/p#/ions/bulkVelocity/(x,y,z)
 * /t#/pl#/p#/ions/charge_flux/(x,y,z)
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
    using Attributes = Super::Attributes;
    using GridLayout = H5Writer::GridLayout;
    using FloatType  = H5Writer::FloatType;
    using ModelView  = std::decay_t<decltype(std::declval<H5Writer&>().modelView())>;
    using Accessors  = ModelView::FluidAccessors;
    using Computers  = ModelView::FluidComputers;

    static constexpr auto dimension = GridLayout::dimension;


    FluidDiagnosticWriter(H5Writer& h5Writer)
        : Super{h5Writer}
    {
    }

    void write(DiagnosticProperties&) override;
    void compute_as_needed(DiagnosticProperties&) override;

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
};




template<typename H5Writer>
void FluidDiagnosticWriter<H5Writer>::compute_as_needed(DiagnosticProperties& diagnostic)
{
    auto& h5Writer  = this->h5Writer_;
    auto& modelView = h5Writer.modelView();

    Computers::getOrCreateFor(modelView).compute(diagnostic.quantity, modelView, h5Writer.minLevel,
                                                 h5Writer.maxLevel, h5Writer.timestamp());
}



template<typename H5Writer>
void FluidDiagnosticWriter<H5Writer>::createFiles(DiagnosticProperties& diagnostic)
{
    for (auto const& pop : this->h5Writer_.modelView().getIons())
    {
        std::string tree{"/ions/pop/" + pop.name() + "/"};
        checkCreateFileFor_(diagnostic, fileData_, tree, "density", "charge_density", "flux",
                            "momentum_tensor");
    }

    std::string tree{"/ions/"};
    checkCreateFileFor_(diagnostic, fileData_, tree, "charge_density", "mass_density",
                        "bulkVelocity", "charge_flux", "momentum_tensor");
}




template<typename H5Writer>
void FluidDiagnosticWriter<H5Writer>::getDataSetInfo(DiagnosticProperties& diagnostic,
                                                     std::size_t iLevel, std::string const& patchID,
                                                     Attributes& patchAttributes)
{
    auto& modelView = this->h5Writer_.modelView();
    std::string lvlPatchID{std::to_string(iLevel) + "_" + patchID};

    auto const infoDS = [](auto& field, std::string const& name, auto& attr) {
        auto ghosts = GridLayout::nDNbrGhosts(field.physicalQuantity());
        for (std::uint8_t i = 1; i < GridLayout::dimension; ++i)
            if (ghosts[i] != ghosts[i - 1])
                throw std::runtime_error("ghosts per direction must be constant");
        attr[name + "_ghosts"] = static_cast<std::size_t>(ghosts[0]);

        // highfive doesn't accept uint32 which ndarray.shape() is
        auto const& shape = field.shape();
        attr[name]        = std::vector<std::size_t>(shape.data(), shape.data() + shape.size());
    };

    auto const action = [&](auto& quantity, std::string const& name, auto& attr) {
        using Quantity = std::decay_t<decltype(quantity)>;
        if constexpr (std::is_same_v<Quantity, typename Accessors::Field_t>)
            infoDS(quantity, name, attr);
        else
            for (auto const& [id, type] : core::Components::componentMap<Quantity::rank>())
                infoDS(quantity.getComponent(type), name + "_" + id, attr);
    };

    Accessors::getOrCreateFor(modelView).dispatch(
        diagnostic.quantity, modelView,
        [&](auto& field, std::string const& name, std::string const& ownerKey) {
            action(field, name, patchAttributes[lvlPatchID][ownerKey]);
        });
}




template<typename H5Writer>
void FluidDiagnosticWriter<H5Writer>::initDataSets(
    DiagnosticProperties& diagnostic,
    std::unordered_map<std::size_t, std::vector<std::string>> const& patchIDs,
    Attributes& patchAttributes, std::size_t maxLevel)
{
    auto& h5Writer  = this->h5Writer_;
    auto& modelView = h5Writer.modelView();
    auto& accessors = Accessors::getOrCreateFor(modelView);
    auto& h5file    = Super::h5FileForQuantity(diagnostic);

    auto const initDS = [&](auto& path, auto& attr, std::string const& key, bool const null) {
        auto dsPath = path + key;
        h5Writer.template createDataSet<FloatType>(
            h5file, dsPath,
            null ? std::vector<std::size_t>(GridLayout::dimension, 0)
                 : attr[key].template to<std::vector<std::size_t>>());
        this->writeGhostsAttr_(h5file, dsPath,
                               null ? 0 : attr[key + "_ghosts"].template to<std::size_t>(), null);
    };

    auto const initPatch = [&](auto& lvl, auto& attr, std::string patchID = "") {
        bool const null  = patchID.empty();
        std::string path = h5Writer.getPatchPathAddTimestamp(lvl, patchID) + "/";

        accessors.dispatch(diagnostic.quantity, modelView,
                           [&](auto& field, std::string const& name, std::string const& ownerKey) {
                               using Field = std::decay_t<decltype(field)>;
                               if constexpr (std::is_same_v<Field, typename Accessors::Field_t>)
                                   initDS(path, attr[ownerKey], name, null);
                               else
                                   for (auto& [id, type] :
                                        core::Components::componentMap<Field::rank>())
                                       initDS(path, attr[ownerKey], name + "_" + id, null);
                           });
    };

    initDataSets_(patchIDs, patchAttributes, maxLevel, initPatch);
}


template<typename H5Writer>
void FluidDiagnosticWriter<H5Writer>::write(DiagnosticProperties& diagnostic)
{
    auto& h5Writer   = this->h5Writer_;
    auto& modelView  = h5Writer.modelView();
    auto& h5file     = Super::h5FileForQuantity(diagnostic);
    std::string path = h5Writer.patchPath() + "/";

    Accessors::getOrCreateFor(modelView).dispatch(
        diagnostic.quantity, modelView, [&](auto& field, std::string const& name, auto const&) {
            using Field = std::decay_t<decltype(field)>;
            if constexpr (std::is_same_v<Field, typename Accessors::Field_t>)
                h5file.template write_data_set_flat<GridLayout::dimension>(path + name,
                                                                           field.data());
            else
                h5Writer.writeTensorFieldAsDataset(h5file, path + name, field);
        });
}


template<typename H5Writer>
void FluidDiagnosticWriter<H5Writer>::writeAttributes(
    DiagnosticProperties& diagnostic, Attributes& fileAttributes,
    std::unordered_map<std::size_t, std::vector<std::pair<std::string, Attributes>>>&
        patchAttributes,
    std::size_t maxLevel)
{
    auto& h5Writer = this->h5Writer_;
    auto& h5file   = Super::h5FileForQuantity(diagnostic);

    auto checkWrite = [&](auto& tree, std::string qty, auto const& pop) {
        if (diagnostic.quantity == tree + qty)
            this->writeIonPopAttributes_(h5file, pop);
    };

    for (auto& pop : h5Writer.modelView().getIons())
    {
        std::string tree = "/ions/pop/" + pop.name() + "/";
        checkWrite(tree, "density", pop);
        checkWrite(tree, "charge_density", pop);
        checkWrite(tree, "flux", pop);
        checkWrite(tree, "momentum_tensor", pop);
    }

    writeAttributes_(diagnostic, h5file, fileAttributes, patchAttributes, maxLevel);
}

} // namespace PHARE::diagnostic::h5

#endif /* PHARE_DIAGNOSTIC_DETAIL_TYPES_FLUID_H */
