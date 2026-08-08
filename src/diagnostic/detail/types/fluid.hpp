#ifndef PHARE_DIAGNOSTIC_DETAIL_TYPES_FLUID_HPP
#define PHARE_DIAGNOSTIC_DETAIL_TYPES_FLUID_HPP

#include "amr/physical_models/hybrid_model.hpp"
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
    using Attributes       = Super::Attributes;
    using ModelViewVariant = Super::ModelViewVariant;
    using ModelMapper_t    = Super::ModelMapper_t;
    using FloatType        = H5Writer::FloatType;

    static constexpr auto dimension = H5Writer::dimension;


    FluidDiagnosticWriter(H5Writer& h5Writer)
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
};


template<typename H5Writer>
void FluidDiagnosticWriter<H5Writer>::compute(DiagnosticProperties& diagnostic)
{
    using HybridModelView_t = std::decay_t<decltype(this->h5Writer_.mapper().hyridModelView())>;
    using Model_t           = HybridModelView_t::Model_t;

    if constexpr (solver::is_hybrid_model_v<Model_t>)
    {
        using Computers = HybridModelView_t::FluidComputers;

        auto const& qty = diagnostic.quantity;
        auto& h5Writer  = this->h5Writer_;
        auto& modelView = h5Writer.mapper().hyridModelView();
        auto& computers = Computers::getOrCreateFor(modelView);

        computers.compute(qty, modelView, h5Writer.minLevel, h5Writer.maxLevel,
                          h5Writer.timestamp());
    }
}



template<typename H5Writer>
void FluidDiagnosticWriter<H5Writer>::createFiles(DiagnosticProperties& diagnostic)
{
    for (auto const& pop : this->h5Writer_.mapper().hyridModelView().getIons())
    {
        std::string tree{"/ions/pop/" + pop.name() + "/"};
        checkCreateFileFor_(diagnostic, fileData_, tree, "density", "charge_density", "flux",
                            "momentum_tensor", "kinetic_energy_flux_vector");
    }

    std::string tree{"/ions/"};
    checkCreateFileFor_(diagnostic, fileData_, tree, "charge_density", "mass_density",
                        "bulkVelocity", "momentum_tensor", "kinetic_energy_flux_vector");
}




template<typename H5Writer>
void FluidDiagnosticWriter<H5Writer>::getDataSetInfo(DiagnosticProperties& diagnostic,
                                                     std::size_t iLevel, std::string const& patchID,
                                                     Attributes& patchAttributes,
                                                     ModelViewVariant& modelViewVariant)
{
    std::visit(
        [&](auto& mv) {
            using ModelView_t = std::decay_t<decltype(mv)>;
            using Model_t     = ModelView_t::Model_t;

            // fluid diagnostics only ever concern the hybrid model at present; modelView can
            // legitimately hold the mhd alternative here (diagnostics_ isn't filtered per-model
            // before dispatch), so this is a real per-instantiation branch, not dead code
            if constexpr (solver::is_hybrid_model_v<Model_t>)
            {
                using GridLayout = ModelView_t::GridLayout;
                using Accessors  = ModelView_t::FluidAccessors;

                auto& model     = mv.model();
                auto& accessors = Accessors::getOrCreateFor(mv);
                auto const& qty = diagnostic.quantity;
                std::string lvlPatchID{std::to_string(iLevel) + "_" + patchID};

                auto const infoDS = [](auto& field, std::string const& name, auto& attr) {
                    auto ghosts = GridLayout::nDNbrGhosts(field.physicalQuantity());
                    for (std::uint8_t i = 1; i < GridLayout::dimension; ++i)
                        if (ghosts[i] != ghosts[i - 1])
                            throw std::runtime_error("ghosts per direction must be constant");
                    attr[name + "_ghosts"] = static_cast<std::size_t>(ghosts[0]);

                    // highfive doesn't accept uint32 which ndarray.shape() is
                    auto const& shape = field.shape();
                    attr[name]
                        = std::vector<std::size_t>(shape.data(), shape.data() + shape.size());
                };

                auto const action = [&](auto& quantity, std::string const& name, auto& attr) {
                    using Quantity = std::decay_t<decltype(quantity)>;
                    if constexpr (std::is_same_v<Quantity, typename Accessors::Field_t>)
                        infoDS(quantity, name, attr);
                    else
                        for (auto const& [id, type] :
                             core::Components::componentMap<Quantity::rank>())
                            infoDS(quantity.getComponent(type), name + "_" + id, attr);
                };

                accessors.dispatch(
                    qty, model,
                    [&](auto& field, std::string const& name, std::string const& ownerKey) {
                        action(field, name, patchAttributes[lvlPatchID][ownerKey]);
                    });
            }
        },
        modelViewVariant);
}




template<typename H5Writer>
void FluidDiagnosticWriter<H5Writer>::initDataSets(
    DiagnosticProperties& diagnostic,
    std::unordered_map<std::size_t, std::vector<std::string>> const& patchIDs,
    Attributes& patchAttributes, std::size_t maxLevel)
{
    using HybridModelView_t = std::decay_t<decltype(this->h5Writer_.mapper().hyridModelView())>;
    using GridLayout        = HybridModelView_t::GridLayout;
    using Accessors         = HybridModelView_t::FluidAccessors;

    auto& h5Writer  = this->h5Writer_;
    auto& modelView = h5Writer.mapper().hyridModelView();
    auto& model     = modelView.model();
    auto& accessors = Accessors::getOrCreateFor(modelView);
    auto& h5file    = Super::h5FileForQuantity(diagnostic);
    auto const& qty = diagnostic.quantity;

    auto const initDS = [&](auto& path, auto& attr, std::string const& key, bool null) {
        auto dsPath = path + key;
        h5Writer.template createDataSet<FloatType>(
            h5file, dsPath,
            null ? std::vector<std::size_t>(GridLayout::dimension, 0)
                 : attr[key].template to<std::vector<std::size_t>>());
        this->writeGhostsAttr_(h5file, dsPath,
                               null ? 0 : attr[key + "_ghosts"].template to<std::size_t>(), null);
    };

    auto const initPatch = [&](auto& lvl, auto& attr, std::string patchID = "") {
        bool null        = patchID.empty();
        std::string path = h5Writer.getPatchPathAddTimestamp(lvl, patchID) + "/";

        auto const action = [&](auto& field, std::string const& name, auto& fieldAttr) {
            using Field = std::decay_t<decltype(field)>;
            if constexpr (std::is_same_v<Field, typename Accessors::Field_t>)
                initDS(path, fieldAttr, name, null);
            else
                for (auto& [id, type] : core::Components::componentMap<Field::rank>())
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
void FluidDiagnosticWriter<H5Writer>::write(DiagnosticProperties& diagnostic)
{
    using HybridModelView_t = std::decay_t<decltype(this->h5Writer_.mapper().hyridModelView())>;
    using GridLayout        = HybridModelView_t::GridLayout;
    using Accessors         = HybridModelView_t::FluidAccessors;

    auto& h5Writer   = this->h5Writer_;
    auto& modelView  = h5Writer.mapper().hyridModelView();
    auto& model      = modelView.model();
    auto& accessors  = Accessors::getOrCreateFor(modelView);
    auto& h5file     = Super::h5FileForQuantity(diagnostic);
    auto const& qty  = diagnostic.quantity;
    std::string path = h5Writer.patchPath() + "/";

    auto const action = [&](auto& field, std::string const& name) {
        using Field = std::decay_t<decltype(field)>;
        if constexpr (std::is_same_v<Field, typename Accessors::Field_t>)
            h5file.template write_data_set_flat<GridLayout::dimension>(path + name, field.data());
        else
            H5Writer::writeTensorFieldAsDataset(h5file, path + name, field);
    };

    accessors.dispatch(qty, model, [&](auto& field, std::string const& name, std::string const&) {
        action(field, name);
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

    for (auto& pop : h5Writer.mapper().hyridModelView().getIons())
    {
        std::string tree = "/ions/pop/" + pop.name() + "/";
        checkWrite(tree, "density", pop);
        checkWrite(tree, "charge_density", pop);
        checkWrite(tree, "flux", pop);
        checkWrite(tree, "momentum_tensor", pop);
        checkWrite(tree, "kinetic_energy_flux_vector", pop);
    }

    writeAttributes_(diagnostic, h5file, fileAttributes, patchAttributes, maxLevel);
}

} // namespace PHARE::diagnostic::h5

#endif /* PHARE_DIAGNOSTIC_DETAIL_TYPES_FLUID_H */
