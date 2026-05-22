#ifndef PHARE_DIAGNOSTIC_DETAIL_TYPES_FLUID_HPP
#define PHARE_DIAGNOSTIC_DETAIL_TYPES_FLUID_HPP

#include "amr/physical_models/hybrid_model.hpp"
#include "core/data/vecfield/vecfield_component.hpp"

#include "diagnostic/detail/h5typewriter.hpp"
#include "diagnostic/computers/diagnostic_computers.hpp"

#include <memory>
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
    using Attributes = Super::Attributes;
    using GridLayout = H5Writer::GridLayout;
    using FloatType  = H5Writer::FloatType;

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
    struct HybridFluidComputers;


    auto isActiveDiag(DiagnosticProperties const& diagnostic, std::string const& tree,
                      std::string var)
    {
        return diagnostic.quantity == tree + var;
    };


    std::unique_ptr<HybridFluidComputers> hybridComputers_;
};

template<typename H5Writer>
struct FluidDiagnosticWriter<H5Writer>::HybridFluidComputers
{
    using Model_t         = H5Writer::ModelView_t::Model_t;
    using IonPopulation_t = Model_t::ions_type::value_type;
    using IonsFunctor     = std::function<void(H5Writer&)>;
    using PopFunctor      = std::function<void(H5Writer&, IonPopulation_t&)>;

    static HybridFluidComputers& getOrCreateFor(auto& fluid_writer)
    {
        if (!fluid_writer.hybridComputers_)
            fluid_writer.hybridComputers_ = std::make_unique<HybridFluidComputers>(fluid_writer);
        return *fluid_writer.hybridComputers_;
    }

    HybridFluidComputers(auto& fluid_writer)
    {
        for (auto& pop : fluid_writer.h5Writer_.modelView().getIons())
        {
            pop_functors["/ions/pop/" + pop.name() + "/momentum_tensor"]
                = [](auto&&... args) { compute_pop_momentum_tensor(args...); };
            pop_functors["/ions/pop/" + pop.name() + "/kinetic_energy_flux_vector"]
                = [](auto&&... args) { compute_pop_kinetic_energy_flux_vector(args...); };
        }
    }

    std::unordered_map<std::string, IonsFunctor> ion_functors{{
        "/ions/momentum_tensor",
        [](H5Writer& h5Writer) { compute_momentum_tensor(h5Writer); },
    }};
    std::unordered_map<std::string, PopFunctor> pop_functors;
};




template<typename H5Writer>
void FluidDiagnosticWriter<H5Writer>::compute(DiagnosticProperties& diagnostic)
{
    using Model_t = H5Writer::ModelView_t::Model_t;

    auto const qty = diagnostic.quantity;
    if constexpr (solver::is_hybrid_model_v<Model_t>)
    {
        auto& computers     = HybridFluidComputers::getOrCreateFor(*this);
        bool const ion_func = computers.ion_functors.contains(qty);
        bool const pop_func = computers.pop_functors.contains(qty);

        if (!ion_func and !pop_func)
            return;

        if (ion_func)
            computers.ion_functors[qty](this->h5Writer_);

        else
            for (auto& pop : this->h5Writer_.modelView().getIons())
                if (auto key = "/" + pop.name() + "/"; qty.find(key) != std::string::npos)
                {
                    computers.pop_functors[qty](this->h5Writer_, pop);
                    return;
                }
    }
}



template<typename H5Writer>
void FluidDiagnosticWriter<H5Writer>::createFiles(DiagnosticProperties& diagnostic)
{
    for (auto const& pop : this->h5Writer_.modelView().getIons())
    {
        std::string tree{"/ions/pop/" + pop.name() + "/"};
        checkCreateFileFor_(diagnostic, fileData_, tree, "density", "charge_density", "flux",
                            "momentum_tensor", "kinetic_energy_flux_vector");
    }

    std::string tree{"/ions/"};
    checkCreateFileFor_(diagnostic, fileData_, tree, "charge_density", "mass_density",
                        "bulkVelocity", "momentum_tensor");
}




template<typename H5Writer>
void FluidDiagnosticWriter<H5Writer>::getDataSetInfo(DiagnosticProperties& diagnostic,
                                                     std::size_t iLevel, std::string const& patchID,
                                                     Attributes& patchAttributes)
{
    auto& h5Writer = this->h5Writer_;
    auto& ions     = h5Writer.modelView().getIons();
    std::string lvlPatchID{std::to_string(iLevel) + "_" + patchID};


    auto const setGhostNbr = [](auto const& field, auto& attr, auto const& name) {
        auto ghosts = GridLayout::nDNbrGhosts(field.physicalQuantity());
        for (std::uint8_t i = 1; i < GridLayout::dimension; ++i)
            if (ghosts[i] != ghosts[i - 1])
                throw std::runtime_error("ghosts per direction must be constant");
        attr[name + "_ghosts"] = static_cast<std::size_t>(ghosts[0]);
    };

    auto const infoDS = [&](auto& field, std::string name, auto& attr) {
        // highfive doesn't accept uint32 which ndarray.shape() is
        auto const& shape = field.shape();
        attr[name]        = std::vector<std::size_t>(shape.data(), shape.data() + shape.size());
        setGhostNbr(field, attr, name);
    };

    auto const infoVF = [&](auto& vecF, std::string name, auto& attr) {
        for (auto const& [id, type] : core::VectorComponents::map())
            infoDS(vecF.getComponent(type), name + "_" + id, attr);
    };

    auto const infoTF = [&](auto& tensorF, std::string name, auto& attr) {
        for (auto const& [id, type] : core::TensorComponents::map())
            infoDS(tensorF.getComponent(type), name + "_" + id, attr);
    };

    for (auto& pop : ions)
    {
        std::string tree{"/ions/pop/" + pop.name() + "/"};
        auto& popAttr = patchAttributes[lvlPatchID]["fluid_" + pop.name()];
        if (isActiveDiag(diagnostic, tree, "density"))
            infoDS(pop.particleDensity(), "density", popAttr);
        if (isActiveDiag(diagnostic, tree, "charge_density"))
            infoDS(pop.chargeDensity(), "charge_density", popAttr);
        if (isActiveDiag(diagnostic, tree, "flux"))
            infoVF(pop.flux(), "flux", popAttr);
        if (isActiveDiag(diagnostic, tree, "momentum_tensor"))
            infoTF(pop.momentumTensor(), "momentum_tensor", popAttr);
        if (isActiveDiag(diagnostic, tree, "kinetic_energy_flux_vector"))
            infoVF(pop.kineticEnergyFlux(), "kinetic_energy_flux_vector", popAttr);
    }

    std::string tree{"/ions/"};
    if (isActiveDiag(diagnostic, tree, "charge_density"))
        infoDS(ions.chargeDensity(), "charge_density", patchAttributes[lvlPatchID]["ion"]);
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
    auto& h5file   = Super::h5FileForQuantity(diagnostic);

    auto const writeGhosts = [&](auto& path, auto& attr, std::string key, auto null) {
        this->writeGhostsAttr_(h5file, path,
                               null ? 0 : attr[key + "_ghosts"].template to<std::size_t>(), null);
    };

    auto const initDS = [&](auto& path, auto& attr, std::string key, auto null) {
        auto dsPath = path + key;
        h5Writer.template createDataSet<FloatType>(
            h5file, dsPath,
            null ? std::vector<std::size_t>(GridLayout::dimension, 0)
                 : attr[key].template to<std::vector<std::size_t>>());
        writeGhosts(dsPath, attr, key, null);
    };
    auto const initVF = [&](auto& path, auto& attr, std::string key, auto null) {
        for (auto& [id, type] : core::Components::componentMap())
            initDS(path, attr, key + "_" + id, null);
    };
    auto const initTF = [&](auto& path, auto& attr, std::string key, auto null) {
        for (auto& [id, type] : core::Components::componentMap<2>())
            initDS(path, attr, key + "_" + id, null);
    };

    auto const initPatch = [&](auto& lvl, auto& attr, std::string patchID = "") {
        bool null        = patchID.empty();
        std::string path = h5Writer.getPatchPathAddTimestamp(lvl, patchID) + "/";

        for (auto& pop : ions)
        {
            std::string popId{"fluid_" + pop.name()};
            std::string tree{"/ions/pop/" + pop.name() + "/"};
            std::string popPath(path + "pop/" + pop.name() + "/");
            if (isActiveDiag(diagnostic, tree, "density"))
                initDS(path, attr[popId], "density", null);
            if (isActiveDiag(diagnostic, tree, "charge_density"))
                initDS(path, attr[popId], "charge_density", null);
            if (isActiveDiag(diagnostic, tree, "flux"))
                initVF(path, attr[popId], "flux", null);
            if (isActiveDiag(diagnostic, tree, "momentum_tensor"))
                initTF(path, attr[popId], "momentum_tensor", null);
            if (isActiveDiag(diagnostic, tree, "kinetic_energy_flux_vector"))
                initVF(path, attr[popId], "kinetic_energy_flux_vector", null);
        }

        std::string tree{"/ions/"};
        if (isActiveDiag(diagnostic, tree, "charge_density"))
            initDS(path, attr["ion"], "charge_density", null);
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
    auto& h5file   = Super::h5FileForQuantity(diagnostic);

    auto writeDS = [&](auto path, auto& field) {
        h5file.template write_data_set_flat<GridLayout::dimension>(path, field.data());
    };
    auto writeTF
        = [&](auto path, auto& vecF) { h5Writer.writeTensorFieldAsDataset(h5file, path, vecF); };

    std::string path = h5Writer.patchPath() + "/";
    for (auto& pop : ions)
    {
        std::string tree{"/ions/pop/" + pop.name() + "/"};
        if (isActiveDiag(diagnostic, tree, "density"))
            writeDS(path + "density", pop.particleDensity());
        if (isActiveDiag(diagnostic, tree, "charge_density"))
            writeDS(path + "charge_density", pop.chargeDensity());
        if (isActiveDiag(diagnostic, tree, "flux"))
            writeTF(path + "flux", pop.flux());
        if (isActiveDiag(diagnostic, tree, "momentum_tensor"))
            writeTF(path + "momentum_tensor", pop.momentumTensor());
        if (isActiveDiag(diagnostic, tree, "kinetic_energy_flux_vector"))
            writeTF(path + "kinetic_energy_flux_vector", pop.kineticEnergyFlux());
    }

    std::string tree{"/ions/"};
    if (isActiveDiag(diagnostic, tree, "charge_density"))
        writeDS(path + "charge_density", ions.chargeDensity());
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
        checkWrite(tree, "kinetic_energy_flux_vector", pop);
    }

    writeAttributes_(diagnostic, h5file, fileAttributes, patchAttributes, maxLevel);
}

} // namespace PHARE::diagnostic::h5

#endif /* PHARE_DIAGNOSTIC_DETAIL_TYPES_FLUID_H */
