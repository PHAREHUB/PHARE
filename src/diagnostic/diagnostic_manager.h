
#ifndef PHARE_Diagnostic_MANAGER_HPP_
#define PHARE_Diagnostic_MANAGER_HPP_

#include "data_provider.h"
#include "utilities/types.h"
#include "physical_models/hybrid_model.h"
#include "resources_manager/resources_manager.h"

namespace PHARE
{
namespace diagnostic
{
    enum class Mode { LIGHT, FULL };

    struct Types
    {
        static std::string_view constexpr rho_s = "rho_s", flux_s = "flux_s", E = "E", B = "B",
                                          space_box = "space_box";
    };
    struct FieldInfo
    {
        double const* const data; // B/E xyz
        size_t const size;        // B/E xyz
        std::string const id;
    };
} /*namespace diagnostic*/

struct Diagnostic
{
    size_t compute_every = 1, write_every = 1;
    size_t start_iteration = 0, end_iteration = -1; /* likely to be time rather than index*/
                                                    /* do we allow ranges?*/
    std::string name, species, type;
};

class DiagnosticWriter
{
public:
    virtual void write(Diagnostic&)   = 0;
    virtual void compute(Diagnostic&) = 0;
    virtual ~DiagnosticWriter() {}
};

template<typename Writer>
class DiagnosticsManager
{
public:
    using DiagnosticWritingList = std::vector<
        std::tuple<std::reference_wrapper<Diagnostic>, std::shared_ptr<DiagnosticWriter>>>;
    DiagnosticsManager(Writer& writer)
        : writer_(writer)
    {
    }

    DiagnosticsManager& addDiagDict(PHARE::initializer::PHAREDict<1>& dict);
    DiagnosticsManager& addDiagDict(PHARE::initializer::PHAREDict<1>&& dict)
    {
        return addDiagDict(dict);
    }

    void dump();

private:
    Writer& writer_;
    std::vector<Diagnostic> diagnostics;

    DiagnosticsManager(const DiagnosticsManager&)             = delete;
    DiagnosticsManager(const DiagnosticsManager&&)            = delete;
    DiagnosticsManager& operator&(const DiagnosticsManager&)  = delete;
    DiagnosticsManager& operator&(const DiagnosticsManager&&) = delete;
};

/*TODO
   iterations
*/
template<typename Writer>
void DiagnosticsManager<Writer>::dump(/*time iteration*/)
{
    size_t iter = 1; // TODO replace with time/iteration idx
    auto active = [](size_t current, size_t start, size_t end, size_t every) {
        return (current >= start && current <= end) && current % every == 0;
    };
    DiagnosticWritingList diagnosticWriters;
    for (auto& diag : diagnostics)
    {
        if (active(iter, diag.start_iteration, diag.end_iteration, diag.compute_every))
        {
            writer_.getWriter(diag.type)->compute(diag);
        }
        if (active(iter, diag.start_iteration, diag.end_iteration, diag.write_every))
        {
            diagnosticWriters.emplace_back(diag, writer_.getWriter(diag.type));
        }
    }
    writer_.dump(diagnosticWriters);
}

template<typename Writer>
DiagnosticsManager<Writer>&
    DiagnosticsManager<Writer>::addDiagDict(PHARE::initializer::PHAREDict<1>& dict)
{
    size_t &compute_every   = dict["diag"]["compute_every"].template to<std::size_t>(),
           &write_every     = dict["diag"]["write_every"].template to<std::size_t>(),
           &start_iteration = dict["diag"]["start_iteration"].template to<std::size_t>(),
           &end_iteration   = dict["diag"]["end_iteration"].template to<std::size_t>();
    std::string &name       = dict["diag"]["name"].template to<std::string>(),
                &species    = dict["diag"]["species"].template to<std::string>(),
                &type       = dict["diag"]["type"].template to<std::string>();

    diagnostics.emplace_back(PHARE::core::aggregate_adapter<Diagnostic>(
        compute_every, write_every, start_iteration, end_iteration, name, species, type));
    return *this;
}

// Generic Template declaration, to override per Concrete model type
template<typename Model, typename ModelParams>
class DiagnosticModelView
{
public:
    using Fields = std::vector<std::shared_ptr<diagnostic::FieldInfo>>;

    std::vector<Fields> getElectromagFields() = 0;
    std::tuple<void> getResources()           = 0;
};

// HybridModel<Args...> specialization
template<typename ModelParams>
class DiagnosticModelView<solver::type_list_to_hybrid_model_t<ModelParams>, ModelParams>
{
public:
    using Model      = solver::type_list_to_hybrid_model_t<ModelParams>;
    using ResMan     = typename Model::resources_manager_type;
    using VecField   = typename Model::vecfield_type;
    using Electromag = typename Model::electromag_type;
    using Ions       = typename Model::ions_type;
    using Grid       = typename Model::gridLayout_type;
    using Guard      = amr::ResourcesGuard<ResMan, Electromag, Ions>; // inherited
    using Resources  = std::tuple<Electromag&, Ions&>;
    using Fields     = std::vector<std::shared_ptr<diagnostic::FieldInfo>>;
    using Attributes = cppdict::Dict<float, double, size_t, std::string>;

    static constexpr auto dimensions = Model::dimension;

    using ParticleInfo = cppdict::Dict<std::array<int, dimensions>, std::array<float, dimensions>,
                                       std::array<double, 3>, size_t, double>;

    DiagnosticModelView(Model& model)
        : model_(model)
    {
    }

    std::vector<Fields> getElectromagFields() const { return {getB(), getE()}; }

    auto& getParticles() const { return model_.state.ions; };

    auto getParticlePacker(std::vector<core::Particle<1>> const&);

    auto getPatchAttributes(Grid& grid);

    Resources getResources()
    {
        return std::forward_as_tuple(model_.state.electromag, model_.state.ions);
    }


protected:
    Fields getB() const { return get(model_.state.electromag.B, "B"); }

    Fields getE() const { return get(model_.state.electromag.E, "E"); }

    auto get(VecField& vecField, std::string&& id) const
    {
        std::vector<std::shared_ptr<diagnostic::FieldInfo>> fInfo;
        for (auto const& key : {"x", "y", "z"})
        {
            auto& field = vecField.getComponent(core::Components::at(key));
            fInfo.emplace_back(std::make_shared<core::aggregate_adapter<diagnostic::FieldInfo>>(
                field.data(), field.size(), id + key));
        }
        return fInfo;
    }

    std::string getPatchOrigin(Grid& grid)
    {
        auto& bieldx = model_.state.electromag.B.getComponent(core::Component::X);
        return grid
            .fieldNodeCoordinates(bieldx, grid.origin(),
                                  grid.ghostStartIndex(bieldx, core::Direction::X))
            .str();
    }

protected:
    Model& model_;
};

template<typename ModelParams>
auto DiagnosticModelView<solver::type_list_to_hybrid_model_t<ModelParams>,
                         ModelParams>::getPatchAttributes(Grid& grid)
{
    Attributes dict;
    dict["id"]     = std::string("id");
    dict["float"]  = 0.0f;
    dict["double"] = 0.0;
    dict["size_t"] = size_t{0};
    dict["origin"] = getPatchOrigin(grid);
    return dict;
}

class ParticularParticlePartPackerPicker
{
public:
    ParticularParticlePartPackerPicker(std::vector<core::Particle<1>> const& particles)
        : particles_{particles}
        , keys_{{"weight", "charge", "iCell", "delta", "v"}}
    {
    }

    auto get(size_t i) const
    {
        auto& particle = particles_[i];
        return std::forward_as_tuple(particle.weight, particle.charge, particle.iCell,
                                     particle.delta, particle.v);
    }

    auto& keys() const { return keys_; }
    bool hasNext() const { return it_ < particles_.size(); }
    auto next() { return get(it_++); }

private:
    size_t it_ = 0;
    std::vector<core::Particle<1>> const& particles_;
    std::vector<std::string> keys_;
};

template<typename PartPicker>
class ParticularPacker
{
public:
    ParticularPacker(std::vector<core::Particle<1>> const& particles)
        : picker_{particles}
    {
    }

    auto& keys() const { return picker_.keys(); }
    bool hasNext() const { return picker_.hasNext(); }
    auto next() { return picker_.next(); }

    auto first() const { return picker_.get(0); }
private:
    PartPicker picker_;
};

template<typename ModelParams>
auto DiagnosticModelView<solver::type_list_to_hybrid_model_t<ModelParams>, ModelParams>::
    getParticlePacker(std::vector<core::Particle<1>> const& particles)
{
    return ParticularPacker<ParticularParticlePartPackerPicker>{particles};
}

} // namespace PHARE

#endif //  PHARE_Diagnostic_MANAGER_HPP_
