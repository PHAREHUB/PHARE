
#ifndef PHARE_DIAGNOSTIC_MANAGER_HPP_
#define PHARE_DIAGNOSTIC_MANAGER_HPP_

#include "core/data/particles/particle_array.h"
#include "initializer/data_provider.h"
#include "solver/physical_models/hybrid_model.h"
#include "diagnostic_dao.h"

#include <utility>

namespace PHARE::diagnostic
{
enum class Mode { LIGHT, FULL };


bool isActive(std::vector<double>& active, double time_stamp, double time_step, size_t& last)
{
    bool is_active = active[last] + time_step > time_stamp;

    if (is_active)
        last++;

    return is_active;
}

template<typename Diag>
bool needsCompute(Diag& diag, double time_stamp, double time_step)
{
    return isActive(diag.computeTimestamps, time_stamp, time_step, diag.lastCompute);
}

template<typename Diag>
bool needsWrite(Diag& diag, double time_stamp, double time_step)
{
    return isActive(diag.writeTimestamps, time_stamp, time_step, diag.lastWrite);
}



template<typename DiagManager>
void handleInputDiagnostics(DiagManager& dMan, PHARE::initializer::PHAREDict& diags)
{
    std::vector<std::string> const keys = {"fluid", "electromag", "particle"};

    for (auto& key : keys)
    {
        size_t it = 0;
        while (diags.contains(key) && diags[key].contains(key + std::to_string(it)))
        {
            std::string path = key + std::to_string(it);
            auto copy        = diags[key][path];
            copy["category"] = key;
            dMan.addDiagDict(copy);
            it++;
        }
    }
}

class IDiagnosticsManager
{
public:
    virtual void dump(double time_stamp, double time_step) = 0;
    virtual ~IDiagnosticsManager();
};
IDiagnosticsManager::~IDiagnosticsManager() {}

template<typename Writer>
class DiagnosticsManager : public IDiagnosticsManager
{
public:
    DiagnosticsManager(Writer& writer)
        : writer_{writer}
    {
    }


    static std::unique_ptr<DiagnosticsManager> from(Writer& writer, initializer::PHAREDict& dict)
    {
        auto dMan = std::make_unique<DiagnosticsManager>(writer);
        handleInputDiagnostics(*dMan, dict);
        return dMan;
    }


    void dump(double time_stamp, double time_step) override;
    DiagnosticsManager& addDiagDict(PHARE::initializer::PHAREDict& dict);
    DiagnosticsManager& addDiagDict(PHARE::initializer::PHAREDict&& dict)
    {
        return addDiagDict(dict);
    }
    void addDiagnostic(DiagnosticDAO& diagnostic) { diagnostics_.emplace_back(diagnostic); }

    auto& diagnostics() const { return diagnostics_; }

protected:
    std::vector<DiagnosticDAO> diagnostics_;

private:
    Writer& writer_;

    DiagnosticsManager(const DiagnosticsManager&)             = delete;
    DiagnosticsManager(const DiagnosticsManager&&)            = delete;
    DiagnosticsManager& operator&(const DiagnosticsManager&)  = delete;
    DiagnosticsManager& operator&(const DiagnosticsManager&&) = delete;
};

template<typename Writer>
DiagnosticsManager<Writer>&
DiagnosticsManager<Writer>::addDiagDict(PHARE::initializer::PHAREDict& dict)
{
    auto& dao    = diagnostics_.emplace_back(DiagnosticDAO{});
    dao.category = dict["category"].template to<std::string>();
    dao.type     = dict["type"].template to<std::string>();

    dao.writeTimestamps   = dict["write_timestamps"].template to<std::vector<double>>();
    dao.computeTimestamps = dict["compute_timestamps"].template to<std::vector<double>>();

    return *this;
}



/*TODO
   iterations
*/
template<typename Writer>
void DiagnosticsManager<Writer>::dump(double time_stamp, double time_step)
{
    std::vector<DiagnosticDAO*> activeDiagnostics;
    for (auto& diag : diagnostics_)
    {
        if (needsCompute(diag, time_stamp, time_step))
        {
            writer_.getDiagnosticWriterForType(diag.category)->compute(diag);
        }
        if (needsWrite(diag, time_stamp, time_step))
        {
            activeDiagnostics.emplace_back(&diag);
        }
    }
    writer_.dump(activeDiagnostics, time_stamp);
}



// Generic Template declaration, to override per Concrete model type
class IDiagnosticModelView
{
public:
    virtual ~IDiagnosticModelView();
};
IDiagnosticModelView::~IDiagnosticModelView() {}

template<typename Model, typename ModelParams>
class DiagnosticModelView : public IDiagnosticModelView
{
};


// HybridModel<Args...> specialization
template<typename ModelParams>
class DiagnosticModelView<solver::type_list_to_hybrid_model_t<ModelParams>, ModelParams>
    : public IDiagnosticModelView
{
public:
    using Model      = solver::type_list_to_hybrid_model_t<ModelParams>;
    using VecField   = typename Model::vecfield_type;
    using GridLayout = typename Model::gridLayout_type;
    using Attributes = cppdict::Dict<float, double, size_t, std::string>;

    static constexpr auto dimension = Model::dimension;

    DiagnosticModelView(Model& model)
        : model_{model}
    {
    }

    std::vector<VecField*> getElectromagFields() const
    {
        return {&model_.state.electromag.B, &model_.state.electromag.E};
    }

    auto& getIons() const { return model_.state.ions; }

    auto getPatchAttributes(GridLayout& grid);
    static auto getEmptyPatchAttributes();


protected:
    Model& model_;
};



template<typename ModelParams>
auto DiagnosticModelView<solver::type_list_to_hybrid_model_t<ModelParams>,
                         ModelParams>::getPatchAttributes(GridLayout& grid)
{
    Attributes dict;
    dict["origin"]   = grid.origin().str();
    dict["nbrCells"] = core::Point<uint32_t, dimension>{grid.nbrCells()}.str();
    return dict;
}

// this function is a for std visit and MPI sychronization
// so we can write attributes for disparate numbers of patches across MPI nodes
// This function should match getPatchAttributes in structure
template<typename ModelParams>
auto DiagnosticModelView<solver::type_list_to_hybrid_model_t<ModelParams>,
                         ModelParams>::getEmptyPatchAttributes()
{
    Attributes dict;
    dict["origin"]   = std::string{""};
    dict["nbrCells"] = std::string{""};
    return dict;
}



template<size_t dim>
class ParticlePacker
{
public:
    ParticlePacker(core::ParticleArray<dim> const& particles)
        : particles_{particles}
    {
    }

    static auto get(core::Particle<dim> const& particle)
    {
        return std::forward_as_tuple(particle.weight, particle.charge, particle.iCell,
                                     particle.delta, particle.v);
    }

    static auto empty()
    {
        core::Particle<dim> particle;
        return get(particle);
    }

    static auto& keys() { return keys_; }

    auto get(size_t i) const { return get(particles_[i]); }
    bool hasNext() const { return it_ < particles_.size(); }
    auto next() { return get(it_++); }

private:
    core::ParticleArray<dim> const& particles_;
    size_t it_ = 0;
    static inline std::array<std::string, 5> keys_{"weight", "charge", "iCell", "delta", "v"};
};



template<std::size_t dim>
struct ContiguousParticles
{
    std::vector<int> iCell;
    std::vector<float> delta;
    std::vector<double> weight, charge, v;
    size_t size_;
    auto size() const { return size_; }
    ContiguousParticles(size_t s)
        : iCell(s * dim)
        , delta(s * dim)
        , weight(s)
        , charge(s)
        , v(s * 3)
        , size_(s)
    {
    }
};



// generic subclass of model specialized superclass
template<typename Hierarchy, typename Model>
class AMRDiagnosticModelView : public DiagnosticModelView<Model, typename Model::type_list>
{
public:
    using Super      = DiagnosticModelView<Model, typename Model::type_list>;
    using ResMan     = typename Model::resources_manager_type;
    using GridLayout = typename Model::gridLayout_type;
    using Super::model_;
    static constexpr auto dimension = Model::dimension;

    AMRDiagnosticModelView(Hierarchy& hierarchy, Model& model)
        : Super{model}
        , hierarchy_{hierarchy}
    {
    }


    template<typename Action, typename... Args>
    void visitHierarchy(Action&& action, int minLevel = 0, int maxLevel = 0)
    {
        auto& resMan = *model_.resourcesManager;
        PHARE::amr::visitHierarchy<GridLayout>(hierarchy_, resMan, std::forward<Action>(action),
                                               minLevel, maxLevel, model_);
    }


    std::string getLayoutTypeString() { return std::string{GridLayout::implT::type}; }


private:
    Hierarchy& hierarchy_;

    AMRDiagnosticModelView(const AMRDiagnosticModelView&)             = delete;
    AMRDiagnosticModelView(const AMRDiagnosticModelView&&)            = delete;
    AMRDiagnosticModelView& operator&(const AMRDiagnosticModelView&)  = delete;
    AMRDiagnosticModelView& operator&(const AMRDiagnosticModelView&&) = delete;
};


} // namespace PHARE::diagnostic

#endif /* PHARE_DIAGNOSTIC_MANAGER_HPP_ */
