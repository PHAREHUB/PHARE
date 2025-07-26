#ifndef PHARE_DIAGNOSTIC_MANAGER_HPP_
#define PHARE_DIAGNOSTIC_MANAGER_HPP_

#include "core/def.hpp"
#include "core/data/particles/particle_array.hpp"
#include "diagnostic/diagnostic_model_view.hpp"
#include "initializer/data_provider.hpp"
#include "diagnostic_props.hpp"

#include <type_traits>
#include <utility>
#include <cmath>
#include <memory>
#include <map>

namespace PHARE::diagnostic
{
enum class Mode { LIGHT, FULL };



template<typename DiagManager>
void registerDiagnostics(DiagManager& dMan, initializer::PHAREDict const& diagsParams)
{
    auto const diagTypes = []() {
        if constexpr (solver::is_hybrid_model_v<typename DiagManager::Model_t>)
        {
            return std::vector<std::string>{"fluid", "electromag", "particles", "meta", "info"};
        }
        else if constexpr (solver::is_mhd_model_v<typename DiagManager::Model_t>)
        {
            return std::vector<std::string>{"mhd", "meta", "electromag"};
        }
        else
        {
            static_assert(false, "Unsupported model type");
            return std::vector<std::string>{};
        }
    }();

    for (auto& diagType : diagTypes)
    {
        // several diags of the same type can be registeroed
        // corresponding to several blocks in user input
        // fluid0, fluid1, electromag0, electromag1, etc.
        std::size_t diagBlockID = 0;
        while (diagsParams.contains(diagType)
               && diagsParams[diagType].contains(diagType + std::to_string(diagBlockID)))
        {
            std::string const diagName = diagType + std::to_string(diagBlockID);
            dMan.addDiagDict(diagsParams[diagType][diagName]);
            ++diagBlockID;
        }
    }
}




class IDiagnosticsManager
{
public:
    virtual bool dump(double timeStamp, double timeStep)         = 0;
    virtual void dump_level(std::size_t level, double timeStamp) = 0;
    inline virtual ~IDiagnosticsManager();
};
IDiagnosticsManager::~IDiagnosticsManager() {}



template<typename Writer>
class DiagnosticsManager : public IDiagnosticsManager
{
public:
    using Model_t = typename Writer::Model_t;

    bool dump(double timeStamp, double timeStep) override;


    void dump_level(std::size_t level, double timeStamp) override;


    DiagnosticsManager(std::unique_ptr<Writer>&& writer_ptr)
        : writer_{std::move(writer_ptr)}
    {
        if (!writer_)
            throw std::runtime_error("Error: DiagnosticsManager received null Writer");
    }


    template<typename Hierarchy, typename Model>
    NO_DISCARD static std::unique_ptr<DiagnosticsManager>
    make_unique(Hierarchy& hier, Model& model, initializer::PHAREDict const& dict)
    {
        auto dMan = std::make_unique<DiagnosticsManager>(Writer::make_unique(hier, model, dict));
        registerDiagnostics(*dMan, dict);
        return dMan;
    }


    DiagnosticsManager& addDiagDict(initializer::PHAREDict const& dict);


    DiagnosticsManager& addDiagDict(initializer::PHAREDict&& dict) { return addDiagDict(dict); }


    NO_DISCARD auto& diagnostics() const { return diagnostics_; }


    NO_DISCARD Writer& writer() { return *writer_.get(); }


    DiagnosticsManager(DiagnosticsManager const&)            = delete;
    DiagnosticsManager(DiagnosticsManager&&)                 = delete;
    DiagnosticsManager& operator=(DiagnosticsManager const&) = delete;
    DiagnosticsManager& operator=(DiagnosticsManager&&)      = delete;

private:
    NO_DISCARD bool needsAction_(double nextTime, double timeStamp, double timeStep)
    {
        // casting to float to truncate double to avoid trailing imprecision
        return static_cast<float>(std::abs(nextTime - timeStamp)) < static_cast<float>(timeStep);
    }


    NO_DISCARD bool needsWrite_(DiagnosticProperties& diag, double timeStamp, double timeStep)
    {
        auto nextWrite = nextWrite_[diag.type + diag.quantity];
        return nextWrite < diag.writeTimestamps.size()
               and needsAction_(diag.writeTimestamps[nextWrite], timeStamp, timeStep);
    }


    NO_DISCARD bool needsCompute_(DiagnosticProperties& diag, double timeStamp, double timeStep)
    {
        auto nextCompute = nextCompute_[diag.type + diag.quantity];
        return nextCompute < diag.computeTimestamps.size()
               and needsAction_(diag.computeTimestamps[nextCompute], timeStamp, timeStep);
    }


    std::vector<DiagnosticProperties> diagnostics_;
    std::unique_ptr<Writer> writer_;
    std::map<std::string, std::size_t> nextCompute_;
    std::map<std::string, std::size_t> nextWrite_;
};




template<typename Writer>
DiagnosticsManager<Writer>&
DiagnosticsManager<Writer>::addDiagDict(initializer::PHAREDict const& diagParams)
{
    auto& diagProps           = diagnostics_.emplace_back(DiagnosticProperties{});
    diagProps.type            = diagParams["type"].template to<std::string>();
    diagProps.quantity        = diagParams["quantity"].template to<std::string>();
    diagProps.writeTimestamps = diagParams["write_timestamps"].template to<std::vector<double>>();
    diagProps["flush_every"]  = diagParams["flush_every"].template to<std::size_t>();

    diagProps.computeTimestamps
        = diagParams["compute_timestamps"].template to<std::vector<double>>();

    diagProps.nAttributes = diagParams["n_attributes"].template to<std::size_t>();
    for (std::size_t i = 0; i < diagProps.nAttributes; ++i)
    {
        std::string idx = std::to_string(i);
        std::string key = diagParams["attribute_" + idx + "_key"].template to<std::string>();
        if (key == "heat_capacity_ratio")
        {
            double val = diagParams["attribute_" + idx + "_value"].template to<double>();
            diagProps.fileAttributes[key] = val;
        }
        else
        {
            std::string val = diagParams["attribute_" + idx + "_value"].template to<std::string>();
            diagProps.fileAttributes[key] = val;
        }
    }

    return *this;
}


template<typename Writer>
void DiagnosticsManager<Writer>::dump_level(std::size_t level, double timeStamp)
{
    std::vector<DiagnosticProperties*> activeDiagnostics;

    for (auto& diag : diagnostics_)
        activeDiagnostics.emplace_back(&diag);

    writer_->dump_level(level, activeDiagnostics, timeStamp);
}


template<typename Writer>
bool DiagnosticsManager<Writer>::dump(double timeStamp, double timeStep)
{
    std::vector<DiagnosticProperties*> activeDiagnostics;
    for (auto& diag : diagnostics_)
    {
        auto diagID = diag.type + diag.quantity;

        if (needsCompute_(diag, timeStamp, timeStep))
        {
            writer_->getDiagnosticWriterForType(diag.type)->compute(diag);
            nextCompute_[diagID]++;
        }
        if (needsWrite_(diag, timeStamp, timeStep))
        {
            activeDiagnostics.emplace_back(&diag);
        }
    }

    if (activeDiagnostics.size() > 0)
    {
        PHARE_LOG_SCOPE(1, "DiagnosticsManager::dump");
        writer_->dump(activeDiagnostics, timeStamp);

        for (auto const* diag : activeDiagnostics)
        {
            nextWrite_[diag->type + diag->quantity]++;
        }
    }

    return activeDiagnostics.size() > 0;
}

} // namespace PHARE::diagnostic

#endif /* PHARE_DIAGNOSTIC_MANAGER_HPP_ */
