
#ifndef PHARE_DIAGNOSTIC_MANAGER_HPP_
#define PHARE_DIAGNOSTIC_MANAGER_HPP_

#include "core/data/particles/particle_array.h"
#include "initializer/data_provider.h"
#include "diagnostic_props.h"

#include <utility>

namespace PHARE::diagnostic
{
enum class Mode { LIGHT, FULL };



template<typename DiagManager>
void registerDiagnostics(DiagManager& dMan, PHARE::initializer::PHAREDict& diagsParams)
{
    std::vector<std::string> const diagTypes = {"fluid", "electromag", "particle"};

    for (auto& diagType : diagTypes)
    {
        size_t diagBlockID = 0;
        while (diagsParams.contains(diagType)
               && diagsParams[diagType].contains(diagType + std::to_string(diagBlockID)))
        {
            const std::string diagName = diagType + std::to_string(diagBlockID);
            dMan.addDiagDict(diagsParams[diagType][diagName]);
            diagBlockID++;
        }
    }
}




class IDiagnosticsManager
{
public:
    virtual void dump(double timeStamp, double timeStep) = 0;
    inline virtual ~IDiagnosticsManager();
};
IDiagnosticsManager::~IDiagnosticsManager() {}



template<typename Writer>
class DiagnosticsManager : public IDiagnosticsManager
{
public:
    DiagnosticsManager(std::unique_ptr<Writer>&& writer_ptr)
        : writer_{std::forward<std::unique_ptr<Writer>>(writer_ptr)}
    {
        if (!writer_)
            throw std::runtime_error("Error: DiagnosticsManager received null Writer");
    }

    template<typename Hierarchy, typename Model>
    static std::unique_ptr<DiagnosticsManager> make_unique(Hierarchy& hier, Model& model,
                                                           initializer::PHAREDict& dict)
    {
        auto dMan = std::make_unique<DiagnosticsManager>(Writer::make_unique(hier, model, dict));
        registerDiagnostics(*dMan, dict);
        return dMan;
    }

    void dump(double timeStamp, double timeStep) override;

    DiagnosticsManager& addDiagDict(PHARE::initializer::PHAREDict& dict);


    DiagnosticsManager& addDiagDict(PHARE::initializer::PHAREDict&& dict)
    {
        return addDiagDict(dict);
    }


    void addDiagnostic(DiagnosticProperties& diagnostic)
    {
        diagnostics_.emplace_back(diagnostic);
        nextWrite_[diagnostic.type + diagnostic.quantity]   = -1;
        nextCompute_[diagnostic.type + diagnostic.quantity] = -1;
    }

    auto& diagnostics() const { return diagnostics_; }


    bool needsWrite(DiagnosticProperties& diag, double timeStamp, double timeStep)
    {
        auto nextWrite = nextWrite_[diag.type + diag.quantity];
        return timeStamp + timeStep > diag.writeTimestamps[nextWrite];
    }

    bool needsCompute(DiagnosticProperties& diag, double timeStamp, double timeStep)
    {
        auto nextCompute = nextCompute_[diag.type + diag.quantity];
        return timeStamp + timeStep > diag.computeTimestamps[nextCompute];
    }

    Writer& writer() { return *writer_.get(); }

protected:
    std::vector<DiagnosticProperties> diagnostics_;

private:
    std::unique_ptr<Writer> writer_;

    std::map<std::string, std::size_t> nextCompute_;
    std::map<std::string, std::size_t> nextWrite_;


    DiagnosticsManager(const DiagnosticsManager&)  = delete;
    DiagnosticsManager(const DiagnosticsManager&&) = delete;
    DiagnosticsManager& operator=(const DiagnosticsManager&) = delete;
    DiagnosticsManager& operator=(const DiagnosticsManager&&) = delete;
};




template<typename Writer>
DiagnosticsManager<Writer>&
DiagnosticsManager<Writer>::addDiagDict(PHARE::initializer::PHAREDict& diagInputs)
{
    auto& diagProps           = diagnostics_.emplace_back(DiagnosticProperties{});
    diagProps.type            = diagInputs["type"].template to<std::string>();
    diagProps.quantity        = diagInputs["quantity"].template to<std::string>();
    diagProps.writeTimestamps = diagInputs["write_timestamps"].template to<std::vector<double>>();
    diagProps.computeTimestamps
        = diagInputs["compute_timestamps"].template to<std::vector<double>>();

    return *this;
}




template<typename Writer>
void DiagnosticsManager<Writer>::dump(double timeStamp, double timeStep)
{
    std::vector<DiagnosticProperties*> activeDiagnostics;
    for (auto& diag : diagnostics_)
    {
        auto diagID = diag.type + diag.quantity;

        if (needsCompute(diag, timeStamp, timeStep))
        {
            writer_->getDiagnosticWriterForType(diag.type)->compute(diag);
            nextCompute_[diagID]++;
        }
        if (needsWrite(diag, timeStamp, timeStep))
        {
            activeDiagnostics.emplace_back(&diag);
        }
    }
    writer_->dump(activeDiagnostics, timeStamp);

    for (auto const* diag : activeDiagnostics)
    {
        nextWrite_[diag->type + diag->quantity]++;
    }
}

} // namespace PHARE::diagnostic

#endif /* PHARE_DIAGNOSTIC_MANAGER_HPP_ */
