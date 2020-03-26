
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
            std::string diagName = diagType + std::to_string(diagBlockID);
            auto copy            = diagsParams[diagType][diagName];
            copy["category"]     = diagType;
            dMan.addDiagDict(copy);
            diagBlockID++;
        }
    }
}




class IDiagnosticsManager
{
public:
    virtual void dump(double timeStamp, double timeStep) = 0;
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


    DiagnosticsManager(const DiagnosticsManager&)  = delete;
    DiagnosticsManager(const DiagnosticsManager&&) = delete;
    DiagnosticsManager& operator=(const DiagnosticsManager&) = delete;
    DiagnosticsManager& operator=(const DiagnosticsManager&&) = delete;


    static std::unique_ptr<DiagnosticsManager> from(Writer& writer, initializer::PHAREDict& dict)
    {
        auto dMan = std::make_unique<DiagnosticsManager>(writer);
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
        lastWrite_[diagnostic.type + diagnostic.quantity]   = 0;
        lastCompute_[diagnostic.type + diagnostic.quantity] = 0;
    }

    auto& diagnostics() const { return diagnostics_; }


    bool needsWrite(DiagnosticProperties& diag, double timeStamp, double timeStep)
    {
        auto lastWrite = lastWrite_[diag.type + diag.quantity];
        return diag.writeTimestamps[lastWrite] + timeStep > timeStamp;
    }

    bool needsCompute(DiagnosticProperties& diag, double timeStamp, double timeStep)
    {
        auto lastCompute = lastCompute_[diag.type + diag.quantity];
        return diag.computeTimestamps[lastCompute] + timeStep > timeStamp;
    }



protected:
    std::vector<DiagnosticProperties> diagnostics_;

private:
    Writer& writer_;


    std::map<std::string, std::size_t> lastCompute_;
    std::map<std::string, std::size_t> lastWrite_;
};




template<typename Writer>
DiagnosticsManager<Writer>&
DiagnosticsManager<Writer>::addDiagDict(PHARE::initializer::PHAREDict& diagInputs)
{
    auto& diagProps           = diagnostics_.emplace_back(DiagnosticProperties{});
    diagProps.type            = diagInputs["category"].template to<std::string>();
    diagProps.quantity        = diagInputs["type"].template to<std::string>();
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
            writer_.getDiagnosticWriterForType(diag.type)->compute(diag);
            lastCompute_[diagID]++;
        }
        if (needsWrite(diag, timeStamp, timeStep))
        {
            activeDiagnostics.emplace_back(&diag);
        }
    }
    writer_.dump(activeDiagnostics, timeStamp);

    for (auto const* diag : activeDiagnostics)
    {
        lastWrite_[diag->type + diag->quantity]++;
    }
}




} // namespace PHARE::diagnostic

#endif /* PHARE_DIAGNOSTIC_MANAGER_HPP_ */
