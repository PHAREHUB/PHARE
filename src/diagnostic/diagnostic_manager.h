
#ifndef PHARE_DIAGNOSTIC_MANAGER_HPP_
#define PHARE_DIAGNOSTIC_MANAGER_HPP_

#include "core/data/particles/particle_array.h"
#include "initializer/data_provider.h"
#include "diagnostic_dao.h"

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


    void addDiagnostic(DiagnosticProperties& diagnostic) { diagnostics_.emplace_back(diagnostic); }

    auto& diagnostics() const { return diagnostics_; }


protected:
    std::vector<DiagnosticProperties> diagnostics_;

private:
    Writer& writer_;

    DiagnosticsManager(const DiagnosticsManager&)             = delete;
    DiagnosticsManager(const DiagnosticsManager&&)            = delete;
    DiagnosticsManager& operator&(const DiagnosticsManager&)  = delete;
    DiagnosticsManager& operator&(const DiagnosticsManager&&) = delete;
};




template<typename Writer>
DiagnosticsManager<Writer>&
DiagnosticsManager<Writer>::addDiagDict(PHARE::initializer::PHAREDict& diagInputs)
{
    auto& diagProps           = diagnostics_.emplace_back(DiagnosticProperties{});
    diagProps.type            = diagInputs["category"].template to<std::string>();
    diagProps.start_iteration = diagInputs["start_iteration"].template to<double>();
    diagProps.last_iteration  = diagInputs["last_iteration"].template to<double>();
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
        if (diag.needsCompute(timeStamp, timeStep))
        {
            writer_.getDiagnosticWriterForType(diag.type)->compute(diag);
            diag.lastCompute++;
        }
        if (diag.needsWrite(timeStamp, timeStep))
        {
            activeDiagnostics.emplace_back(&diag);
        }
    }
    writer_.dump(activeDiagnostics, timeStamp);
}




} // namespace PHARE::diagnostic

#endif /* PHARE_DIAGNOSTIC_MANAGER_HPP_ */
