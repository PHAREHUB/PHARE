#ifndef PHARE_DIAGNOSTIC_HI5_DIAGNOSTIC_HPP_
#define PHARE_DIAGNOSTIC_HI5_DIAGNOSTIC_HPP_

#include "diagnosticmanager.h"

namespace PHARE
{
class Hi5DiagnosticsManager : public DiagnosticsManager
{
public:
    Hi5DiagnosticsManager() {}

    void init() override;
    void log(Diagnostic&) override;

protected:
    virtual ~Hi5DiagnosticsManager() {}

private:
    Hi5DiagnosticsManager(const Hi5DiagnosticsManager&)             = delete;
    Hi5DiagnosticsManager(const Hi5DiagnosticsManager&&)            = delete;
    Hi5DiagnosticsManager& operator&(const Hi5DiagnosticsManager&)  = delete;
    Hi5DiagnosticsManager& operator&(const Hi5DiagnosticsManager&&) = delete;
};

void Hi5DiagnosticsManager::init() {}

void Hi5DiagnosticsManager::log(Diagnostic&) {}

} // namespace PHARE


#endif //  PHARE_DIAGNOSTIC_HI5_DIAGNOSTIC_HPP_