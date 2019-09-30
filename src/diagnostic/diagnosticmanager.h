


#ifndef PHARE_DIAGNOSTIC_MANAGER_HPP_
#define PHARE_DIAGNOSTIC_MANAGER_HPP_


namespace PHARE
{
class Diagnostic
{
public:
    Diagnostic() {}
    virtual void dump() = 0;

protected:
    virtual ~Diagnostic() {}

private:
    Diagnostic(const Diagnostic&)             = delete;
    Diagnostic(const Diagnostic&&)            = delete;
    Diagnostic& operator&(const Diagnostic&)  = delete;
    Diagnostic& operator&(const Diagnostic&&) = delete;
};

class DiagnosticsManager
{
public:
    DiagnosticsManager() {}
    virtual void init()           = 0;
    virtual void log(Diagnostic&) = 0;

protected:
    virtual ~DiagnosticsManager() {}

private:
    DiagnosticsManager(const DiagnosticsManager&)             = delete;
    DiagnosticsManager(const DiagnosticsManager&&)            = delete;
    DiagnosticsManager& operator&(const DiagnosticsManager&)  = delete;
    DiagnosticsManager& operator&(const DiagnosticsManager&&) = delete;
};


} // namespace PHARE



#endif //  PHARE_DIAGNOSTIC_MANAGER_HPP_