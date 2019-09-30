
#ifndef PHARE_DIAGNOSTIC_MANAGER_HPP_
#define PHARE_DIAGNOSTIC_MANAGER_HPP_

#include <unordered_set>
#include "data_provider.h"

namespace PHARE
{
namespace diagnostic
{
    enum class Level { INFO, DBG };
    enum class Mode { LIGHT, FULL };

    template<class T> // this is so we can use struct {} initializtion with shared_ptrs/forwarding
    struct aggregate_adapter : public T
    {
        template<class... Args>
        aggregate_adapter(Args&&... args)
            : T{std::forward<Args>(args)...}
        {
        }
    };
} // namespace diagnostic

struct Diagnostic
{
    size_t compute_every = 0, write_every = 0;
    size_t start_iteration = 0, end_iteration = 0; // likely to be time rather than index
                                                   // do we allow ranges?
    std::string name, species, type;
};

template<typename MODEL>
class ModelDiagnostic
{
public:
    ModelDiagnostic(MODEL& model)
        : model_(model)
    {
    }

protected:
    MODEL& model_;
};


struct DiagnosticStaticInfo
{
    static std::unordered_set<std::string> const types;
};
std::unordered_set<std::string> const DiagnosticStaticInfo::types{
    {"rho_s", "flux_s", "E", "B", "space_box"}};

template<typename DIAGNOSTIC, diagnostic::Level level = diagnostic::Level::INFO>
class DiagnosticsManager
{
public:
    using DAGG = diagnostic::aggregate_adapter<Diagnostic>;
    DiagnosticsManager(DIAGNOSTIC& ddd)
        : ddd_(ddd)
    {
    }
    void dump() { ddd_.template dump<level>(diagnostics); }

    DiagnosticsManager& addDiagnostic(std::shared_ptr<Diagnostic> const& diagnostic)
    {
        diagnostics.emplace_back(diagnostic);
        return *this;
    }

    DiagnosticsManager& addDiagDict(PHARE::initializer::PHAREDict<1>& dict);

    bool validateDiagDict(PHARE::initializer::PHAREDict<1>& dict);

    DiagnosticsManager& addDiagDict() // TODO REMOVE / this is garbage data
    {
        auto DIAG_DICT_IDX = diagnostics.size();
        PHARE::initializer::PHAREDict<1> dict;
        dict["diag"]["name"]            = std::string{"NAME_" + std::to_string(DIAG_DICT_IDX)};
        dict["diag"]["species"]         = std::string{"SPECIES_" + std::to_string(DIAG_DICT_IDX)};
        dict["diag"]["type"]            = std::string{"TYPE_" + std::to_string(DIAG_DICT_IDX)};
        dict["diag"]["compute_every"]   = std::size_t{5};
        dict["diag"]["write_every"]     = std::size_t{10};
        dict["diag"]["start_iteration"] = std::size_t{1};
        dict["diag"]["end_iteration"]   = std::size_t{100};
        DIAG_DICT_IDX++;
        return addDiagDict(dict);
    }

private:
    DIAGNOSTIC& ddd_;
    std::vector<std::shared_ptr<Diagnostic>> diagnostics;

    DiagnosticsManager(const DiagnosticsManager&)             = delete;
    DiagnosticsManager(const DiagnosticsManager&&)            = delete;
    DiagnosticsManager& operator&(const DiagnosticsManager&)  = delete;
    DiagnosticsManager& operator&(const DiagnosticsManager&&) = delete;
};

template<typename D, diagnostic::Level L>
DiagnosticsManager<D, L>&
    DiagnosticsManager<D, L>::addDiagDict(PHARE::initializer::PHAREDict<1>& dict)
{
    size_t &compute_every   = dict["diag"]["compute_every"].template to<std::size_t>(),
           &write_every     = dict["diag"]["write_every"].template to<std::size_t>(),
           &start_iteration = dict["diag"]["start_iteration"].template to<std::size_t>(),
           &end_iteration   = dict["diag"]["end_iteration"].template to<std::size_t>();
    std::string &name       = dict["diag"]["name"].template to<std::string>(),
                &species    = dict["diag"]["species"].template to<std::string>(),
                &type       = dict["diag"]["type"].template to<std::string>();

    diagnostics.emplace_back(std::make_shared<DAGG>(compute_every, write_every, start_iteration,
                                                    end_iteration, name, species, type));
    return *this;
}

template<typename D, diagnostic::Level L>
bool DiagnosticsManager<D, L>::validateDiagDict(PHARE::initializer::PHAREDict<1>& dict)
{
    return /*ddd_.validate(dict) &&*/ true;
}

} // namespace PHARE

#endif //  PHARE_DIAGNOSTIC_MANAGER_HPP_
