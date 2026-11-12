#ifndef PHARE_SIMULATOR_SIMULATOR_DIAGNOSTICS_HPP
#define PHARE_SIMULATOR_SIMULATOR_DIAGNOSTICS_HPP

#include "simulator/simulator_decl.hpp"

#include "diagnostic/diagnostics.hpp"

namespace PHARE
{

template<auto opts>
void Simulator<opts>::diagnostics_init(initializer::PHAREDict const& dict, auto& model)
{
    dMan = PHARE::diagnostic::DiagnosticsManagerResolver::make_unique(*hierarchy_, model, dict);

    if (dict.contains("fine_dump_lvl_max"))
    {
        auto fine_dump_lvl_max = dict["fine_dump_lvl_max"].template to<int>();

        if (fine_dump_lvl_max > 0)
        { // copy for later
            this->fineDumpLvlMax                  = static_cast<std::size_t>(fine_dump_lvl_max);
            functors_["pre_advance"]["fine_dump"] = [&](SimFunctorParams const& params) {
                std::size_t level_nbr = params["level_nbr"].template to<int>();
                auto timestamp        = params["timestamp"].template to<double>();

                if (this->fineDumpLvlMax >= level_nbr)
                    this->dMan->dump_level(level_nbr, timestamp);
            };
        }
    }

    this->allowEmergencyDumps = cppdict::get_value(dict, "allow_emergency_dumps", false);
}


template<auto opts>
void Simulator<opts>::hybrid_diagnostics_init(initializer::PHAREDict const& dict)
{
    if constexpr (has_hybrid_v<opts>)
        diagnostics_init(dict, *hyb_.model_);
}


template<auto opts>
void Simulator<opts>::mhd_diagnostics_init(initializer::PHAREDict const& dict)
{
    if constexpr (has_mhd_v<opts>)
        diagnostics_init(dict, *mhd_.model_);
}

} // namespace PHARE

#endif /*PHARE_SIMULATOR_SIMULATOR_DIAGNOSTICS_HPP*/
