#ifndef PHARE_CORE_MODELS_EXTERNAL_FIELD_UPDATER_DEFS_HPP
#define PHARE_CORE_MODELS_EXTERNAL_FIELD_UPDATER_DEFS_HPP

#include <cstdint>

namespace PHARE::core
{
enum class ExternalFieldUpdaterType : std::uint16_t { Zero, Dipole, UserDefined };

} // namespace PHARE::core

#endif // PHARE_CORE_MODELS_EXTERNAL_FIELD_UPDATER_DEFS_HPP
