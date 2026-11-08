// IWYU pragma: private, include "core/logger.hpp"

#ifndef PHARE_CORE_UTILITIES_LOGGER_LOGGER_DEFAULTS_HPP
#define PHARE_CORE_UTILITIES_LOGGER_LOGGER_DEFAULTS_HPP

#include "core/def/phlop.hpp"

#ifndef PHARE_SCOPE_TIMER
#define PHARE_SCOPE_TIMER(str) // nothing
#endif                         // PHARE_SCOPE_TIMER

#ifndef PHARE_SCOPE
#define PHARE_SCOPE(str, ...) // nothing
#endif                        // PHARE_SCOPE


#if PHARE_LOG_LEVEL >= 1
#define PHARE_LOG_SCOPE_1(str) PHARE_SCOPE(str, PHARE_LOG_SCOPE_FLAGS_1)
#else
#define PHARE_LOG_SCOPE_1(str)
#endif // LOG_LEVEL >= 1


#if PHARE_LOG_LEVEL >= 2
#define PHARE_LOG_SCOPE_2(str) PHARE_SCOPE(str, PHARE_LOG_SCOPE_FLAGS_2)
#else
#define PHARE_LOG_SCOPE_2(str)
#endif // LOG_LEVEL >= 2


#if PHARE_LOG_LEVEL == 3
#define PHARE_LOG_SCOPE_3(str) PHARE_SCOPE(str, PHARE_LOG_SCOPE_FLAGS_3)
#else
#define PHARE_LOG_SCOPE_3(str)
#endif // LOG_LEVEL == 3


#define PHARE_LOG_START(lvl, str)
#define PHARE_LOG_STOP(lvl, str)
#define PHARE_LOG_SCOPE(lvl, str) PHARE_STR_CAT(PHARE_LOG_SCOPE_, lvl)(str)


#endif /* PHARE_CORE_UTILITIES_LOGGER_LOGGER_DEFAULTS_HPP */
