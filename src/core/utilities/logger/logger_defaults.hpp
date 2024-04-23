#ifndef PHARE_CORE_LOGGER_DEFAULTS_HPP
#define PHARE_CORE_LOGGER_DEFAULTS_HPP

#include "core/utilities/run_timer.hpp"


#if PHARE_LOG_LEVEL == 5

#define PHARE_LOG_SCOPE_5(str)

#else

#define PHARE_LOG_SCOPE_5(str)

#endif // LOG_LEVEL == 5

#if PHARE_LOG_LEVEL >= 1

#define PHARE_LOG_SCOPE_1(str) RUN_TIMER_SCOPE(str)

#else

#define PHARE_LOG_SCOPE_1(str)

#endif // LOG_LEVEL >= 1


#define PHARE_LOG_START(lvl, str)
#define PHARE_LOG_STOP(lvl, str)
#define PHARE_LOG_SCOPE(lvl, str) STR_CAT(PHARE_LOG_SCOPE_, lvl)(str)


#endif /* PHARE_CORE_LOGGER_DEFAULTS_HPP */
