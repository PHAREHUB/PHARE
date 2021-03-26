
#ifndef PHARE_CORE_LOGGER_H
#define PHARE_CORE_LOGGER_H

#if PHARE_WITH_CALIPER
#include "caliper/cali.h"

#define PHARE_LOG_START(str) CALI_MARK_BEGIN(str);
#define PHARE_LOG_STOP(str) CALI_MARK_END(str);
#define PHARE_LOG_SCOPE(str) CALI_CXX_MARK_FUNCTION

#else

#define PHARE_LOG_START(str)
#define PHARE_LOG_STOP(str)
#define PHARE_LOG_SCOPE(str)

#endif

#endif /* PHARE_CORE_LOGGER_H */
