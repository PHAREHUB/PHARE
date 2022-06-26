
#ifndef PHARE_CORE_LOGGER_HPP
#define PHARE_CORE_LOGGER_HPP


#if !defined(NDEBUG) || defined(PHARE_FORCE_DEBUG_DO)
#define PHARE_LOG_LINE_STR(str)                                                                    \
    std::cout << __FILE__ << ":" << __LINE__ << " - " << str << std::endl;
#else
#define PHARE_LOG_LINE_STR(str)
#endif
#define PHARE_LOG_LINE PHARE_LOG_LINE_STR("")




#if PHARE_WITH_CALIPER
#include "caliper/cali.h"

#define PHARE_LOG_START(str) CALI_MARK_BEGIN(str)
#define PHARE_LOG_STOP(str) CALI_MARK_END(str)
#define PHARE_LOG_SCOPE(str) PHARE::scope_log __phare_scope##__line__(str)

#else

#define PHARE_LOG_START(str)
#define PHARE_LOG_STOP(str)
#define PHARE_LOG_SCOPE(str)

#endif

#include <string>
#include <utility>

namespace PHARE
{
struct scope_log
{
    scope_log(std::string&& str)
        : key{std::move(str)}
    {
        PHARE_LOG_START(key.c_str());
    }
    ~scope_log() { PHARE_LOG_STOP(key.c_str()); }

    std::string key;
};
} // namespace PHARE

#endif /* PHARE_CORE_LOGGER_H */
