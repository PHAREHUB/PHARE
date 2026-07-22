#ifndef PHARE_CORE_LOGGER_HPP
#define PHARE_CORE_LOGGER_HPP

#include "core/def.hpp" // IWYU pragma: keep

#include <string>
#include <cstdint>
#include <utility>
#include <memory>
#include <fstream>
#include <iostream>

#if !defined(PHARE_LOG_LEVEL)
#define PHARE_LOG_LEVEL 0 // 0 == off
#endif

namespace PHARE
{
constexpr static std::uint8_t LOG_LEVEL = PHARE_LOG_LEVEL;
}

#if !defined(NDEBUG) || defined(PHARE_FORCE_DEBUG_DO) || defined(PHARE_FORCE_LOG_LINE)
#include <sstream>  // IWYU pragma: keep
#include <iostream> // IWYU pragma: keep
#define PHARE_LOG_LINE_STR(str)                                                                    \
    std::cout << __FILE__ << ":" << __LINE__ << " - " << str << std::endl;
#define PHARE_LOG_LINE_SS(s) PHARE_LOG_LINE_STR((std::stringstream{} << s).str());
#else
#define PHARE_LOG_LINE_STR(str)
#define PHARE_LOG_LINE_SS(str)
#endif
#define PHARE_LOG_LINE PHARE_LOG_LINE_STR("")

#if PHARE_WITH_CALIPER
#include "caliper/cali.h" // IWYU pragma: keep

#define PHARE_LOG_START(lvl, str) CALI_MARK_BEGIN(str)
#define PHARE_LOG_STOP(lvl, str) CALI_MARK_END(str)
#define PHARE_LOG_SCOPE(lvl, str) PHARE::scope_log PHARE_STR_CAT(__phare_scope, __LINE__)(lvl, str)

#else // !PHARE_WITH_CALIPER

#include "core/utilities/logger/logger_defaults.hpp" // IWYU pragma: keep


#endif // PHARE_WITH_CALIPER

namespace PHARE
{
// owns log_out and redirects std::cout to it for this object's lifetime, restoring
// the original streambuf on destruction. Ownership of the stream lives here (rather
// than in a separate member this class merely references) so the redirected buffer
// can never outlive the object that installed it. Declared as a class member (not
// left to a destructor body) so restoration still happens if a later member throws
// mid-construction - a throwing constructor never runs the class's own destructor,
// only already-constructed members' destructors.
class CoutRedirect
{
public:
    explicit CoutRedirect(std::unique_ptr<std::ofstream> log_out)
        : log_out_{std::move(log_out)}
    {
        if (log_out_)
            saved_ = std::cout.rdbuf(log_out_->rdbuf());
    }

    ~CoutRedirect()
    {
        if (saved_)
            std::cout.rdbuf(saved_);
    }

    CoutRedirect(CoutRedirect const&)            = delete;
    CoutRedirect& operator=(CoutRedirect const&) = delete;

private:
    std::unique_ptr<std::ofstream> log_out_;
    std::streambuf* saved_ = nullptr;
};

struct scope_log
{
    scope_log(int&& i_, std::string&& str)
        : i{i_}
        , key{std::move(str)}
    {
        if (i <= LOG_LEVEL)
        {
            PHARE_LOG_START(i, key.c_str());
        }
    }
    ~scope_log()
    {
        if (i <= LOG_LEVEL)
        {
            PHARE_LOG_STOP(i, key.c_str());
        }
    }

    int i;
    std::string key;
};
} // namespace PHARE

#endif /* PHARE_CORE_LOGGER_H */
