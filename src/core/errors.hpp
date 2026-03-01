#ifndef PHARE_CORE_ERRORS_HPP
#define PHARE_CORE_ERRORS_HPP

#include "core/def.hpp"

#include <cfenv>
#include <atomic>
#include <string>
#include <csignal>
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <unordered_map>


#include "dict.hpp"


namespace PHARE::core
{
class DictionaryException : public std::exception
{
public:
    DictionaryException() = default;
    DictionaryException(auto const& k, auto const& v) { (*this)(k, v); }

    using Dict_t = cppdict::Dict<std::string>;

    auto& operator[](std::string const& key) { return dict_[key]; }
    auto& operator[](std::string const& key) const { return dict_[key]; }
    auto& operator()(std::string const key, std::string const val)
    {
        dict_[key] = val;
        return *this;
    }

    std::string operator()() const
    {
        std::stringstream ss;
        dict_.visit(
            [&](std::string const& key, auto const& val) { ss << key << " " << val << std::endl; });
        return ss.str();
    }

    char const* what() const noexcept override
    {
        static thread_local std::string msg;
        return (msg = (*this)()).c_str();
    }

    std::string id() const
    {
        if (!dict_.contains("ID"))
            throw std::runtime_error("DictionaryException has no ID");
        return (*this)["ID"].template to<std::string>();
    }

private:
    Dict_t dict_;
};


class Errors
{
public:
    NO_DISCARD static Errors& instance()
    {
        static Errors i;
        return i;
    }

    NO_DISCARD bool any() { return errors.size() > 0; }

    void log(std::string const key, std::string const val)
    {
        if (!errors.count(key))
        {
            std::cerr << key << " :" << val << std::endl;
            errors.emplace(key, val);
            error_count.emplace(key, 0);
        }
        error_count[key] = error_count[key] + 1;
    }


private:
    std::unordered_map<std::string, std::string> errors;
    std::unordered_map<std::string, std::size_t> error_count;
};

} // namespace PHARE::core


#if !defined(PHARE_LOG_ERROR)
#define PHARE_LOG_ERROR(x)                                                                         \
    PHARE::core::Errors::instance().log(std::string{__FILE__} + ":" + std::to_string(__LINE__), x);
#endif



#if PHARE_CATCH_FPE

namespace PHARE
{

int static constexpr default_fpe_exceptions()
{
#if defined(PHARE_FPE_ALL)
    return FE_ALL_EXCEPT;
#elif defined(PHARE_FPE_CUSTOM)
    return PHARE_FPE_CUSTOM;
#else // default
    return FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW;
#endif
}

struct FPEWatcher
{
    static void enable()
    {
        ++depth;
        if (depth > 1)
            return;
        feclearexcept(FE_ALL_EXCEPT);
        feenableexcept(exceptions);
    }

    static void disable()
    {
        --depth;
        if (depth > 0)
            return;

        fedisableexcept(exceptions);
    }

    static inline std::atomic<std::size_t> depth = 0;
    static int constexpr exceptions              = default_fpe_exceptions();

    FPEWatcher() { enable(); }
    ~FPEWatcher() { disable(); }
};

} // namespace PHARE



#define PHARE_FPE_START PHARE::FPEWatcher::enable()
#define PHARE_FPE_STOP PHARE::FPEWatcher::disable()
#define PHARE_FPE_SCOPE                                                                            \
    PHARE::FPEWatcher PHARE_STR_CAT(__phare_fpe, __LINE__) {}

#else // !PHARE_CATCH_FPE

#define PHARE_FPE_START
#define PHARE_FPE_STOP
#define PHARE_FPE_SCOPE

#endif // PHARE_FPE

#endif /* PHARE_CORE_ERRORS_H */
