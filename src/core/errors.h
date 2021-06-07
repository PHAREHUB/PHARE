
#ifndef PHARE_CORE_ERRORS_H
#define PHARE_CORE_ERRORS_H

#include <vector>
#include <unordered_map>


namespace PHARE::core
{
class Errors
{
public:
    static Errors& instance()
    {
        static Errors i;
        return i;
    }

    bool any() { return errors.size() > 0; }

    void log(std::string const& key, std::string const& val)
    {
        if (!errors.count(key))
        {
            std::cerr << key << " :" << val << std::endl;
            errors.emplace(key, val);
            error_count.emplace(key, 0);
        }
        error_count[key] = error_count[key] + 1;
    }

    void log(std::string&& key, std::string&& val) { log(key, val); }


private:
    std::unordered_map<std::string, std::string> errors;
    std::unordered_map<std::string, std::size_t> error_count;
};

} // namespace PHARE::core


#if !defined(PHARE_LOG_ERROR)
#define PHARE_LOG_ERROR(x)                                                                         \
    PHARE::core::Errors::instance().log(std::string{__FILE__} + ":" + std::to_string(__LINE__), x);
#endif


#endif /* PHARE_CORE_ERRORS_H */
