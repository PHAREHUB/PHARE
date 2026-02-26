#ifndef PHARE_CORE_ERRORS_HPP
#define PHARE_CORE_ERRORS_HPP

#include "core/def.hpp"

#include <string>
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




#endif /* PHARE_CORE_ERRORS_H */
