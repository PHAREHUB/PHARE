#ifndef PHARE_CORE_UTILITIES_STRINGS_HPP
#define PHARE_CORE_UTILITIES_STRINGS_HPP

#include "core/def.hpp"

#include <string>
#include <sstream>
#include <iomanip>
#include <stdexcept>

namespace PHARE::core
{
template<typename T>
NO_DISCARD std::string to_string_with_precision(T const& a_value, std::size_t const len)
{
    std::ostringstream out;
    out.precision(len);
    out << std::fixed << a_value;
    auto str = out.str();
    // last digit may be non 0 because of rounding
    // and we know at that decimal it should be so we force it
    str.replace(str.end() - 1, str.end(), 1, '0');
    return out.str();
}

template<typename T>
NO_DISCARD auto to_string_fixed_width(T const& value, std::size_t const& precision,
                                      std::size_t const& width, char const& fill = '0')
{
    std::ostringstream out;
    out.width(width);
    out.precision(precision);
    out << std::setfill(fill) << std::fixed << value;
    return out.str();
}


template<typename T>
NO_DISCARD T from_string(std::string const& s)
{
    T t;
    std::stringstream ss(s);
    ss >> t;
    if (ss.fail())
        throw std::runtime_error("Cannot Parse T from: " + s);
    return t;
}

} // namespace PHARE::core

#endif // PHARE_CORE_UTILITIES_STRINGS_HPP
