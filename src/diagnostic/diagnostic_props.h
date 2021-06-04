#ifndef DIAGNOSTIC_DAO_H
#define DIAGNOSTIC_DAO_H

#include <string>
#include <vector>
#include <cstddef>

#include "cppdict/include/dict.hpp"

namespace PHARE::diagnostic
{
struct DiagnosticProperties
{
    // Types limited to actual need, no harm to modify
    using Params         = cppdict::Dict<std::size_t>;
    using FileAttributes = cppdict::Dict<std::string>;

    std::vector<double> writeTimestamps, computeTimestamps;
    std::string type, quantity;

    Params params{}; // supports arbitrary values for specific diagnostic writers
                     // for instance "flushEvery" for H5 file writers

    auto& operator[](std::string const& paramKey) { return params[paramKey]; }

    template<typename T>
    auto const& param(std::string const& paramKey) const
    {
        return params[paramKey].template to<T>();
    }

    std::size_t nAttributes = 0, dumpIdx = 0;
    FileAttributes fileAttributes{};
};

} // namespace PHARE::diagnostic

#endif // DIAGNOSTIC_DAO_H
