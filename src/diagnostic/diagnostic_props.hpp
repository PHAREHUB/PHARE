#ifndef DIAGNOSTIC_DAO_HPP
#define DIAGNOSTIC_DAO_HPP

#include "core/def.hpp"
#include "core/utilities/types.hpp"

#include <string>
#include <type_traits>
#include <vector>
#include <cstddef>

#include "dict.hpp"

namespace PHARE::diagnostic
{
struct DiagnosticProperties
{
    // Types limited to actual need, no harm to modify
    using Params         = cppdict::Dict<std::size_t>;
    using FileAttributes = cppdict::Dict<std::string, double>;

    std::vector<double> writeTimestamps, computeTimestamps, elapsedTimestamps;
    std::string type, quantity;

    Params params{}; // supports arbitrary values for specific diagnostic writers
                     // for instance "flushEvery" for H5 file writers

    NO_DISCARD auto& operator[](std::string const& paramKey) { return params[paramKey]; }

    template<typename T>
    NO_DISCARD auto const& param(std::string const& paramKey) const
    {
        return params[paramKey].template to<T>();
    }

    void forward_file_attribute(std::string const& key, auto& dict_node);

    std::size_t nAttributes = 0, dumpIdx = 0;
    FileAttributes fileAttributes{};
};

template<typename... Ts0>
void forward_to_file_attributes(cppdict::Dict<Ts0...>& dst, std::string const& key, auto& dict_node)
{
    std::visit(
        [&](auto const& val) {
            using Type = std::decay_t<decltype(val)>;
            if constexpr (core::is_any_of<Type, Ts0...>())
                dst[key] = val;
        },
        dict_node.data);
}

void DiagnosticProperties::forward_file_attribute(std::string const& key, auto& dict_node)
{
    forward_to_file_attributes(fileAttributes, key, dict_node);
}


// fluid quantities the time loop does not keep up to date: computed when written, only when all
// levels are at the same time (amr::HybridLazyMoments)
NO_DISCARD inline bool isLazyFluidQuantity(DiagnosticProperties const& diagnostic)
{
    if (diagnostic.type != "fluid")
        return false;
    auto const& q       = diagnostic.quantity;
    auto const endsWith = [&](std::string const& end) {
        return q.size() >= end.size() and q.compare(q.size() - end.size(), end.size(), end) == 0;
    };
    // "/ions/pop/<name>/flux", "/ions/charge_flux", "/ions/bulkVelocity", and the momentum tensors
    return endsWith("/flux") or endsWith("/charge_flux") or endsWith("/bulkVelocity")
           or endsWith("/momentum_tensor");
}



} // namespace PHARE::diagnostic

#endif // DIAGNOSTIC_DAO_HPP
