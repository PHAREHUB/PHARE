#ifndef PHARE_CORE_UTILITIES_RESOURCES_VARIANTS_HPP
#define PHARE_CORE_UTILITIES_RESOURCES_VARIANTS_HPP

#include "core/utilities/types.hpp"

#include <core/logger.hpp>
#include <tuple>
#include <variant>
#include <stdexcept>

namespace PHARE::core
{
template<typename T>
auto decay_to_ptr()
{
    return [](T& arg) mutable -> void* { return const_cast<std::decay_t<T>*>(&arg); };
}

template<typename... Ts>
struct varient_visitor_overloads : Ts...
{
    using Ts::operator()...;
};

template<typename... Ts>
varient_visitor_overloads(Ts&&...) -> varient_visitor_overloads<std::decay_t<Ts>...>;


template<typename... Args>
auto constexpr _visit_ptr_overloads(std::tuple<Args...>*)
{
    return varient_visitor_overloads{decay_to_ptr<Args>()...,
                                     [](auto&) mutable -> void* { return nullptr; }};
}


template<typename T, typename... Ts>
struct unique : std::type_identity<T>
{
};

template<typename... Ts, typename U, typename... Us>
struct unique<std::tuple<Ts...>, U, Us...>
    : std::conditional_t<(std::is_same_v<U, Ts> || ...), unique<std::tuple<Ts...>, Us...>,
                         unique<std::tuple<Ts..., U>, Us...>>
{
};

template<typename... Ts>
using unique_tuple = typename unique<std::tuple<>, Ts...>::type;


template<typename... Args>
auto constexpr visit_ptr_overloads()
{
    return _visit_ptr_overloads(static_cast<unique_tuple<Args...>*>(nullptr));
}


template<typename Type, typename Variants>
auto& get_as_ref_or_throw(Variants& variants, std::size_t const start = 0)
{
    for (std::size_t idx = start; idx < variants.size(); ++idx)
        if (auto type = std::visit(visit_ptr_overloads<Type>(), variants[idx]))
            return *reinterpret_cast<Type*>(type);

    throw std::runtime_error("No element in variant for type");
}


template<typename Type>
auto& get_from_variants(auto& variants, Type& arg)
{
    std::size_t start = 0;

    while (start < variants.size())
    {
        if (auto& res = get_as_ref_or_throw<Type>(variants, start); res.name() == arg.name())
            return res;
        ++start;
    }

    throw std::runtime_error("Required name not found in variants: " + arg.name());
}


template<typename... Args>
auto get_from_variants(auto& variants, Args&... args)
    requires(sizeof...(Args) > 1)
{
    return std::forward_as_tuple(get_from_variants(variants, args)...);
}


} // namespace PHARE::core


#endif /*PHARE_CORE_UTILITIES_RESOURCES_VARIANTS_HPP*/
