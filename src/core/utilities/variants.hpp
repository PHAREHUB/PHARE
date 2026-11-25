#ifndef PHARE_CORE_UTILITIES_RESOURCES_VARIANTS_HPP
#define PHARE_CORE_UTILITIES_RESOURCES_VARIANTS_HPP

#include <tuple>
#include <string>
#include <variant>
#include <stdexcept>
#include <type_traits>



namespace PHARE::core
{
template<typename T>
concept std_variant
    = requires(std::decay_t<T>& t) { []<typename... Ts>(std::variant<Ts...>&) {}(t); };


template<typename... Ts>
struct variant_visitor_overloads : Ts...
{
    using Ts::operator()...;
};

template<typename... Ts>
variant_visitor_overloads(Ts&&...) -> variant_visitor_overloads<std::decay_t<Ts>...>;


// first element of type Type in variants, from index start
template<typename Type, typename Variants>
auto& get_as_ref_or_throw(Variants& variants, std::size_t const start = 0)
{
    for (std::size_t idx = start; idx < variants.size(); ++idx)
        if (auto ptr = std::get_if<std::decay_t<Type>>(&variants[idx]))
            return *ptr;

    throw std::runtime_error("No element in variant for type");
}


// the element of type Type in variants with the same name as arg
template<typename Type>
auto& get_from_variants(auto& variants, Type& arg)
{
    for (std::size_t idx = 0; idx < variants.size(); ++idx)
        if (auto ptr = std::get_if<std::decay_t<Type>>(&variants[idx]);
            ptr and ptr->name() == arg.name())
            return *ptr;

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
