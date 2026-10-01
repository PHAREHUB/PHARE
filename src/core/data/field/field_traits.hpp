#ifndef PHARE_CORE_DATA_FIELD_FIELD_TRAITS
#define PHARE_CORE_DATA_FIELD_FIELD_TRAITS

#include <concepts>
#include <cstddef>
#include <string>

namespace PHARE::core
{

/**
 * @brief Define the requirements for a Field type.
 *
 * A Field must provide static metadata, identification methods,
 * memory access, and dimension-specific indexing operators.
 */
template<typename T>
concept IsField = requires(T field) {
    { T::dimension } -> std::convertible_to<std::size_t>;
    typename T::value_type;
    typename T::physical_quantity_type;

    { field.name() } -> std::convertible_to<std::string const&>;
    { field.physicalQuantity() } -> std::same_as<typename T::physical_quantity_type const&>;

    { field.data() } -> std::same_as<typename T::value_type*>; // Inherited from NdArrayView

    requires((T::dimension == 1 && requires(T f) {
                 { f(std::declval<std::size_t>()) } -> std::same_as<typename T::value_type&>;
             }) || (T::dimension == 2 && requires(T f) {
                 {
                     f(std::declval<std::size_t>(), std::declval<std::size_t>())
                 } -> std::same_as<typename T::value_type&>;
             }) || (T::dimension == 3 && requires(T f) {
                 {
                     f(std::declval<std::size_t>(), std::declval<std::size_t>(),
                       std::declval<std::size_t>())
                 } -> std::same_as<typename T::value_type&>;
             }));
};

} // namespace PHARE::core

#endif // PHARE_CORE_DATA_FIELD_FIELD_TRAITS
