#ifndef PHARE_CORE_DATA_TENSOR_FIELD_TRAITS
#define PHARE_CORE_DATA_TENSOR_FIELD_TRAITS

#include "core/data/field/field_traits.hpp"
#include "core/data/vecfield/vecfield_component.hpp"

#include <concepts>
#include <string>

namespace PHARE::core
{
/**
 * @brief Define the requirements for a tensor field type.
 *
 * @see TensorField
 */
template<typename T>
concept IsTensorField = requires(T tf, T const ctf, T const& crtf, Component component, size_t i) {
    requires IsField<typename T::field_type>;
    typename T::value_type;
    typename T::tensor_t;

    requires std::same_as<decltype(T::dimension), std::size_t const>;
    requires std::same_as<decltype(T::rank), std::size_t const>;
    { T::size() } -> std::convertible_to<std::size_t>;
    requires std::bool_constant<(T::size() >= 0)>::value;

    { tf.name() } -> std::same_as<std::string const&>;
    { tf.getComponent(component) } -> std::same_as<typename T::field_type&>;
    { ctf.getComponent(component) } -> std::same_as<typename T::field_type const&>;
    { tf(component) } -> std::same_as<typename T::field_type&>;
    { tf.getComponentName(component) } -> std::same_as<std::string>;
    { tf[i] } -> std::same_as<typename T::field_type&>;
    // missing 'components' overloads
    { tf.copyData(crtf) } -> std::same_as<void>;

    {
        tf.begin()
    } -> std::same_as<typename std::array<typename T::field_type, T::size()>::iterator>;
    { tf.end() } -> std::same_as<decltype(tf.begin())>;
    {
        ctf.begin()
    } -> std::same_as<typename std::array<typename T::field_type, T::size()>::const_iterator>;
    { ctf.end() } -> std::same_as<decltype(ctf.begin())>;

    { ctf.componentNames() } -> std::same_as<std::array<std::string, T::size()> const&>;
    { ctf.physicalQuantity() } -> std::same_as<typename T::tensor_t const&>;
};

/**
 * @brief A type verifying this concept is either a Field or a TensorField.
 */
template<typename T>
concept IsScalarOrTensorField = IsField<T> || IsTensorField<T>;


/**
 * @brief Select the physical quantity type based on tensoriality.
 */
template<typename ScalarOrTensorFieldT, bool is_scalar>
struct PhysicalQuantityTypeSelector;
/** @brief Specialization for scalar fields */
template<typename ScalarOrTensorFieldT>
struct PhysicalQuantityTypeSelector<ScalarOrTensorFieldT, true>
{
    using type = ScalarOrTensorFieldT::physical_quantity_type;
};
/** @brief Specialization for tensor fields */
template<typename ScalarOrTensorFieldT>
struct PhysicalQuantityTypeSelector<ScalarOrTensorFieldT, false>
{
    using type = ScalarOrTensorFieldT::tensor_t;
};


/**
 * @brief Select the underlying field type based on tensoriality.
 */
template<typename ScalarOrTensorFieldT, bool is_scalar>
struct FieldTypeSelector;
/** @brief Specialization for scalar fields */
template<typename ScalarOrTensorFieldT>
struct FieldTypeSelector<ScalarOrTensorFieldT, true>
{
    using type = ScalarOrTensorFieldT;
};
/** @brief Specialization for tensor fields */
template<typename ScalarOrTensorFieldT>
struct FieldTypeSelector<ScalarOrTensorFieldT, false>
{
    using type = ScalarOrTensorFieldT::field_type;
};


/**
 * @brief Select the underlying field type based on tensoriality.
 */
template<typename ScalarOrTensorFieldT, bool is_scalar>
struct NumberOfComponentsSelector;
/** @brief Specialization for scalar fields */
template<typename ScalarOrTensorFieldT>
struct NumberOfComponentsSelector<ScalarOrTensorFieldT, true>
{
    static constexpr size_t value = 1;
};
/** @brief Specialization for tensor fields */
template<typename ScalarOrTensorFieldT>
struct NumberOfComponentsSelector<ScalarOrTensorFieldT, false>
{
    static constexpr size_t value = ScalarOrTensorFieldT::size();
};

} // namespace PHARE::core

#endif // PHARE_CORE_DATA_TENSOR_FIELD_TRAITS
