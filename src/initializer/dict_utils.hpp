#ifndef PHARE_INITIALIZER_DICT_UTILS_HPP
#define PHARE_INITIALIZER_DICT_UTILS_HPP

#include <array>
#include <string>

#include "initializer/data_provider.hpp"

namespace PHARE::initializer
{
/**
 * @brief Fill an array of @p dimension elements by reading sub-keys "x", "y", "z"
 *        from @p dict[@p key].
 *
 * Only the components up to @p dimension are read; the rest of the array is
 * left untouched.  This overload writes directly into a caller-supplied buffer.
 *
 * @tparam Type      Scalar type to extract (e.g. @c int, @c double, @c std::string).
 * @tparam dimension Number of spatial dimensions (1, 2, or 3).
 * @param dict  Source dictionary containing the nested x/y/z sub-dict.
 * @param key   Key within @p dict whose value holds the x/y/z sub-dict.
 * @param arr   Output buffer of at least @p dimension elements.
 */
template<typename Type, std::size_t dimension>
void parseDimXYZType(PHAREDict const& dict, std::string key, Type* arr)
{
    arr[0] = dict[key]["x"].template to<Type>();
    if constexpr (dimension > 1)
        arr[1] = dict[key]["y"].template to<Type>();
    if constexpr (dimension > 2)
        arr[2] = dict[key]["z"].template to<Type>();
}

/**
 * @brief Read sub-keys "x", "y", "z" from @p dict[@p key] into a
 *        @c std::array of @p dimension elements.
 *
 * @tparam Type      Scalar type to extract (e.g. @c int, @c double, @c std::string).
 * @tparam dimension Number of spatial dimensions (1, 2, or 3).
 * @param dict  Source dictionary containing the nested x/y/z sub-dict.
 * @param key   Key within @p dict whose value holds the x/y/z sub-dict.
 * @return A @c std::array<Type, dimension> filled with the x[/y[/z]] values.
 */
template<typename Type, std::size_t dimension>
auto parseDimXYZType(PHAREDict const& dict, std::string key)
{
    std::array<Type, dimension> arr;
    parseDimXYZType<Type, dimension>(dict, key, arr.data());
    return arr;
}

} // namespace PHARE::initializer

#endif // PHARE_INITIALIZER_DICT_UTILS_HPP
