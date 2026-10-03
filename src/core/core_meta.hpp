#ifndef PHARE_CORE_CORE_META_HPP
#define PHARE_CORE_CORE_META_HPP

namespace PHARE::core
{
template<typename T>
concept IsTensorField = requires(T const& t) { t.components(); };

// TensorField also has physicalQuantity(), so exclude it explicitly - IsField means
//  "scalar field", not merely "anything with a physicalQuantity()"
template<typename T>
concept IsField = requires(T const& t) { t.physicalQuantity(); } && !IsTensorField<T>;

template<typename T>
constexpr bool is_field_v = IsField<T>;

template<typename T>
constexpr bool is_tensor_field_v = IsTensorField<T>;

} // namespace PHARE::core

#endif // PHARE_CORE_CORE_META_HPP
