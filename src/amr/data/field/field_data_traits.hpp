#ifndef PHARE_SRC_AMR_FIELD_FIELD_DATA_TRAITS_HPP
#define PHARE_SRC_AMR_FIELD_FIELD_DATA_TRAITS_HPP

#include <SAMRAI/hier/Patch.h>
#include <SAMRAI/hier/PatchData.h>

#include <concepts>

namespace PHARE::amr
{
/**
 * @brief Concept ensuring a type satisfies the PHARE FieldData interface.
 */
template<typename T>
concept IsFieldData
    = std::derived_from<T, SAMRAI::hier::PatchData>
      && requires(T a, T const ca, SAMRAI::hier::Patch const& patch) {
             // Type aliases
             typename T::gridlayout_type;
             typename T::grid_type;
             typename T::physical_quantity_type;

             // Static constexpr variables
             requires std::same_as<decltype(T::dimension), std::size_t const>;
             requires std::same_as<decltype(T::interp_order), std::size_t const>;

             // Public member variables
             requires std::same_as<decltype(a.gridLayout), typename T::gridlayout_type>;
             requires std::same_as<decltype(a.field), typename T::grid_type>;

             // API requirements
             { a.getPointer() } -> std::same_as<typename T::grid_type::field_type*>;
             { T::getLayout(patch, 0) } -> std::same_as<typename T::gridlayout_type const&>;
             { T::getField(patch, 0) } -> std::same_as<typename T::grid_type&>;
         };

} // namespace PHARE::amr

#endif // PHARE_SRC_AMR_FIELD_FIELD_DATA_TRAITS_HPP
