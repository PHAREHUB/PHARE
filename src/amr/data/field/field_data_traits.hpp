#ifndef PHARE_SRC_AMR_FIELD_FIELD_DATA_TRAITS_HPP
#define PHARE_SRC_AMR_FIELD_FIELD_DATA_TRAITS_HPP

#include "core/data/field/field_tiles.hpp"

#include <SAMRAI/hier/Patch.h>
#include <SAMRAI/hier/PatchData.h>

#include <concepts>
#include <type_traits>

namespace PHARE::amr
{
// what FieldData::getPointer returns: Grid's operator& yields its field view, a tile set is
// returned as is (its field_type is the per-tile field)
template<typename Grid_t>
using field_data_pointer_t = std::conditional_t<core::is_field_tile_set_v<Grid_t>, Grid_t*,
                                                typename Grid_t::field_type*>;

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

             // Public member variables
             requires std::same_as<decltype(a.gridLayout), typename T::gridlayout_type>;
             requires std::same_as<decltype(a.field), typename T::grid_type>;

             // API requirements
             { a.getPointer() } -> std::same_as<field_data_pointer_t<typename T::grid_type>>;
             { T::getLayout(patch, 0) } -> std::same_as<typename T::gridlayout_type const&>;
             { T::getField(patch, 0) } -> std::same_as<typename T::grid_type&>;
         };

} // namespace PHARE::amr

#endif // PHARE_SRC_AMR_FIELD_FIELD_DATA_TRAITS_HPP
