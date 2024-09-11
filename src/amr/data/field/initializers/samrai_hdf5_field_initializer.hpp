#ifndef _PHARE_AMR_DATA_FIELD_INITIAZILIZERS_SAMRAI_HDF5_INITIALIZER_HPP_
#define _PHARE_AMR_DATA_FIELD_INITIAZILIZERS_SAMRAI_HDF5_INITIALIZER_HPP_

#include <cassert>

#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/data/ndarray/ndarray_vector.hpp"

#include "core/utilities/types.hpp"
#include "core/utilities/point/point.hpp"

#include "amr/data/initializers/samrai_hdf5_initializer.hpp"

namespace PHARE::amr
{


template<typename Field_t, typename GridLayout>
class SamraiHDF5FieldInitializer // no virtual classes needed (yet)
{
public:
    static constexpr auto dimension = GridLayout::dimension;

    SamraiHDF5FieldInitializer() {}

    void load(Field_t& Field, GridLayout const& layout) const;
};


template<typename Field_t, typename GridLayout>
void SamraiHDF5FieldInitializer<Field_t, GridLayout>::load(Field_t& field,
                                                           GridLayout const& layout) const
{
    auto const local_cell
        = [&](auto const& box, auto const& point) { return layout.AMRToLocal(point, box); };

    auto const& dest_box  = layout.AMRBox();
    auto const& centering = layout.centering(field.physicalQuantity());
    auto const& overlaps  = SamraiH5Interface<GridLayout>::INSTANCE().box_intersections(dest_box);
    for (auto const& [overlap_box, h5FilePtr, pdataptr] : overlaps)
    {
        auto& h5File       = *h5FilePtr;
        auto& pdata        = *pdataptr;
        auto const src_box = pdata.box;
        auto const data    = h5File.template read_data_set_flat<double>(
            pdata.base_path + "/" + field.name() + "##default/field_" + field.name());
        core::Box<std::uint32_t, GridLayout::dimension> const lcl_src_box{
            core::Point{core::ConstArray<std::uint32_t, GridLayout::dimension>()},
            core::Point{
                core::for_N<GridLayout::dimension, core::for_N_R_mode::make_array>([&](auto i) {
                    return static_cast<std::uint32_t>(
                        src_box.upper[i] - src_box.lower[i] + (GridLayout::nbrGhosts() * 2)
                        + (centering[i] == core::QtyCentering::primal ? 1 : 0));
                })}};
        auto data_view = core::make_array_view(data.data(), *lcl_src_box.shape());
        for (auto const& point : overlap_box)
            field(local_cell(dest_box, point)) = data_view(local_cell(src_box, point));
    }
}



} // namespace PHARE::amr



#endif
