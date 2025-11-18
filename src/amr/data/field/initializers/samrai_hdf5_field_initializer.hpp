#ifndef _PHARE_AMR_DATA_FIELD_INITIALIZERS_SAMRAI_HDF5_INITIALIZER_HPP_
#define _PHARE_AMR_DATA_FIELD_INITIALIZERS_SAMRAI_HDF5_INITIALIZER_HPP_


#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/data/ndarray/ndarray_vector.hpp"

#include "core/utilities/types.hpp"
#include "core/utilities/point/point.hpp"

#if !PHARE_HAS_HIGHFIVE
#error // HIGHFIVE REQUIRED!
#endif // PHARE_HAS_HIGHFIVE

#include "amr/data/initializers/samrai_hdf5_initializer.hpp"

#include <cassert>


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
        core::Box<std::uint32_t, GridLayout::dimension> const lcl_src_gbox{
            core::Point{core::ConstArray<std::uint32_t, GridLayout::dimension>()},
            core::Point{
                core::for_N<GridLayout::dimension, core::for_N_R_mode::make_array>([&](auto i) {
                    return static_cast<std::uint32_t>(
                        src_box.upper[i] - src_box.lower[i] + (GridLayout::nbrGhosts() * 2)
                        + (centering[i] == core::QtyCentering::primal ? 1 : 0));
                })}};
        auto const data_view   = core::make_array_view(data.data(), *lcl_src_gbox.shape());
        auto const overlap_gb  = grow(overlap_box, GridLayout::nbrGhosts());
        auto const lcl_src_box = layout.AMRToLocal(overlap_gb, src_box);
        auto const lcl_dst_box = layout.AMRToLocal(overlap_gb, dest_box);
        auto src_it            = lcl_src_box.begin();
        auto dst_it            = lcl_dst_box.begin();
        for (; src_it != lcl_src_box.end(); ++src_it, ++dst_it)
            field(*dst_it) = data_view(*src_it);
    }
}



} // namespace PHARE::amr



#endif /*_PHARE_AMR_DATA_FIELD_INITIALIZERS_SAMRAI_HDF5_INITIALIZER_HPP_*/
