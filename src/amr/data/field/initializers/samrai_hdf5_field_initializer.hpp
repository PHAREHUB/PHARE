#ifndef _PHARE_AMR_DATA_FIELD_INITIAZILIZERS_SAMRAI_HDF5_INITIALIZER_HPP_
#define _PHARE_AMR_DATA_FIELD_INITIAZILIZERS_SAMRAI_HDF5_INITIALIZER_HPP_

#include <memory>
#include <random>
#include <cassert>
#include <functional>

#include "core/def.hpp"
#include "core/logger.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/data/ndarray/ndarray_vector.hpp"
#include "core/hybrid/hybrid_quantities.hpp"
#include "core/utilities/types.hpp"
#include "core/data/ions/particle_initializers/particle_initializer.hpp"
#include "core/data/particles/particle.hpp"
#include "initializer/data_provider.hpp"
#include "core/utilities/point/point.hpp"

#include "hdf5/detail/h5/h5_file.hpp"


#include "SAMRAI/hier/PatchDataRestartManager.h"

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
    bool static constexpr c_ordering = false;
    // using Viewer                     = core::NdArrayViewer<Field_t::dimension>;

    auto const local_cell = [&](auto const& box, auto const& point) {
        core::Point<std::uint32_t, dimension> localPoint;
        auto localStart = layout.physicalStartIndex(core::QtyCentering::dual, core::Direction::X);
        for (std::size_t i = 0; i < dimension; ++i)
            localPoint[i] = point[i] - (box.lower[i] - localStart);
        return localPoint;
    };

    // PHARE_LOG_LINE_STR("SamraiHDF5FieldInitializer::load");

    auto const& centering = layout.centering(field.physicalQuantity());
    auto const& dest_box  = layout.AMRBox(); // grow(layout.AMRBox(), GridLayout::nbrGhosts());
    auto const& overlaps  = SamraiH5Interface<GridLayout>::INSTANCE().box_intersections(dest_box);

    // PHARE_LOG_LINE_STR(layout.AMRBox());
    for (auto const& [h5FilePtr, pdataptr] : overlaps)
    {
        auto& h5File = *h5FilePtr;
        auto& pdata  = *pdataptr;
        // PHARE_LOG_LINE_STR(pdata.box);

        auto const src_box = [&]() {
            auto copy = pdata.box;
            // for (std::size_t i = 0; i < Field_t::dimension; ++i)
            //     copy.upper[i] += centering[i] == core::QtyCentering::primal ? 1 : 0;
            return copy;
        }();

        std::vector<double> data;
        std::string const fieldpath
            = pdata.base_path + "/" + field.name() + "##default/field_" + field.name();
        h5File.file().getDataSet(fieldpath).read(data);

        // PHARE_LOG_LINE_STR(data.size());

        core::Box<std::uint32_t, GridLayout::dimension> const lcl_src_box{
            core::Point{core::ConstArray<std::uint32_t, GridLayout::dimension>()},
            core::Point{
                core::for_N<GridLayout::dimension, core::for_N_R_mode::make_array>([&](auto i) {
                    return static_cast<std::uint32_t>(
                        src_box.upper[i] - src_box.lower[i] + (GridLayout::nbrGhosts() * 2)
                        + (centering[i] == core::QtyCentering::primal ? 1 : 0));
                })}};

        // PHARE_LOG_LINE_STR(lcl_src_box.size());

        assert(data.size() == lcl_src_box.size());

        auto data_view = core::make_array_view<c_ordering>(data.data(), *lcl_src_box.shape());
        auto dst_iter  = dest_box.begin();
        auto src_iter  = src_box.begin();
        for (; src_iter != src_box.end(); ++src_iter)
        {
            // for (auto const& el : **src_iter)
            //     if (el < 0)
            //         continue;
            // if (core::any(*src_iter, [](auto const& el) { el < 0; }))
            //     continue; // skip ghosts

            if (isIn(core::Point{*src_iter}, dest_box))
            {
                while (*dst_iter != *src_iter and dst_iter != dest_box.end())
                    ++dst_iter;

                if (dst_iter == dest_box.end() or not isIn(core::Point{*dst_iter}, dest_box))
                    break;

                // auto src_idx = Viewer::idx(lcl_src_box.shape(), *src_iter);
                // auto dst_idx = Viewer::idx(dest_box.shape(), *dst_iter);

                // PHARE_LOG_LINE_STR(src_idx);
                // PHARE_LOG_LINE_STR(dst_idx);

                field(local_cell(dest_box, *dst_iter)) = data_view(local_cell(src_box, *src_iter));

                // PHARE_LOG_LINE_STR(*src_iter);
                // PHARE_LOG_LINE_STR(*dst_iter);

                // PHARE_LOG_LINE_STR(field(local_cell(dest_box, *dst_iter)));
                // PHARE_LOG_LINE_STR(data_view(local_cell(src_box, *src_iter)));
            }
        }
    }
}



} // namespace PHARE::amr


#endif
