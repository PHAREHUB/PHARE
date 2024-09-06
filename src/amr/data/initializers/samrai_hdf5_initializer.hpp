#ifndef _PHARE_AMR_DATA_INITIAZILIZERS_SAMRAI_HDF5_INITIALIZER_HPP_
#define _PHARE_AMR_DATA_INITIAZILIZERS_SAMRAI_HDF5_INITIALIZER_HPP_

#include <memory>
#include <random>
#include <cassert>
#include <functional>

#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/hybrid/hybrid_quantities.hpp"
#include "core/utilities/types.hpp"
#include "core/data/ions/particle_initializers/particle_initializer.hpp"
#include "core/data/particles/particle.hpp"
#include "initializer/data_provider.hpp"
#include "core/utilities/point/point.hpp"
#include "core/def.hpp"
#include "core/logger.hpp"

#include "hdf5/detail/h5/h5_file.hpp"
#include "amr/utilities/box/amr_box.hpp"

#include "SAMRAI/hier/PatchDataRestartManager.h"


namespace PHARE::amr
{

template<std::size_t dim>
struct SamraiH5PatchDataInfo
{
    using Box_t = core::Box<int, dim>;

    SamraiH5PatchDataInfo(Box_t const& b, std::string const& p)
        : box{b}
        , base_path{p}
    {
    }


    Box_t const box;
    std::string const base_path;
};



template<typename GridLayout>
class SamraiH5Interface
{
    struct SamraiHDF5File;

public:
    using Box_t = core::Box<int, GridLayout::dimension>;

    static SamraiH5Interface& INSTANCE()
    {
        static SamraiH5Interface i;
        return i;
    }

    void populate_from(std::string const& dir, int const& idx, int const& mpi_size);

    NO_DISCARD auto static getRestartFileFullPath(std::string path, int const& idx,
                                                  int const& mpi_size, int const& rank)
    {
        return path                                                          //
               + "/restore." + SAMRAI::tbox::Utilities::intToString(idx, 6)  //
               + "/nodes." + SAMRAI::tbox::Utilities::nodeToString(mpi_size) //
               + "/proc." + SAMRAI::tbox::Utilities::processorToString(rank);
    }

    auto box_intersections(Box_t const& box)
    {
        using Pair = std::pair<SamraiHDF5File const* const,
                               SamraiH5PatchDataInfo<GridLayout::dimension> const* const>;
        std::vector<Pair> overlaps;
        for (auto const& h5File : restart_files)
            for (auto const& patch : h5File->patches)
                if (auto const intersection = box * patch.box)
                    overlaps.emplace_back(h5File.get(), &patch);
        return overlaps;
    }



private:
    std::vector<std::unique_ptr<SamraiHDF5File>> restart_files;
    std::unordered_map<std::string, std::string> box2dataset;
};

template<typename GridLayout>
struct SamraiH5Interface<GridLayout>::SamraiHDF5File : public hdf5::h5::HighFiveFile
{
    using Super = hdf5::h5::HighFiveFile;

    SamraiHDF5File(std::string const& filepath)
        : Super{filepath, HighFive::File::ReadOnly, /*para=*/false}
    {
        // auto box_type = box_compound_type();
        // box_type.commit(file(), "d_box");

        // getBoxFromPath("/PHARE_hierarchy/level_0000/level_0000-patch_0000000-block_0000000/"
        //                "EM_B_y##default/d_box");
    }



    /*
       DATASET "d_box" {
          DATATYPE  H5T_COMPOUND {
             H5T_STD_I32BE "dim";
             H5T_ARRAY { [3] H5T_STD_I32BE } "lo";
             H5T_ARRAY { [3] H5T_STD_I32BE } "hi";
          }
          DATASPACE  SIMPLE { ( 1 ) / ( 1 ) }
          DATA {
          (0): {
                1,
                [ 0, 0, 0 ],
                [ 199, 0, 0 ]
             }
          }
          ATTRIBUTE "Type" {
             DATATYPE  H5T_STD_I8BE
             DATASPACE  SCALAR
             DATA {
             (0): 2
             }
          }
       }
    */



    struct BoxData
    {
        std::int32_t dim;
        std::array<std::int32_t, 3> lo;
        std::array<std::int32_t, 3> hi;
    };

    // struct BoxData
    // {
    //     int dim;
    //     int lo0, lo1, lo2;
    //     int hi0, hi1, hi2;
    // };

    HighFive::CompoundType static box_compound_type()
    {
        return {{"dim", HighFive::create_datatype<std::int32_t>()},
                {"lo", HighFive::create_datatype<std::array<std::int32_t, 3>>()},
                {"hi", HighFive::create_datatype<std::array<std::int32_t, 3>>()}};
    }

    auto getBoxFromPath(std::string const& path) const
    {
        // auto const& data = Super::read_data_set<int>(path);
        // PHARE_LOG_LINE_STR(data.size());
        std::vector<BoxData> boxes;
        Super::file().getDataSet(path).read(boxes);
        assert(boxes.size() == 1);

        // auto const& bx = boxes[0];
        // return Box_t{core::as_sized_array<GridLayout::dimension>(bx.lo0, bx.lo1, bx.lo2),
        //              core::as_sized_array<GridLayout::dimension>(bx.hi0, bx.hi1, bx.hi2)};

        return Box_t{core::sized_array<GridLayout::dimension>(boxes[0].lo),
                     core::sized_array<GridLayout::dimension>(boxes[0].hi)};
    }

    std::vector<SamraiH5PatchDataInfo<GridLayout::dimension>> patches;
};




template<typename GridLayout>
void SamraiH5Interface<GridLayout>::populate_from(std::string const& dir, int const& idx,
                                                  int const& mpi_size)
{
    Box_t const mock{{0}, {99}};
    for (int rank = 0; rank < mpi_size; ++rank)
    {
        auto const hdf5_filepath = getRestartFileFullPath(dir, idx, mpi_size, rank);
        auto& h5File = *restart_files.emplace_back(std::make_unique<SamraiHDF5File>(hdf5_filepath));
        for (auto const& group : h5File.scan_for_groups({"level_0000", "field_EM_B_x"}))
        {
            auto const em_path = group.substr(0, group.rfind("/"));
            h5File.patches.emplace_back(mock, em_path.substr(0, em_path.rfind("/")));
        }
    }
}



} // namespace PHARE::amr


#endif
