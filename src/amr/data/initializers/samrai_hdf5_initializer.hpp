#ifndef _PHARE_AMR_DATA_INITIALIZERS_SAMRAI_HDF5_INITIALIZER_HPP_
#define _PHARE_AMR_DATA_INITIALIZERS_SAMRAI_HDF5_INITIALIZER_HPP_

#include <memory>
#include <cassert>

#include "core/def.hpp"
#include "core/utilities/types.hpp"
#include "core/utilities/box/box.hpp"
#include "hdf5/detail/h5/h5_file.hpp"

#include "SAMRAI/tbox/HDFDatabase.h"


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

    void populate_from(std::string const& dir, int const& idx, int const& mpi_size,
                       std::string const& field_name = "field_EM_B_x");

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
        using Tup = std::tuple<Box_t, SamraiHDF5File const* const,
                               SamraiH5PatchDataInfo<GridLayout::dimension> const* const>;
        std::vector<Tup> overlaps;
        for (auto const& h5File : restart_files)
            for (auto const& patch : h5File->patches)
                if (auto const intersection = box * patch.box)
                    overlaps.emplace_back(*intersection, h5File.get(), &patch);
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
        , filepath_{filepath}
    {
    }

    auto getBoxFromPath(std::string const& path) const
    {
        std::size_t constexpr samrai_dim = 3; // always 3!
        auto constexpr _to_std_array     = [](auto& i) {
            return core::sized_array<GridLayout::dimension>(
                *reinterpret_cast<std::array<int, samrai_dim> const*>(&i));
        };

        SAMRAI::tbox::HDFDatabase db{"db"};
        db.open(filepath_);
        auto const boxes = db.getDatabaseBoxVector(path);
        return Box_t{_to_std_array(boxes[0].d_data.d_lo), _to_std_array(boxes[0].d_data.d_hi)};
    }

    std::string const filepath_;
    std::vector<SamraiH5PatchDataInfo<GridLayout::dimension>> patches;
};



template<typename GridLayout>
void SamraiH5Interface<GridLayout>::populate_from(std::string const& dir, int const& idx,
                                                  int const& mpi_size,
                                                  std::string const& field_name)
{
    if (restart_files.size()) // executed per pop, but we only need to run this once
        return;

    for (int rank = 0; rank < mpi_size; ++rank)
    {
        auto const hdf5_filepath = getRestartFileFullPath(dir, idx, mpi_size, rank);
        auto& h5File = *restart_files.emplace_back(std::make_unique<SamraiHDF5File>(hdf5_filepath));
        for (auto const& group : h5File.scan_for_groups({"level_0000", field_name}))
        {
            auto const field_path = group.substr(0, group.rfind("/"));
            auto const& field_box = h5File.getBoxFromPath(field_path + "/d_box");
            h5File.patches.emplace_back(field_box, field_path.substr(0, field_path.rfind("/")));
        }
    }
}



} // namespace PHARE::amr


#endif /*_PHARE_AMR_DATA_INITIALIZERS_SAMRAI_HDF5_INITIALIZER_HPP_*/
