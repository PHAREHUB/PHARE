#ifndef PHARE_DIAGNOSTIC_DETAIL_VTK_H5_TYPE_WRITER_HPP
#define PHARE_DIAGNOSTIC_DETAIL_VTK_H5_TYPE_WRITER_HPP


#include "core/utilities/box/box.hpp"
#include "core/utilities/algorithm.hpp"
#include "core/utilities/mpi_utils.hpp"
#include "core/data/tensorfield/tensorfield.hpp"

#include "core/utilities/types.hpp"
#include "diagnostic/diagnostic_writer.hpp"

#include "hdf5/detail/h5/h5_file.hpp"

#include <string>
#include <stdexcept>
#include <unordered_map>



namespace PHARE::diagnostic::vtkh5::detail
{

// some testing shows 32 to be the max value supported, something to do with a value `H5S_MAX_RANK`
static inline auto const CHUNK_SIZE = core::get_env_as("PHARE_VTK_H5_CHUNK_SIZE", std::size_t{32});

} // namespace PHARE::diagnostic::vtkh5::detail


namespace PHARE::diagnostic::vtkh5
{

using namespace hdf5::h5;

template<typename Writer>
class H5TypeWriter : public PHARE::diagnostic::TypeWriter
{
    using FloatType                      = float; // Writer::FloatType;
    std::string static inline base       = "/VTKHDF/";
    std::string static inline level_base = base + "Level";
    std::string static inline step_level = base + "Steps/Level";
    using physical_quantity_type         = Writer::ModelView::Field::physical_quantity_type;

public:
    static constexpr auto dimension = Writer::dimension;
    using GridLayout                = Writer::GridLayout;
    using Box_t                     = core::Box<int, dimension>;

    H5TypeWriter(Writer& h5Writer)
        : h5Writer_{h5Writer}
    {
        if constexpr (dimension < 2)
            throw std::runtime_error("VTK Diagnostics do not support 1D");
    }

protected:
    auto static flatten_boxes(auto const& boxes)
    {
        std::vector<std::array<int, dimension * 2>> data(boxes.size());
        std::size_t pos = 0;
        for (auto const& box : boxes)
        {
            for (std::uint16_t i = 0; i < dimension; ++i)
            {
                data[pos][0 + i * 2] = box.lower[i];
                data[pos][1 + i * 2] = box.upper[i];
            }
            ++pos;
        }
        return data;
    }


    // main entry point from quantity specific writers
    struct VTKFileWriter;

    // Mostly gridlayout and box info
    struct VTKFileFieldInfo;

    // Writes fields to HDF5 in VTK format
    struct VTKFileFieldWriter;

    // Writes tensor fields to HDF5 in VTK format
    template<std::size_t rank>
    struct VTKFileTensorFieldWriter;

    auto& getOrCreateH5File(DiagnosticProperties const& diagnostic)
    {
        if (!fileData_.count(diagnostic.quantity))
            fileData_.emplace(diagnostic.quantity, this->h5Writer_.makeFile(diagnostic));
        return *fileData_.at(diagnostic.quantity);
    }


    Writer& h5Writer_;
    std::unordered_map<std::string, std::unique_ptr<HighFiveFile>> fileData_;
};



template<typename Writer>
struct H5TypeWriter<Writer>::VTKFileFieldInfo
{
    auto static constexpr primal_qty = physical_quantity_type::rho;

    VTKFileFieldInfo(auto const& lyout)
        : lvl{std::to_string(lyout.levelNumber())}
        , layout{lyout}
        , ghost_box{layout.AMRGhostBoxFor(primal_qty)}
        , local_box{layout.AMRToLocal(core::grow(ghost_box, -1 * layout.nbrGhosts()))}
    {
    }


    std::string lvl;
    GridLayout const& layout;
    Box_t const ghost_box;
    core::Box<std::uint32_t, dimension> const local_box;
    std::uint32_t const primal_row_len = local_box.shape(dimension - 1);
    std::string const path             = level_base + lvl + "/PointData/data";
};



template<typename Writer>
struct H5TypeWriter<Writer>::VTKFileFieldWriter
{
    void write2D(auto const& field)
    {
        auto ds            = fw->h5file.getDataSet(finfo.path);
        auto const lcl_box = finfo.local_box;
        auto const write   = [&]() {
            for (std::uint32_t i = lcl_box.lower[0]; i <= lcl_box.upper[0];
                 ++i, data_offset += finfo.primal_row_len)
                ds.select({data_offset}, {finfo.primal_row_len})
                    .write_raw(&field(i, lcl_box.lower[1]));
        };
        write();
        write();
    }


    void write3D(auto const& field)
    {
        auto ds            = fw->h5file.getDataSet(finfo.path);
        auto const lcl_box = finfo.local_box;
        for (std::uint32_t i = lcl_box.lower[0]; i <= lcl_box.upper[0]; ++i)
            for (std::uint32_t j = lcl_box.lower[1]; j <= lcl_box.upper[1];
                 ++j, data_offset += finfo.primal_row_len)
                ds.select({data_offset}, {finfo.primal_row_len})
                    .write_raw(&field(i, j, lcl_box.lower[2]));
    }


    void operator()(auto const& field)
    {
        auto& tmp         = fw->typewriter->h5Writer_.modelView().tmpField();
        auto const frimal = core::convert_to_fortran_primal(tmp, field, finfo.layout);
        if constexpr (dimension == 2)
            write2D(frimal);
        if constexpr (dimension == 3)
            write3D(frimal);
    }

    VTKFileWriter* fw;
    VTKFileFieldInfo finfo;
    std::size_t& data_offset = fw->data_offset; // NOT OWNED HERE
};

template<typename Writer>
template<std::size_t rank>
struct H5TypeWriter<Writer>::VTKFileTensorFieldWriter
{
    auto static constexpr N = core::detail::tensor_field_dim_from_rank<rank>();

    void write2D(auto const& tf)
    {
        auto ds            = fw->h5file.getDataSet(finfo.path);
        auto const lcl_box = finfo.local_box;
        auto const write   = [&]() {
            for (std::uint32_t i = lcl_box.lower[0]; i <= lcl_box.upper[0];
                 ++i, data_offset += finfo.primal_row_len)
                for (std::uint32_t c = 0; c < N; ++c)
                    ds.select({data_offset, c}, {finfo.primal_row_len, 1})
                        .write_raw(&tf[c](i, lcl_box.lower[1]));
        };
        write();
        write();
    }


    void write3D(auto const& tf)
    {
        auto ds            = fw->h5file.getDataSet(finfo.path);
        auto const lcl_box = finfo.local_box;
        for (std::uint32_t i = lcl_box.lower[0]; i <= lcl_box.upper[0]; ++i)
            for (std::uint32_t j = lcl_box.lower[1]; j <= lcl_box.upper[1];
                 ++j, data_offset += finfo.primal_row_len)
                for (std::uint32_t c = 0; c < N; ++c)
                    ds.select({data_offset, c}, {finfo.primal_row_len, 1})
                        .write_raw(&tf[c](i, j, lcl_box.lower[2]));
    }


    void operator()(auto const& tf)
    {
        auto& tmp         = fw->typewriter->h5Writer_.modelView().template tmpTensorField<rank>();
        auto const frimal = core::convert_to_fortran_primal(tmp, tf, finfo.layout);
        if constexpr (dimension == 2)
            write2D(frimal);
        if constexpr (dimension == 3)
            write3D(frimal);
    }

    VTKFileWriter* fw;
    VTKFileFieldInfo finfo;
    std::size_t& data_offset = fw->data_offset; // NOT OWNED HERE
};



template<typename Writer>
struct H5TypeWriter<Writer>::VTKFileWriter
{
    std::size_t static constexpr boxValsIn3D = 6; // lo0, up0, lo1, up1, lo2, up2

    VTKFileWriter(auto& prop, auto tw)
        : diagnostic{prop}
        , typewriter{tw}
        , h5file{tw->getOrCreateH5File(prop)}
    {
        {
            initDSDefault(base + "/Steps/Values");
            auto steps_group = h5file.file().getGroup(base + "/Steps");
            if (!steps_group.hasAttribute("NSteps"))
                steps_group.template createAttribute<int>("NSteps", HighFive::DataSpace::From(0))
                    .write(0);
            auto steps_attr = steps_group.getAttribute("NSteps");
            steps_attr.write(steps_attr.template read<int>() + 1);
        }

        {
            auto const& timestamp = tw->h5Writer_.timestamp();
            auto ds               = h5file.getDataSet(base + "/Steps/Values");
            auto const old_size   = ds.getDimensions()[0];
            ds.resize({old_size + 1});
            ds.select({old_size}, {1}).write(timestamp);
        }

        auto root = h5file.file().getGroup(base);

        if (!root.hasAttribute("Version"))
            root.template createAttribute<std::array<int, 2>>("Version", std::array<int, 2>{2, 2});

        if (!root.hasAttribute("Origin"))
            root.template createAttribute<std::array<int, 3>>("Origin",
                                                              std::array<int, 3>{0, 0, 0});

        if (!root.hasAttribute("Type"))
            root.template createAttribute<std::string>("Type", "OverlappingAMR");

        if (!root.hasAttribute("GridDescription"))
            root.template createAttribute<std::string>("GridDescription", "XYZ");
    }


    auto level_spacing(std::size_t const lvl) const
    {
        auto const mesh_size = typewriter->h5Writer_.modelView().cellWidth();
        return core::for_N_make_array<dimension>(
            [&](auto i) { return static_cast<float>(mesh_size[i] / std::pow(2, lvl)); });
    }

    void writeField(auto const& field, auto const& layout)
    {
        VTKFileFieldWriter{this, {layout}}(field);
    }

    template<std::size_t rank = 2> // TODO convert duals to primals
    void writeTensorField(auto const& tf, auto const& layout)
    {
        VTKFileTensorFieldWriter<rank>{this, {layout}}(tf);
    }

    template<typename T = FloatType>
    void initDS(auto const& path, auto const& ds) const
    {
        h5file.template create_chunked_data_set<T>(path, std::vector<hsize_t>{detail::CHUNK_SIZE},
                                                   ds);
    }

    template<typename T = FloatType>
    void initDSDefault(auto const& path) const
    {
        initDS<T>(path, HighFive::DataSpace({0}, {HighFive::DataSpace::UNLIMITED}));
    }

    void initFileLevel(int const level) const
    {
        auto const lvl = std::to_string(level);

        {
            auto const path = level_base + lvl + "/AMRBox";
            h5file.template create_chunked_data_set<int>(
                path, std::vector<hsize_t>{detail::CHUNK_SIZE, boxValsIn3D},
                HighFive::DataSpace({0, boxValsIn3D},
                                    {HighFive::DataSpace::UNLIMITED, boxValsIn3D}));
        }

        auto const chucked_int_datasets
            = std::array{step_level + lvl + "/PointDataOffset/data",
                         step_level + lvl + "/AMRBoxOffset", step_level + lvl + "/NumberOfAMRBox"};

        for (auto const& path : chucked_int_datasets)
        {
            h5file.template create_chunked_data_set<int>(
                path, std::vector<hsize_t>{detail::CHUNK_SIZE},
                HighFive::DataSpace({0}, {HighFive::DataSpace::UNLIMITED}));
        }

        auto level_group = h5file.file().getGroup(level_base + lvl);
        if (!level_group.hasAttribute("Spacing"))
        {
            level_group.template createAttribute<std::array<float, dimension>>(
                "Spacing", level_spacing(level));

            level_group.createGroup("CellData");
            level_group.createGroup("FieldData");

            auto steps_group = h5file.file().getGroup(step_level + lvl);
            steps_group.createGroup("CellDataOffset");
            steps_group.createGroup("FieldDataOffset");
        }
    }

    void initFieldFileLevel(int const level, auto& boxes)
    {
        initDSDefault(level_base + std::to_string(level) + "/PointData/data");
        resize(level, boxes);
    }

    template<std::size_t rank = 2>
    void initTensorFieldFileLevel(auto const level, auto& boxes)
    {
        auto constexpr N = core::detail::tensor_field_dim_from_rank<rank>();
        auto const path  = level_base + std::to_string(level) + "/PointData/data";
        h5file.template create_chunked_data_set<FloatType>(
            path, std::vector<hsize_t>{detail::CHUNK_SIZE, N},
            HighFive::DataSpace({0, N}, {HighFive::DataSpace::UNLIMITED, N}));
        resize<N>(level, boxes);
    }



    template<std::size_t N = 1>
    void resize_data(auto const ilvl, auto const& boxes)
    {
        // data is per face of a cube (cell)
        //  so 2d data is * 2
        //  and 1d data is * 4
        constexpr static auto X_TIMES = std::array{4, 2, /* 3d noop */ 1}[dimension - 1];

        auto const lvl       = std::to_string(ilvl);
        auto const data_path = level_base + lvl + "/PointData/data";
        auto point_data_ds   = h5file.getDataSet(data_path);
        data_offset          = point_data_ds.getDimensions()[0];

        {
            auto ds             = h5file.getDataSet(step_level + lvl + "/PointDataOffset/data");
            auto const old_size = ds.getDimensions()[0];
            ds.resize({old_size + 1});
            ds.select({old_size}, {1}).write(data_offset);
        }

        auto const rank_data_size = core::mpi::collect(
            core::sum_from(boxes, [](auto const& b) { return b.size() * X_TIMES; }));
        auto const new_size = data_offset + core::sum(rank_data_size);

        if constexpr (N == 1)
            point_data_ds.resize({new_size});
        else
            point_data_ds.resize({new_size, N});
        for (int i = 0; i < core::mpi::rank(); ++i)
            data_offset += rank_data_size[i];
    }


    void resize_boxes(auto const ilvl, auto const& boxes)
    {
        auto const lvl           = std::to_string(ilvl);
        auto const rank_box_size = core::mpi::collect(boxes.size());
        auto const total_boxes   = core::sum(rank_box_size);

        {
            auto ds             = h5file.getDataSet(step_level + lvl + "/NumberOfAMRBox");
            auto const old_size = ds.getDimensions()[0];
            ds.resize({old_size + 1});
            ds.select({old_size}, {1}).write(total_boxes);
        }

        auto amrbox_ds = h5file.getDataSet(level_base + lvl + "/AMRBox");
        box_offset     = amrbox_ds.getDimensions()[0];

        {
            auto ds             = h5file.getDataSet(step_level + lvl + "/AMRBoxOffset");
            auto const old_size = ds.getDimensions()[0];
            ds.resize({old_size + 1});
            ds.select({old_size}, {1}).write(box_offset);
        }

        amrbox_ds.resize({box_offset + total_boxes, boxValsIn3D});
        for (int i = 0; i < core::mpi::rank(); ++i)
            box_offset += rank_box_size[i];

        auto const vtk_boxes = flatten_boxes(boxes);
        amrbox_ds.select({box_offset, 0}, {boxes.size(), dimension * 2}).write(vtk_boxes);
    }

    template<std::size_t N = 1>
    void resize(auto const ilvl, auto& boxes)
    {
        resize_boxes(ilvl, boxes);
        for (auto& box : boxes)
            box.upper += 1; // primal
        resize_data<N>(ilvl, boxes);
    }


    DiagnosticProperties& diagnostic;
    H5TypeWriter<Writer>* typewriter;
    HighFiveFile& h5file;
    std::size_t box_offset = 0, data_offset = 0;
};


} // namespace PHARE::diagnostic::vtkh5

#endif // PHARE_DIAGNOSTIC_DETAIL_VTK_H5_TYPE_WRITER_HPP
