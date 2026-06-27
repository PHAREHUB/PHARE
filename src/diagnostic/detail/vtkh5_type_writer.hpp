#ifndef PHARE_DIAGNOSTIC_DETAIL_VTK_H5_TYPE_WRITER_HPP
#define PHARE_DIAGNOSTIC_DETAIL_VTK_H5_TYPE_WRITER_HPP


#include "core/logger.hpp"
#include "core/utilities/types.hpp"
#include "core/utilities/algorithm.hpp"
#include "core/utilities/mpi_utils.hpp"
#include "core/data/tensorfield/tensorfield.hpp"

#include "amr/resources_manager/amr_utils.hpp"

#include "diagnostic/diagnostic_writer.hpp"

#include "hdf5/detail/h5/h5_file.hpp"

#include <string>
#include <unordered_map>

#if !defined(PHARE_DIAG_DOUBLES)
#error // PHARE_DIAG_DOUBLES not defined
#endif


namespace PHARE::diagnostic::vtkh5
{

using namespace hdf5::h5;


template<std::size_t dim>
struct HierarchyData
{
    // data is duplicated in lower dimensions to fit a 3d view
    //  so 2d data is * 2
    //  and 1d data is * 4
    constexpr static auto X_TIMES = std::array{4, 2, 1}[dim - 1];

    static auto& INSTANCE()
    {
        static HierarchyData data;
        return data;
    }

    static void reset(auto& h5Writer) // called once per dump
    {
        PHARE_LOG_SCOPE(3, "HierarchyData<H5Writer::reset");

        auto const maxLevels = h5Writer.modelView().maxLevel() + 1;
        auto& data           = INSTANCE();

        data.level_data_size.resize(maxLevels);
        data.n_boxes_per_level.resize(maxLevels);
        data.level_boxes_per_rank.resize(maxLevels);
        data.level_rank_data_size.resize(maxLevels);
        data.flattened_lcl_level_boxes.resize(maxLevels);

        h5Writer.modelView().onLevels(
            [&](auto const& level) {
                auto const ilvl = level.getLevelNumber();

                data.n_boxes_per_level[ilvl]    = 0;
                data.level_boxes_per_rank[ilvl] = amr::boxesPerRankOn<dim>(level);
                data.flattened_lcl_level_boxes[ilvl]
                    = flatten_boxes(data.level_boxes_per_rank[ilvl][core::mpi::rank()]);

                for (std::size_t i = 0; i < data.level_boxes_per_rank[ilvl].size(); ++i)
                {
                    data.level_rank_data_size[ilvl].resize(core::mpi::size());

                    data.level_rank_data_size[ilvl][i] = 0;
                    for (auto box : data.level_boxes_per_rank[ilvl][i])
                    {
                        box.upper += 1;                                             // all primal
                        data.level_rank_data_size[ilvl][i] += box.size() * X_TIMES; // no ghosts
                    }
                    data.n_boxes_per_level[ilvl] += data.level_boxes_per_rank[ilvl][i].size();
                }

                data.level_data_size[ilvl] = core::sum(data.level_rank_data_size[ilvl]);
            },
            [&](int const ilvl) { // ilvl does not exist currently
                data.level_data_size[ilvl]   = 0;
                data.n_boxes_per_level[ilvl] = 0;
                data.level_rank_data_size[ilvl].clear();
                data.level_rank_data_size[ilvl].resize(core::mpi::size());
                data.level_boxes_per_rank[ilvl].clear();
                data.level_boxes_per_rank[ilvl].resize(core::mpi::size());
                data.flattened_lcl_level_boxes[ilvl].clear();
            },
            h5Writer.minLevel, h5Writer.maxLevel);
    }


    auto static flatten_boxes(auto const& boxes)
    {
        std::vector<std::array<int, dim * 2>> data(boxes.size());
        std::size_t pos = 0;
        for (auto const& box : boxes)
        {
            for (std::size_t i = 0; i < dim; ++i)
            {
                data[pos][0 + i * 2] = box.lower[i];
                data[pos][1 + i * 2] = box.upper[i];
            }
            ++pos;
        }
        return data;
    }

    std::vector<std::vector<std::array<int, dim * 2>>> flattened_lcl_level_boxes;
    std::vector<std::vector<std::vector<core::Box<int, dim>>>> level_boxes_per_rank;
    std::vector<std::vector<int>> level_rank_data_size;
    std::vector<int> n_boxes_per_level, level_data_size;
};


template<typename Writer>
class H5TypeWriter : public PHARE::diagnostic::TypeWriter
{
    using FloatType                      = std::conditional_t<PHARE_DIAG_DOUBLES, double, float>;
    using HierData                       = HierarchyData<Writer::dimension>;
    using ModelView                      = Writer::ModelView;
    using physical_quantity_type         = ModelView::physical_quantity_type;
    using Attributes                     = Writer::Attributes;
    std::string static inline const base = "/VTKHDF";
    std::string static inline const level_base = base + "/Level";
    std::string static inline const step_level = base + "/Steps/Level";
    auto static inline const level_data_path
        = [](auto const ilvl) { return level_base + std::to_string(ilvl) + "/PointData/data"; };


public:
    static constexpr auto dimension = Writer::dimension;

    H5TypeWriter(Writer& h5Writer)
        : h5Writer_{h5Writer}
    {
    }

    virtual void setup(DiagnosticProperties&) = 0;

    void writeFileAttributes(DiagnosticProperties const&, Attributes const&);

    //------ valid for all h5type writers -------------------------------------
    void finalize(DiagnosticProperties& diagnostic)
    {
        // we close the file by removing the associated file
        // from the map. This is done only at flush time otherwise
        ++diagnostic.dumpIdx;

        assert(diagnostic.params.contains("flush_every"));

        std::size_t flushEvery = diagnostic.param<std::size_t>("flush_every");

        if (flushEvery != Writer::flush_never and diagnostic.dumpIdx % flushEvery == 0)
        {
            fileData_.erase(diagnostic.fileKey());
            assert(fileData_.count(diagnostic.fileKey()) == 0);
        }
    }

protected:
    class VTKFileInitializer;

    class VTKFileWriter;


    auto& getOrCreateH5File(DiagnosticProperties const& diagnostic)
    {
        if (!fileExistsFor(diagnostic))
            fileData_.emplace(diagnostic.fileKey(), this->h5Writer_.makeFile(diagnostic));
        return *fileData_.at(diagnostic.fileKey());
    }

    bool fileExistsFor(DiagnosticProperties const& diagnostic) const
    {
        return fileData_.count(diagnostic.fileKey());
    }


    Writer& h5Writer_;
    std::unordered_map<std::string, std::unique_ptr<HighFiveFile>> fileData_;
};


template<typename Writer>
void H5TypeWriter<Writer>::writeFileAttributes(DiagnosticProperties const& prop,
                                               Attributes const& attrs)
{
    if (fileExistsFor(prop)) // otherwise qty not supported (yet)
    {
        h5Writer_.writeGlobalAttributeDict(*fileData_.at(prop.fileKey()), attrs, "/");
        if (prop.nAttributes > 0)
            h5Writer_.writeGlobalAttributeDict(*fileData_.at(prop.fileKey()),
                                               prop.fileAttributes, "/");
    }
}


template<typename Writer>
class H5TypeWriter<Writer>::VTKFileInitializer
{
    auto static constexpr primal_qty         = physical_quantity_type::all_primal_field;
    std::size_t static constexpr boxValsIn3D = 6; // lo0, up0, lo1, up1, lo2, up2
public:
    VTKFileInitializer(DiagnosticProperties const& prop, H5TypeWriter<Writer>* const typewriter);

    void initFileLevel(int const ilvl);

    std::size_t initFieldFileLevel(auto const ilvl) { return initAnyFieldLevel(ilvl); }

    template<std::size_t rank = 2>
    std::size_t initTensorFieldFileLevel(auto const ilvl)
    {
        return initAnyFieldLevel<core::detail::tensor_field_dim_from_rank<rank>()>(ilvl);
    }

    template<std::size_t rank = 2>
    std::size_t initTensorFieldFileLevelWithSlice(auto const ilvl,
                                                   auto const& amr_slice_box)
    {
        return initAnyFieldLevelSliced<core::detail::tensor_field_dim_from_rank<rank>()>(
            ilvl, amr_slice_box);
    }

private:
    auto level_spacing(int const lvl) const
    {
        auto const mesh_size = typewriter->h5Writer_.modelView().cellWidth();
        return core::for_N_make_array<dimension>(
            [&](auto i) { return static_cast<float>(mesh_size[i] / std::pow(2, lvl)); });
    }

    template<std::size_t N = 1>
    std::size_t initAnyFieldLevel(auto const ilvl)
    {
        h5file.create_resizable_2d_data_set<FloatType, N>(level_data_path(ilvl));
        initFileLevel(ilvl);
        resize_boxes(ilvl);
        resize_data<N>(ilvl);
        return data_offset;
    }

    template<std::size_t N = 1>
    std::size_t initAnyFieldLevelSliced(auto const ilvl,
                                         auto const& amr_slice_box)
    {
        h5file.create_resizable_2d_data_set<FloatType, N>(level_data_path(ilvl));
        initFileLevel(ilvl);
        resize_boxes_sliced(ilvl, amr_slice_box);
        resize_data_sliced<N>(ilvl, amr_slice_box);
        return data_offset;
    }

    template<std::size_t N = 1>
    void resize_data(int const ilvl);

    template<std::size_t N = 1>
    void resize_data_sliced(int const ilvl, auto const& amr_slice_box);

    void resize_boxes(int const ilvl);

    void resize_boxes_sliced(int const ilvl, auto const& amr_slice_box);

    bool const newFile;
    DiagnosticProperties const& diagnostic;
    H5TypeWriter<Writer>* const typewriter;
    HighFiveFile& h5file;
    std::size_t data_offset;
    ModelView& modelView = typewriter->h5Writer_.modelView();
};


template<typename Writer>
class H5TypeWriter<Writer>::VTKFileWriter
{
    auto static constexpr primal_qty         = physical_quantity_type::all_primal_field;
    std::size_t static constexpr boxValsIn3D = 6; // lo0, up0, lo1, up1, lo2, up2


public:
    VTKFileWriter(DiagnosticProperties const& prop, H5TypeWriter<Writer>* const typewriter,
                  std::size_t const offset);

    void writeField(auto const& field, auto const& layout);

    template<std::size_t rank = 2>
    void writeTensorField(auto const& tf, auto const& layout);

    template<std::size_t rank = 2>
    void writeTensorFieldSlice(auto const& tf, auto const& layout,
                                auto const& amr_slice_box);



private:
    auto local_box(auto const& layout) const
    {
        return layout.AMRToLocal(layout.AMRBoxFor(primal_qty));
    }

    DiagnosticProperties const& diagnostic;
    H5TypeWriter<Writer>* const typewriter;
    HighFiveFile& h5file;
    std::size_t data_offset = 0;
    ModelView& modelView    = typewriter->h5Writer_.modelView();
};


template<typename Writer>
H5TypeWriter<Writer>::VTKFileInitializer::VTKFileInitializer(DiagnosticProperties const& prop,
                                                             H5TypeWriter<Writer>* const tw)
    : newFile{!tw->fileExistsFor(prop)}
    , diagnostic{prop}
    , typewriter{tw}
    , h5file{tw->getOrCreateH5File(prop)}
{
    {
        PHARE_LOG_SCOPE(3, "VTKFileInitializer::VTKFileInitializer::0");
        h5file.create_resizable_1d_data_set<FloatType>(base + "/Steps/Values");
        auto steps_group = h5file.file().getGroup(base + "/Steps");
        if (!steps_group.hasAttribute("NSteps"))
            steps_group.template createAttribute<int>("NSteps", HighFive::DataSpace::From(0))
                .write(0);
        auto steps_attr = steps_group.getAttribute("NSteps");
        steps_attr.write(steps_attr.template read<int>() + 1);
    }

    {
        PHARE_LOG_SCOPE(3, "VTKFileInitializer::VTKFileInitializer::1");
        auto const& timestamp = typewriter->h5Writer_.timestamp();
        auto ds               = h5file.getDataSet(base + "/Steps/Values");
        auto const old_size   = ds.getDimensions()[0];
        ds.resize({old_size + 1});
        ds.select({old_size}, {1}).write(timestamp);
    }

    if (newFile)
    {
        PHARE_LOG_SCOPE(3, "VTKFileInitializer::VTKFileInitializer::2");
        auto root = h5file.file().getGroup(base);

        if (!root.hasAttribute("Version"))
            root.template createAttribute<std::array<int, 2>>("Version", std::array<int, 2>{2, 2});

        if (!root.hasAttribute("Origin"))
            root.template createAttribute<std::array<float, 3>>("Origin",
                                                                std::array<float, 3>{0, 0, 0});

        if (!root.hasAttribute("Type"))
            root.template createAttribute<std::string>("Type", "OverlappingAMR");

        if (!root.hasAttribute("GridDescription"))
            root.template createAttribute<std::string>("GridDescription", "XYZ");
    }
}


template<typename Writer>
H5TypeWriter<Writer>::VTKFileWriter::VTKFileWriter(DiagnosticProperties const& prop,
                                                   H5TypeWriter<Writer>* const tw,
                                                   std::size_t const offset)
    : diagnostic{prop}
    , typewriter{tw}
    , h5file{tw->getOrCreateH5File(prop)}
    , data_offset{offset}
{
}


template<typename Writer>
void H5TypeWriter<Writer>::VTKFileWriter::writeField(auto const& field, auto const& layout)
{
    PHARE_LOG_SCOPE(3, "VTKFileWriter::writeField");

    auto const& frimal = core::convert_to_fortran_primal(modelView.tmpField(), field, layout);
    auto const size    = local_box(layout).size();
    auto ds            = h5file.getDataSet(level_data_path(layout.levelNumber()));

    PHARE_LOG_SCOPE(3, "VTKFileWriter::writeField::0");
    for (std::uint16_t i = 0; i < HierData::X_TIMES; ++i)
    {
        ds.select({data_offset, 0}, {size, 1}).write_raw(frimal.data());
        data_offset += size;
    }
}


template<typename Writer>
template<std::size_t rank>
void H5TypeWriter<Writer>::VTKFileWriter::writeTensorField(auto const& tf, auto const& layout)
{
    PHARE_LOG_SCOPE(3, "VTKFileWriter::writeTensorField");

    auto static constexpr N = core::detail::tensor_field_dim_from_rank<rank>();
    auto const& frimal
        = core::convert_to_fortran_primal(modelView.template tmpTensorField<rank>(), tf, layout);
    auto const size = local_box(layout).size();
    auto ds         = h5file.getDataSet(level_data_path(layout.levelNumber()));

    PHARE_LOG_SCOPE(3, "VTKFileWriter::writeTensorField::0");
    for (std::uint16_t i = 0; i < HierData::X_TIMES; ++i)
    {
        for (std::uint32_t c = 0; c < N; ++c)
            ds.select({data_offset, c}, {size, 1}).write_raw(frimal[c].data());
        data_offset += size;
    }
}


template<typename Writer>
void H5TypeWriter<Writer>::VTKFileInitializer::initFileLevel(int const ilvl)
{
    PHARE_LOG_SCOPE(3, "VTKFileInitializer::initFileLevel");

    auto const lvl = std::to_string(ilvl);

    h5file.create_resizable_2d_data_set<int, boxValsIn3D>(level_base + lvl + "/AMRBox");
    h5file.create_resizable_1d_data_set<int>(step_level + lvl + "/AMRBoxOffset");
    h5file.create_resizable_1d_data_set<int>(step_level + lvl + "/NumberOfAMRBox");
    h5file.create_resizable_1d_data_set<int>(step_level + lvl + "/PointDataOffset/data");

    auto level_group = h5file.file().getGroup(level_base + lvl);
    if (!level_group.hasAttribute("Spacing"))
    {
        PHARE_LOG_SCOPE(3, "VTKFileInitializer::initFileLevel::0");

        level_group.template createAttribute<std::array<float, dimension>>("Spacing",
                                                                           level_spacing(ilvl));

        level_group.createGroup("CellData");
        level_group.createGroup("FieldData");

        auto steps_group = h5file.file().getGroup(step_level + lvl);
        steps_group.createGroup("CellDataOffset");
        steps_group.createGroup("FieldDataOffset");
    }
}


template<typename Writer>
template<std::size_t N>
void H5TypeWriter<Writer>::VTKFileInitializer::resize_data(int const ilvl)
{
    PHARE_LOG_SCOPE(3, "VTKFileInitializer::resize_data");

    auto const& hier_data = HierData::INSTANCE();
    auto const lvl        = std::to_string(ilvl);
    auto const data_path  = level_data_path(ilvl);
    auto point_data_ds    = h5file.getDataSet(data_path);
    data_offset           = point_data_ds.getDimensions()[0];

    {
        PHARE_LOG_SCOPE(3, "VTKFileInitializer::resize_data::0");
        auto ds             = h5file.getDataSet(step_level + lvl + "/PointDataOffset/data");
        auto const old_size = ds.getDimensions()[0];
        ds.resize({old_size + 1});
        ds.select({old_size}, {1}).write(data_offset);
    }

    PHARE_LOG_SCOPE(3, "VTKFileInitializer::resize_data::1");
    auto const& rank_data_sizes = hier_data.level_rank_data_size[ilvl];
    auto const& level_data_size = hier_data.level_data_size[ilvl];
    auto const new_size         = data_offset + level_data_size;

    point_data_ds.resize({new_size, N});
    for (int i = 0; i < core::mpi::rank(); ++i)
        data_offset += rank_data_sizes[i];
}


template<typename Writer>
void H5TypeWriter<Writer>::VTKFileInitializer::resize_boxes(int const ilvl)
{
    PHARE_LOG_SCOPE(3, "VTKFileInitializer::resize_boxes");

    auto const lvl          = std::to_string(ilvl);
    auto const& hier_data   = HierData::INSTANCE();
    auto const& rank_boxes  = hier_data.level_boxes_per_rank[ilvl];
    auto const& total_boxes = hier_data.n_boxes_per_level[ilvl];

    {
        PHARE_LOG_SCOPE(3, "VTKFileInitializer::resize_boxes::0");
        auto ds             = h5file.getDataSet(step_level + lvl + "/NumberOfAMRBox");
        auto const old_size = ds.getDimensions()[0];
        ds.resize({old_size + 1});
        ds.select({old_size}, {1}).write(total_boxes);
    }

    auto amrbox_ds  = h5file.getDataSet(level_base + lvl + "/AMRBox");
    auto box_offset = amrbox_ds.getDimensions()[0];

    {
        PHARE_LOG_SCOPE(3, "VTKFileInitializer::resize_boxes::1");
        auto ds             = h5file.getDataSet(step_level + lvl + "/AMRBoxOffset");
        auto const old_size = ds.getDimensions()[0];
        ds.resize({old_size + 1});
        ds.select({old_size}, {1}).write(box_offset);
    }

    PHARE_LOG_SCOPE(3, "VTKFileInitializer::resize_boxes::2");
    amrbox_ds.resize({box_offset + total_boxes, boxValsIn3D});
    for (int i = 0; i < core::mpi::rank(); ++i)
        box_offset += rank_boxes[i].size();


    PHARE_LOG_SCOPE(3, "VTKFileInitializer::resize_boxes::3");
    amrbox_ds.select({box_offset, 0}, {rank_boxes[core::mpi::rank()].size(), dimension * 2})
        .write(hier_data.flattened_lcl_level_boxes[ilvl]);
}


template<typename Writer>
template<std::size_t N>
void H5TypeWriter<Writer>::VTKFileInitializer::resize_data_sliced(
    int const ilvl, auto const& amr_slice_box)
{
    PHARE_LOG_SCOPE(3, "VTKFileInitializer::resize_data_sliced");

    auto const& hier_data  = HierData::INSTANCE();
    auto const lvl         = std::to_string(ilvl);
    auto const data_path   = level_data_path(ilvl);
    auto point_data_ds     = h5file.getDataSet(data_path);
    data_offset            = point_data_ds.getDimensions()[0];

    {
        PHARE_LOG_SCOPE(3, "VTKFileInitializer::resize_data_sliced::0");
        auto ds             = h5file.getDataSet(step_level + lvl + "/PointDataOffset/data");
        auto const old_size = ds.getDimensions()[0];
        ds.resize({old_size + 1});
        ds.select({old_size}, {1}).write(data_offset);
    }

    auto const& rank_boxes = hier_data.level_boxes_per_rank[ilvl];
    std::vector<int> rank_data_sizes(core::mpi::size(), 0);
    for (int r = 0; r < core::mpi::size(); ++r)
        for (auto box : rank_boxes[r])
        {
            box.upper += 1; // all primal
            if (auto const isect = box * amr_slice_box)
            {
                // AMRBox uses cell indices; a degenerate primal dim (lo==hi) maps to 1 cell = 2 pts
                std::size_t primal_count = 1;
                for (std::size_t i = 0; i < dimension; ++i)
                    primal_count *= static_cast<std::size_t>(
                        std::max(isect->lower[i] + 1, isect->upper[i]) - isect->lower[i] + 1);
                rank_data_sizes[r] += static_cast<int>(primal_count) * HierData::X_TIMES;
            }
        }

    PHARE_LOG_SCOPE(3, "VTKFileInitializer::resize_data_sliced::1");
    auto const level_data_size = core::sum(rank_data_sizes);
    auto const new_size        = data_offset + static_cast<std::size_t>(level_data_size);
    point_data_ds.resize({new_size, N});
    for (int r = 0; r < core::mpi::rank(); ++r)
        data_offset += static_cast<std::size_t>(rank_data_sizes[r]);
}


template<typename Writer>
void H5TypeWriter<Writer>::VTKFileInitializer::resize_boxes_sliced(
    int const ilvl, auto const& amr_slice_box)
{
    PHARE_LOG_SCOPE(3, "VTKFileInitializer::resize_boxes_sliced");

    constexpr auto dim     = Writer::dimension;
    auto const& hier_data  = HierData::INSTANCE();
    auto const lvl         = std::to_string(ilvl);
    auto const& rank_boxes = hier_data.level_boxes_per_rank[ilvl];

    auto primal_isect = [&](auto box) { return (box.upper += 1, box) * amr_slice_box; };

    std::vector<int> n_isect_per_rank(core::mpi::size(), 0);
    for (int r = 0; r < core::mpi::size(); ++r)
        for (auto const& box : rank_boxes[r])
            if (primal_isect(box))
                ++n_isect_per_rank[r];
    auto const total_isect = core::sum(n_isect_per_rank);

    {
        PHARE_LOG_SCOPE(3, "VTKFileInitializer::resize_boxes_sliced::0");
        auto ds             = h5file.getDataSet(step_level + lvl + "/NumberOfAMRBox");
        auto const old_size = ds.getDimensions()[0];
        ds.resize({old_size + 1});
        ds.select({old_size}, {1}).write(total_isect);
    }

    auto amrbox_ds  = h5file.getDataSet(level_base + lvl + "/AMRBox");
    auto box_offset = static_cast<int>(amrbox_ds.getDimensions()[0]);

    {
        PHARE_LOG_SCOPE(3, "VTKFileInitializer::resize_boxes_sliced::1");
        auto ds             = h5file.getDataSet(step_level + lvl + "/AMRBoxOffset");
        auto const old_size = ds.getDimensions()[0];
        ds.resize({old_size + 1});
        ds.select({old_size}, {1}).write(box_offset);
    }

    PHARE_LOG_SCOPE(3, "VTKFileInitializer::resize_boxes_sliced::2");
    amrbox_ds.resize({static_cast<std::size_t>(box_offset + total_isect), boxValsIn3D});
    for (int r = 0; r < core::mpi::rank(); ++r)
        box_offset += n_isect_per_rank[r];

    PHARE_LOG_SCOPE(3, "VTKFileInitializer::resize_boxes_sliced::3");
    std::vector<std::array<int, dim * 2>> isect_boxes;
    for (auto const& box : rank_boxes[core::mpi::rank()])
        if (auto const isect = primal_isect(box))
        {
            std::array<int, dim * 2> arr{};
            for (std::size_t i = 0; i < dim; ++i)
            {
                // AMRBox uses cell indices: cell_hi = primal_hi - 1
                // For degenerate dims (primal lo==hi) clamp to lo so cell box has 1 cell
                arr[i * 2 + 0] = (*isect).lower[i];
                arr[i * 2 + 1] = std::max((*isect).lower[i], (*isect).upper[i] - 1);
            }
            isect_boxes.push_back(arr);
        }

    if (!isect_boxes.empty())
        amrbox_ds
            .select({static_cast<std::size_t>(box_offset), 0}, {isect_boxes.size(), dim * 2})
            .write(isect_boxes);
}


template<typename Writer>
template<std::size_t rank>
void H5TypeWriter<Writer>::VTKFileWriter::writeTensorFieldSlice(
    auto const& tf, auto const& layout, auto const& amr_slice_box)
{
    PHARE_LOG_SCOPE(3, "VTKFileWriter::writeTensorFieldSlice");

    constexpr auto dim      = Writer::dimension;
    auto static constexpr N = core::detail::tensor_field_dim_from_rank<rank>();

    auto const amr_patch_box = layout.AMRBoxFor(primal_qty);
    auto const maybe_isect   = amr_patch_box * amr_slice_box;
    if (!maybe_isect)
        return;

    auto const& frimal
        = core::convert_to_fortran_primal(modelView.template tmpTensorField<rank>(), tf, layout);

    auto const lcl_full  = layout.AMRToLocal(amr_patch_box);
    auto const lcl_isect = layout.AMRToLocal(*maybe_isect);
    auto const shape     = *lcl_full.shape();
    auto const zero_lo   = lcl_isect.lower - lcl_full.lower;
    // Cell-align: a degenerate primal dim (lo==hi) maps to 1 cell, so extend hi by 1
    auto zero_up = lcl_isect.upper - lcl_full.lower;
    for (std::size_t i = 0; i < dim; ++i)
        if (zero_up[i] == zero_lo[i])
            zero_up[i] = zero_lo[i] + 1;
    std::size_t size = 1;
    for (std::size_t i = 0; i < dim; ++i)
        size *= static_cast<std::size_t>(zero_up[i] - zero_lo[i] + 1);

    auto ds = h5file.getDataSet(level_data_path(layout.levelNumber()));

    PHARE_LOG_SCOPE(3, "VTKFileWriter::writeTensorFieldSlice::0");
    for (std::uint16_t i = 0; i < HierData::X_TIMES; ++i)
    {
        for (std::uint32_t c = 0; c < N; ++c)
        {
            auto const fview = core::make_array_view<false>(frimal[c].data(), shape);
            std::vector<FloatType> sliced(size);
            std::size_t idx = 0;

            if constexpr (dim == 1)
                for (auto ix = zero_lo[0]; ix <= zero_up[0]; ++ix)
                    sliced[idx++] = static_cast<FloatType>(fview(ix));
            else if constexpr (dim == 2)
                for (auto iy = zero_lo[1]; iy <= zero_up[1]; ++iy)
                    for (auto ix = zero_lo[0]; ix <= zero_up[0]; ++ix)
                        sliced[idx++] = static_cast<FloatType>(fview(ix, iy));
            else
                for (auto iz = zero_lo[2]; iz <= zero_up[2]; ++iz)
                    for (auto iy = zero_lo[1]; iy <= zero_up[1]; ++iy)
                        for (auto ix = zero_lo[0]; ix <= zero_up[0]; ++ix)
                            sliced[idx++] = static_cast<FloatType>(fview(ix, iy, iz));

            ds.select({data_offset, c}, {size, 1}).write_raw(sliced.data());
        }
        data_offset += size;
    }
}


} // namespace PHARE::diagnostic::vtkh5

#endif // PHARE_DIAGNOSTIC_DETAIL_VTK_H5_TYPE_WRITER_HPP
