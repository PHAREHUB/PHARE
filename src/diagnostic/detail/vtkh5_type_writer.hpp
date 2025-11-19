#ifndef PHARE_DIAGNOSTIC_DETAIL_VTK_H5_TYPE_WRITER_HPP
#define PHARE_DIAGNOSTIC_DETAIL_VTK_H5_TYPE_WRITER_HPP


#include "core/logger.hpp"
#include "core/utilities/types.hpp"
#include "core/utilities/box/box.hpp"
#include "core/utilities/algorithm.hpp"
#include "core/utilities/mpi_utils.hpp"
#include "core/data/tensorfield/tensorfield.hpp"

#include "diagnostic/diagnostic_writer.hpp"

#include "hdf5/detail/h5/h5_file.hpp"

#include <string>
#include <unordered_map>



namespace PHARE::diagnostic::vtkh5
{

using namespace hdf5::h5;

template<typename Writer>
class H5TypeWriter : public PHARE::diagnostic::TypeWriter
{
    using FloatType                            = float;
    using ModelView                            = Writer::ModelView;
    using physical_quantity_type               = ModelView::Field::physical_quantity_type;
    std::string static inline const base       = "/VTKHDF/";
    std::string static inline const level_base = base + "Level";
    std::string static inline const step_level = base + "Steps/Level";
    auto static inline const level_data_path
        = [](auto const ilvl) { return level_base + std::to_string(ilvl) + "/PointData/data"; };

public:
    static constexpr auto dimension = Writer::dimension;

    H5TypeWriter(Writer& h5Writer)
        : h5Writer_{h5Writer}
    {
    }

protected:
    // data is duplicated in lower dimensions to fit a 3d view
    //  so 2d data is * 2
    //  and 1d data is * 4
    constexpr static auto X_TIMES = std::array{4, 2, 1}[dimension - 1];

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
    class VTKFileWriter;


    auto& getOrCreateH5File(DiagnosticProperties const& diagnostic)
    {
        if (!fileExistsFor(diagnostic))
            fileData_.emplace(diagnostic.quantity, this->h5Writer_.makeFile(diagnostic));
        return *fileData_.at(diagnostic.quantity);
    }

    bool fileExistsFor(DiagnosticProperties const& diagnostic) const
    {
        return fileData_.count(diagnostic.quantity);
    }


    Writer& h5Writer_;
    std::unordered_map<std::string, std::unique_ptr<HighFiveFile>> fileData_;
};




template<typename Writer>
class H5TypeWriter<Writer>::VTKFileWriter
{
    auto static constexpr primal_qty         = physical_quantity_type::rho;
    std::size_t static constexpr boxValsIn3D = 6; // lo0, up0, lo1, up1, lo2, up2


public:
    VTKFileWriter(DiagnosticProperties const& prop, H5TypeWriter<Writer>* const typewriter);

    void writeField(auto const& field, auto const& layout);

    template<std::size_t rank = 2>
    void writeTensorField(auto const& tf, auto const& layout);

    void initFileLevel(int const ilvl) const;

    void initFieldFileLevel(auto const& level) { initAnyFieldLevel(level); }

    template<std::size_t rank = 2>
    void initTensorFieldFileLevel(auto const& level)
    {
        initAnyFieldLevel<core::detail::tensor_field_dim_from_rank<rank>()>(level);
    }

private:
    auto level_spacing(int const lvl) const
    {
        auto const mesh_size = typewriter->h5Writer_.modelView().cellWidth();
        return core::for_N_make_array<dimension>(
            [&](auto i) { return static_cast<float>(mesh_size[i] / std::pow(2, lvl)); });
    }

    auto local_box(auto const& layout) const
    {
        return layout.AMRToLocal(layout.AMRBoxFor(primal_qty));
    }

    template<std::size_t N = 1>
    void initAnyFieldLevel(auto const& level)
    {
        h5file.create_resizable_2d_data_set<FloatType, N>(level_data_path(level.getLevelNumber()));
        resize<N>(level);
    }


    template<std::size_t N = 1>
    void resize_data(int const ilvl, auto const& boxes);

    void resize_boxes(auto const& level, auto const& boxes);

    template<std::size_t N = 1>
    void resize(auto const& level)
    {
        auto const ilvl = level.getLevelNumber();
        auto boxes      = modelView.localLevelBoxes(ilvl);
        resize_boxes(level, boxes);
        for (auto& box : boxes)
            box.upper += 1; // primal
        resize_data<N>(ilvl, boxes);
    }

    bool newFile;
    DiagnosticProperties const& diagnostic;
    H5TypeWriter<Writer>* const typewriter;
    HighFiveFile& h5file;
    std::size_t data_offset = 0;
    ModelView& modelView    = typewriter->h5Writer_.modelView();
};



template<typename Writer>
H5TypeWriter<Writer>::VTKFileWriter::VTKFileWriter(DiagnosticProperties const& prop,
                                                   H5TypeWriter<Writer>* const tw)
    : newFile{!tw->fileExistsFor(prop)}
    , diagnostic{prop}
    , typewriter{tw}
    , h5file{tw->getOrCreateH5File(prop)}
{
    PHARE_LOG_SCOPE(1, "VTKFileWriter::VTKFileWriter");

    // set global per file attributes and datasets
    {
        h5file.create_resizable_1d_data_set<FloatType>(base + "/Steps/Values");

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

    if (newFile)
    {
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
void H5TypeWriter<Writer>::VTKFileWriter::writeField(auto const& field, auto const& layout)
{
    PHARE_LOG_SCOPE(1, "VTKFileWriter::writeField");

    auto const& frimal = core::convert_to_fortran_primal(modelView.tmpField(), field, layout);
    auto const size    = local_box(layout).size();
    auto ds            = h5file.getDataSet(level_data_path(layout.levelNumber()));

    for (std::uint16_t i = 0; i < X_TIMES; ++i)
    {
        ds.select({data_offset, 0}, {size, 1}).write_raw(frimal.data());
        data_offset += size;
    }
}

template<typename Writer>
template<std::size_t rank>
void H5TypeWriter<Writer>::VTKFileWriter::writeTensorField(auto const& tf, auto const& layout)
{
    PHARE_LOG_SCOPE(1, "VTKFileWriter::writeTensorField");

    auto static constexpr N = core::detail::tensor_field_dim_from_rank<rank>();
    auto const& frimal
        = core::convert_to_fortran_primal(modelView.template tmpTensorField<rank>(), tf, layout);
    auto const size = local_box(layout).size();
    auto ds         = h5file.getDataSet(level_data_path(layout.levelNumber()));
    for (std::uint16_t i = 0; i < X_TIMES; ++i)
    {
        for (std::uint32_t c = 0; c < N; ++c)
            ds.select({data_offset, c}, {size, 1}).write_raw(frimal[c].data());
        data_offset += size;
    }
}



template<typename Writer>
void H5TypeWriter<Writer>::VTKFileWriter::initFileLevel(int const ilvl) const
{
    PHARE_LOG_SCOPE(1, "VTKFileWriter::initFileLevel");

    auto const lvl = std::to_string(ilvl);

    h5file.create_resizable_2d_data_set<int, boxValsIn3D>(level_base + lvl + "/AMRBox");
    h5file.create_resizable_1d_data_set<int>(step_level + lvl + "/AMRBoxOffset");
    h5file.create_resizable_1d_data_set<int>(step_level + lvl + "/NumberOfAMRBox");
    h5file.create_resizable_1d_data_set<int>(step_level + lvl + "/PointDataOffset/data");

    auto level_group = h5file.file().getGroup(level_base + lvl);
    if (!level_group.hasAttribute("Spacing"))
    {
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
void H5TypeWriter<Writer>::VTKFileWriter::resize_data(int const ilvl, auto const& boxes)
{
    PHARE_LOG_SCOPE(1, "VTKFileWriter::resize_data");

    auto const lvl       = std::to_string(ilvl);
    auto const data_path = level_data_path(ilvl);
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

    point_data_ds.resize({new_size, N});
    for (int i = 0; i < core::mpi::rank(); ++i)
        data_offset += rank_data_size[i];
}


template<typename Writer>
void H5TypeWriter<Writer>::VTKFileWriter::resize_boxes(auto const& level, auto const& boxes)
{
    PHARE_LOG_SCOPE(1, "VTKFileWriter::resize_boxes");

    auto const ilvl          = level.getLevelNumber();
    auto const lvl           = std::to_string(ilvl);
    auto const rank_box_size = core::mpi::collect(boxes.size());
    auto const total_boxes   = core::sum(rank_box_size);

    {
        auto ds             = h5file.getDataSet(step_level + lvl + "/NumberOfAMRBox");
        auto const old_size = ds.getDimensions()[0];
        ds.resize({old_size + 1});
        ds.select({old_size}, {1}).write(total_boxes);
    }

    auto amrbox_ds  = h5file.getDataSet(level_base + lvl + "/AMRBox");
    auto box_offset = amrbox_ds.getDimensions()[0];

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


} // namespace PHARE::diagnostic::vtkh5

#endif // PHARE_DIAGNOSTIC_DETAIL_VTK_H5_TYPE_WRITER_HPP
