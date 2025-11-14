#ifndef PHARE_DIAGNOSTIC_DETAIL_VTK_H5_WRITER_HPP
#define PHARE_DIAGNOSTIC_DETAIL_VTK_H5_WRITER_HPP


#include "core/utilities/types.hpp"
#include "core/utilities/mpi_utils.hpp"
#include "core/data/tensorfield/tensorfield.hpp"
#include "core/utilities/meta/meta_utilities.hpp"
#include "core/data/vecfield/vecfield_component.hpp"

#include "initializer/data_provider.hpp"

#include "hdf5/detail/h5/h5_file.hpp"

#include "diagnostic/diagnostic_props.hpp"
#include "diagnostic/detail/vtkh5_type_writer.hpp"
#include "diagnostic/detail/vtk_types/electromag.hpp"
#include "diagnostic/detail/vtk_types/fluid.hpp"
#include "diagnostic/detail/vtk_types/info.hpp"
#include "diagnostic/detail/vtk_types/meta.hpp"
#include "diagnostic/detail/vtk_types/particle.hpp"


#if !defined(PHARE_DIAG_DOUBLES)
#error // PHARE_DIAG_DOUBLES not defined
#endif


namespace PHARE::diagnostic::vtkh5
{
using namespace hdf5::h5;


template<typename _ModelView>
class H5Writer
{
    using FloatType = std::conditional_t<PHARE_DIAG_DOUBLES, double, float>;

public:
    using ModelView  = _ModelView;
    using This       = H5Writer<ModelView>;
    using GridLayout = ModelView::GridLayout;
    using Attributes = ModelView::PatchProperties;

    static constexpr auto dimension   = GridLayout::dimension;
    static constexpr auto interpOrder = GridLayout::interp_order;
    static constexpr auto READ_WRITE  = HiFile::AccessMode::OpenOrCreate;

    // flush_never: disables manual file closing, but still occurrs via RAII
    static constexpr std::size_t flush_never = 0;

    template<typename Hierarchy, typename Model>
    H5Writer(Hierarchy& hier, Model& model, std::string const hifivePath, HiFile::AccessMode _flags)
        : flags{_flags}
        , filePath_{hifivePath}
        , modelView_{hier, model}
    {
    }
    ~H5Writer() {}
    H5Writer(H5Writer const&)            = delete;
    H5Writer(H5Writer&&)                 = delete;
    H5Writer& operator=(H5Writer const&) = delete;
    H5Writer& operator=(H5Writer&&)      = delete;


    template<typename Hierarchy, typename Model>
    static auto make_unique(Hierarchy& hier, Model& model, initializer::PHAREDict const& dict)
    {
        std::string filePath     = dict["filePath"].template to<std::string>();
        HiFile::AccessMode flags = READ_WRITE;
        if (dict.contains("mode") and dict["mode"].template to<std::string>() == "overwrite")
            flags |= HiFile::Truncate;
        return std::make_unique<This>(hier, model, filePath, flags);
    }


    void dump(std::vector<DiagnosticProperties*> const&, double current_timestamp);
    void dump_level(std::size_t level, std::vector<DiagnosticProperties*> const& diagnostics,
                    double timestamp);

    template<typename String>
    auto getDiagnosticWriterForType(String& type)
    {
        return typeWriters_.at(type);
    }


    static std::string fileString(std::string fileStr)
    {
        if (fileStr[0] == '/')
            fileStr = fileStr.substr(1);
        std::replace(fileStr.begin(), fileStr.end(), '/', '_');
        return fileStr + ".vtkhdf";
    }


    auto makeFile(DiagnosticProperties const& diagnostic)
    {
        return std::make_unique<HighFiveFile>(filePath_ + "/" + fileString(diagnostic.quantity),
                                              file_flags[diagnostic.type + diagnostic.quantity]);
    }


    static std::string getFullPatchPath(std::string timestamp, int iLevel, std::string globalCoords)
    {
        return "/t/" + timestamp + "/pl" + std::to_string(iLevel) + "/p" + globalCoords;
    }

    template<typename Type, typename Size>
    static void createDataSet(HighFiveFile& h5, std::string const& path, Size const& size)
    {
        if constexpr (std::is_same_v<Type, double>) // force doubles for floats for storage
            h5.create_data_set_per_mpi<FloatType>(path, size);
        else
            h5.create_data_set_per_mpi<Type>(path, size);
    }



    template<typename TensorField>
    static void writeTensorFieldAsDataset(HighFiveFile& h5, std::string path, TensorField& tField)
    {
        for (auto& [id, type] : core::Components::componentMap<TensorField::rank>())
            h5.write_data_set_flat<dimension>(path + "_" + id, tField.getComponent(type).data());
    }

    auto& modelView() { return modelView_; }
    auto timestamp() const { return timestamp_; }

    std::size_t minLevel = 0, maxLevel = 10; // TODO hard-coded to be parametrized somehow
    HiFile::AccessMode flags;


private:
    double timestamp_ = 0;
    std::string filePath_;
    ModelView modelView_;
    Attributes fileAttributes_;

    std::unordered_map<std::string, HiFile::AccessMode> file_flags;

    std::unordered_map<std::string, std::shared_ptr<H5TypeWriter<This>>> typeWriters_{
        {"info", make_writer<InfoDiagnosticWriter<This>>()},
        {"meta", make_writer<MetaDiagnosticWriter<This>>()},
        {"fluid", make_writer<FluidDiagnosticWriter<This>>()},
        {"electromag", make_writer<ElectromagDiagnosticWriter<This>>()},
        {"particle", make_writer<ParticlesDiagnosticWriter<This>>()} //
    };

    template<typename Writer>
    std::shared_ptr<H5TypeWriter<This>> make_writer()
    {
        return std::make_shared<Writer>(*this);
    }




    //  State of this class is controlled via "dump()"
    //  block public access to internal state
    friend class FluidDiagnosticWriter<This>;
    friend class ElectromagDiagnosticWriter<This>;
    friend class ParticlesDiagnosticWriter<This>;
    friend class MetaDiagnosticWriter<This>;
    friend class InfoDiagnosticWriter<This>;
    friend class H5TypeWriter<This>;
};



template<typename ModelView>
void H5Writer<ModelView>::dump(std::vector<DiagnosticProperties*> const& diagnostics,
                               double timestamp)
{
    timestamp_ = timestamp;

    for (auto* diagnostic : diagnostics)
        if (!file_flags.count(diagnostic->type + diagnostic->quantity))
            file_flags[diagnostic->type + diagnostic->quantity] = this->flags;

    for (auto* diagnostic : diagnostics)
        typeWriters_.at(diagnostic->type)->write(*diagnostic);

    for (auto* diagnostic : diagnostics) // don't truncate past first dump
        file_flags[diagnostic->type + diagnostic->quantity] = READ_WRITE;
}

template<typename ModelView>
void H5Writer<ModelView>::dump_level(std::size_t level,
                                     std::vector<DiagnosticProperties*> const& diagnostics,
                                     double timestamp)
{
    std::size_t _minLevel = this->minLevel;
    std::size_t _maxLevel = this->maxLevel;

    this->minLevel = level;
    this->maxLevel = level;

    this->dump(diagnostics, timestamp);

    this->minLevel = _minLevel;
    this->maxLevel = _maxLevel;
}



} // namespace PHARE::diagnostic::vtkh5

#endif /* PHARE_DIAGNOSTIC_DETAIL_VTK_H5_WRITER_HPP */
