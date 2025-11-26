#ifndef PHARE_DIAGNOSTIC_DETAIL_VTK_H5_WRITER_HPP
#define PHARE_DIAGNOSTIC_DETAIL_VTK_H5_WRITER_HPP


#include "core/logger.hpp"
#include "core/utilities/mpi_utils.hpp"

#include "initializer/data_provider.hpp"

#include "diagnostic/diagnostic_props.hpp"
#include "diagnostic/detail/vtkh5_type_writer.hpp"

#include "diagnostic/detail/vtk_types/fluid.hpp"
#include "diagnostic/detail/vtk_types/electromag.hpp"

#include "hdf5/detail/h5/h5_file.hpp"



namespace PHARE::diagnostic::vtkh5
{
using namespace hdf5::h5;



template<typename _ModelView>
class H5Writer
{
    constexpr std::size_t static MAX_LEVEL = 10;

    struct NullTypeWriter : public H5TypeWriter<H5Writer<_ModelView>>
    {
        NullTypeWriter(auto& h5Writer)
            : H5TypeWriter<H5Writer<_ModelView>>{h5Writer}
        {
        }

        void setup(DiagnosticProperties& prop) {}
        void write(DiagnosticProperties& prop)
        {
            if (core::mpi::rank() == 0)
            {
                PHARE_LOG_LINE_SS( //
                    "No diagnostic writer found for " + prop.type + ":" + prop.quantity);
            }
        }
        void compute(DiagnosticProperties&) {}
    };

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


    auto makeFile(DiagnosticProperties const& diagnostic)
    {
        return std::make_unique<HighFiveFile>(filePath_ + "/" + fileString(diagnostic.quantity),
                                              file_flags[diagnostic.type + diagnostic.quantity]);
    }


    auto& modelView() { return modelView_; }
    auto timestamp() const { return timestamp_; }

    std::size_t minLevel = 0, maxLevel = MAX_LEVEL; // TODO hard-coded to be parametrized somehow
    HiFile::AccessMode flags;


private:
    double timestamp_ = 0;
    std::string filePath_;
    ModelView modelView_;
    Attributes fileAttributes_;

    std::unordered_map<std::string, HiFile::AccessMode> file_flags;

    std::unordered_map<std::string, std::shared_ptr<H5TypeWriter<This>>> typeWriters_{
        {"info", make_writer<NullTypeWriter>()},
        {"meta", make_writer<NullTypeWriter>()},
        {"fluid", make_writer<FluidDiagnosticWriter<This>>()},
        {"electromag", make_writer<ElectromagDiagnosticWriter<This>>()},
        {"particle", make_writer<NullTypeWriter>()} //
    };

    template<typename Writer>
    std::shared_ptr<H5TypeWriter<This>> make_writer()
    {
        return std::make_shared<Writer>(*this);
    }


    static std::string fileString(std::string fileStr)
    {
        if (fileStr[0] == '/')
            fileStr = fileStr.substr(1);
        std::replace(fileStr.begin(), fileStr.end(), '/', '_');
        return fileStr + ".vtkhdf";
    }


    //  State of this class is controlled via "dump()"
    //  block public access to internal state
    friend class H5TypeWriter<This>;
    friend class FluidDiagnosticWriter<This>;
    friend class ElectromagDiagnosticWriter<This>;
};



template<typename ModelView>
void H5Writer<ModelView>::dump(std::vector<DiagnosticProperties*> const& diagnostics,
                               double timestamp)
{
    timestamp_ = timestamp;

    for (auto* diagnostic : diagnostics)
        if (!file_flags.count(diagnostic->type + diagnostic->quantity))
            file_flags[diagnostic->type + diagnostic->quantity] = this->flags;

    for (auto* diagnostic : diagnostics) // all collective calls first!
        typeWriters_.at(diagnostic->type)->setup(*diagnostic);

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
