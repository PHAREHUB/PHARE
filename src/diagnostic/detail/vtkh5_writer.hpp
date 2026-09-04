#ifndef PHARE_DIAGNOSTIC_DETAIL_VTK_H5_WRITER_HPP
#define PHARE_DIAGNOSTIC_DETAIL_VTK_H5_WRITER_HPP

#include "core/logger.hpp"
#include "mpi/mpi_utils.hpp"

#include "initializer/data_provider.hpp"

#include "diagnostic/diagnostic_props.hpp"
#include "diagnostic/detail/vtkh5_type_writer.hpp"

#include "diagnostic/detail/vtk_types/fluid.hpp"
#include "diagnostic/detail/vtk_types/electromag.hpp"

#include "hdf5/detail/h5/h5_file.hpp"

#include <memory>
#include <string>
#include <algorithm>
#include <unordered_map>

namespace PHARE::diagnostic::vtkh5
{
using namespace hdf5::h5;

template<typename ModelMapper_t>
class H5Writer
{
    struct NullTypeWriter : public H5TypeWriter<H5Writer<ModelMapper_t>>
    {
        NullTypeWriter(auto& h5Writer)
            : H5TypeWriter<H5Writer<ModelMapper_t>>{h5Writer}
        {
        }

        void setup(DiagnosticProperties& prop) {}
        void write(DiagnosticProperties& prop)
        {
            if (mpi::rank() == 0)
            {
                PHARE_LOG_LINE_SS( //
                    "No diagnostic writer found for " + prop.type + ":" + prop.quantity);
            }
        }
        void compute(DiagnosticProperties&) {}
    };

public:
    // using ModelView  = _ModelView;
    using This       = H5Writer<ModelMapper_t>;
    using Attributes = PatchProperties;

    static constexpr auto dimension  = ModelMapper_t::dimension;
    static constexpr auto READ_WRITE = HiFile::AccessMode::OpenOrCreate;

    // flush_never: disables manual file closing, but still occurrs via RAII
    static constexpr std::size_t flush_never = 0;

    struct Config
    {
        std::string filePath;
        HiFile::AccessMode flags;

        static Config FROM(initializer::PHAREDict const& dict)
        {
            auto flags = READ_WRITE;
            if (dict.contains("mode") and dict["mode"].template to<std::string>() == "overwrite")
                flags |= HiFile::Truncate;
            return {dict["filePath"].template to<std::string>(), flags};
        }
    };
    using Config_t = Config;

    H5Writer(auto& hier, Config const& config, auto&... models)
        : mapper_{hier, models...}
        , config_{config}
    {
        if constexpr (ModelMapper_t::has_hybrid_model)
        {
            typeWriters_.emplace("info", make_writer<NullTypeWriter>());
            typeWriters_.emplace("meta", make_writer<NullTypeWriter>());
            typeWriters_.emplace("fluid", make_writer<FluidDiagnosticWriter<This>>());
            typeWriters_.emplace("electromag", make_writer<ElectromagDiagnosticWriter<This>>());
            typeWriters_.emplace("particle", make_writer<NullTypeWriter>());
        }
        if constexpr (ModelMapper_t::has_mhd_model)
        {
            typeWriters_.emplace("meta", make_writer<NullTypeWriter>());
            typeWriters_.emplace("mhd", make_writer<FluidDiagnosticWriter<This>>());
            typeWriters_.emplace("electromag", make_writer<ElectromagDiagnosticWriter<This>>());
        }
        if constexpr (!ModelMapper_t::has_hybrid_model and !ModelMapper_t::has_mhd_model)
        {
            // MacOS clang unhappy with static_assert(false), requires a dependency on Model
            static_assert(core::dependent_false_v<ModelMapper_t>,
                          "Unsupported model type in H5Writer");
        }
    }

    ~H5Writer() {}
    H5Writer(H5Writer const&)            = delete;
    H5Writer(H5Writer&&)                 = delete;
    H5Writer& operator=(H5Writer const&) = delete;
    H5Writer& operator=(H5Writer&&)      = delete;


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
        return std::make_unique<HighFiveFile>(config_.filePath + "/"
                                                  + fileString(diagnostic.quantity),
                                              file_flags[diagnostic.type + diagnostic.quantity]);
    }

    template<typename Dict>
    static void writeGlobalAttributeDict(HighFiveFile& h5, Dict const& dict, std::string path)
    {
        dict.visit(
            [&](std::string const& key, auto const& val) { h5.write_attribute(path, key, val); });
    }


    auto& mapper() { return mapper_; }
    auto timestamp() const { return timestamp_; }

private:
    ModelMapper_t mapper_;
    Config config_;

public:
    std::size_t minLevel = 0, maxLevel = mapper_.maxLevel();


private:
    double timestamp_ = 0;
    Attributes fileAttributes_;

    std::unordered_map<std::string, HiFile::AccessMode> file_flags;

    std::unordered_map<std::string, std::shared_ptr<H5TypeWriter<This>>> typeWriters_;

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



template<typename ModelMapper_t>
void H5Writer<ModelMapper_t>::dump(std::vector<DiagnosticProperties*> const& diagnostics,
                                   double timestamp)
{
    timestamp_                   = timestamp;
    fileAttributes_["dimension"] = dimension;
    if constexpr (ModelMapper_t::has_hybrid_model)
        fileAttributes_["interpOrder"]
            = ModelMapper_t::HasModels_t::HybridModel_t::GridLayoutT::options.interp_order;
    fileAttributes_["domain_box"]          = mapper_.hierarchy.domainBox();
    fileAttributes_["boundary_conditions"] = mapper_.hierarchy.boundaryConditions();

    HierarchyData<dimension>::reset(*this);

    for (auto* diagnostic : diagnostics)
        if (!file_flags.count(diagnostic->type + diagnostic->quantity))
            file_flags[diagnostic->type + diagnostic->quantity] = config_.flags;

    for (auto* diagnostic : diagnostics) // all collective calls first!
    {
        auto& type_writer = *typeWriters_.at(diagnostic->type);
        type_writer.setup(*diagnostic);
        type_writer.writeFileAttributes(*diagnostic, fileAttributes_);
    }

    for (auto* diagnostic : diagnostics)
    {
        auto& typeWriter = *typeWriters_.at(diagnostic->type);
        typeWriter.compute(*diagnostic); // compute to temporaries then write immediately!
        typeWriter.write(*diagnostic);
    }

    for (auto* diagnostic : diagnostics)
    {
        typeWriters_.at(diagnostic->type)->finalize(*diagnostic);
        // don't truncate past first dump
        file_flags[diagnostic->type + diagnostic->quantity] = READ_WRITE;
    }
}

template<typename ModelMapper_t>
void H5Writer<ModelMapper_t>::dump_level(std::size_t level,
                                         std::vector<DiagnosticProperties*> const& diagnostics,
                                         double timestamp)
{
    throw std::runtime_error("!!DOES NOT WORK!!");
}



} // namespace PHARE::diagnostic::vtkh5

#endif /* PHARE_DIAGNOSTIC_DETAIL_VTK_H5_WRITER_HPP */
