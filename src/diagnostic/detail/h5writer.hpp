#ifndef PHARE_DETAIL_DIAGNOSTIC_HIGHFIVE_HPP
#define PHARE_DETAIL_DIAGNOSTIC_HIGHFIVE_HPP

#include "core/utilities/types.hpp"
#include "core/data/vecfield/vecfield_component.hpp"

#include "mpi/mpi_utils.hpp"

#include "initializer/data_provider.hpp"

#include "diagnostic/diagnostic_props.hpp"
#include "diagnostic/diagnostic_model_view.hpp"

#include "hdf5/detail/h5/h5_file.hpp"

#include "diagnostic/diagnostic_props.hpp"
#include "diagnostic/detail/h5typewriter.hpp"
#include "diagnostic/detail/types/info.hpp"
#include "diagnostic/detail/types/meta.hpp"
#include "diagnostic/detail/types/fluid.hpp"
#include "diagnostic/detail/types/particle.hpp"
#include "diagnostic/detail/types/electromag.hpp"
#include "diagnostic/detail/types/mhd.hpp"


#if !defined(PHARE_DIAG_DOUBLES)
#error // PHARE_DIAG_DOUBLES not defined
#endif


namespace PHARE::diagnostic::h5
{
using namespace hdf5::h5;



template<typename ModelMapper_>
class H5Writer
{
    using FloatType = std::conditional_t<PHARE_DIAG_DOUBLES, double, float>;

    static constexpr std::size_t timestamp_precision = 10;

public:
    using This             = H5Writer<ModelMapper_>;
    using ModelMapper_t    = ModelMapper_;
    using ModelViewVariant = typename ModelMapper_t::ModelViewVariant;
    using Attributes       = PatchProperties;


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
        // several models can share this writer (e.g. MHD + Hybrid) - register the type
        // writers applicable to whichever models are actually present; shared types
        // (meta/electromag) are simply not re-inserted if already registered.
        core::apply(mapper_.models, [&](auto& model) {
            using Model_t = std::decay_t<decltype(model)>;
            if constexpr (solver::is_hybrid_model_v<Model_t>)
            {
                typeWriters_.emplace("info", make_writer<InfoDiagnosticWriter<This>>());
                typeWriters_.emplace("meta", make_writer<MetaDiagnosticWriter<This>>());
                typeWriters_.emplace("fluid", make_writer<FluidDiagnosticWriter<This>>());
                typeWriters_.emplace("electromag", make_writer<ElectromagDiagnosticWriter<This>>());
                typeWriters_.emplace("particle", make_writer<ParticlesDiagnosticWriter<This>>());
            }
            else if constexpr (solver::is_mhd_model_v<Model_t>)
            {
                typeWriters_.emplace("meta", make_writer<MetaDiagnosticWriter<This>>());
                typeWriters_.emplace("mhd", make_writer<MHDDiagnosticWriter<This>>());
                typeWriters_.emplace("electromag", make_writer<ElectromagDiagnosticWriter<This>>());
            }
            else
                static_assert(core::dependent_false_v<Model_t>,
                              "Unsupported model type in H5Writer");
        });
    }

    ~H5Writer() {}

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
        return fileStr + ".h5";
    }

    auto makeFile(std::string const filename, HiFile::AccessMode const file_flag)
    {
        return std::make_unique<HighFiveFile>(config_.filePath + "/" + filename, file_flag);
    }

    auto makeFile(DiagnosticProperties const& diagnostic)
    {
        return makeFile(fileString(diagnostic.quantity),
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


    // global function when all path+key are the same
    template<typename Data>
    static void writeAttribute(HighFiveFile& h5, std::string const& path, std::string const& key,
                               Data const& value)
    {
        h5.write_attribute(path, key, value);
    }


    template<typename Dict>
    static void writeAttributeDict(HighFiveFile& h5, Dict dict, std::string path)
    {
        dict.visit([&](std::string const& key, auto const& val) {
            auto constexpr static unsupported
                = std::is_same_v<std::decay_t<decltype(val)>, std::vector<std::string>>;

            // the dict might have types that are not supported, but not actually contain
            //  any of the types at runtime
            if constexpr (!unsupported)
                h5.write_attributes_per_mpi(path, key, val);

            // runtime detection of unsupported types
            if (unsupported)
                throw std::runtime_error("Unsupported operation: Cannot write attribute: path("
                                         + path + "), key(" + key + ")");
        });
    }

    template<typename Dict>
    static void writeGlobalAttributeDict(HighFiveFile& h5, Dict dict, std::string path)
    {
        dict.visit(
            [&](std::string const& key, auto const& val) { h5.write_attribute(path, key, val); });
    }



    template<typename TensorField>
    static void writeTensorFieldAsDataset(HighFiveFile& h5, std::string path, TensorField& tField)
    {
        for (auto& [id, type] : core::Components::componentMap<TensorField::rank>())
            h5.write_data_set_flat<dimension>(path + "_" + id, tField.getComponent(type).data());
    }

    auto& mapper() { return mapper_; }
    auto timestamp() const { return timestamp_; }

    // the model view for whichever patch is currently being written - set by writeDatasets_
    // right before write() is called, since write() itself has no patch/level parameter
    auto& currentModelView() { return *currentModelView_; }

private:
    ModelMapper_t mapper_;
    Config config_;
    ModelViewVariant* currentModelView_ = nullptr;

public:
    std::size_t minLevel = 0, maxLevel = mapper_.maxLevel();


private:
    double timestamp_ = 0;
    std::string patchPath_; // is passed around as "virtual write()" has no parameters
    Attributes fileAttributes_;

    std::unordered_map<std::string, HiFile::AccessMode> file_flags;

    std::unordered_map<std::string, std::shared_ptr<H5TypeWriter<This>>> typeWriters_;

    template<typename Writer>
    std::shared_ptr<H5TypeWriter<This>> make_writer()
    {
        return std::make_shared<Writer>(*this);
    }


    void initializeDatasets_(std::vector<DiagnosticProperties*> const& diagnotics);
    void writeDatasets_(std::vector<DiagnosticProperties*> const& diagnotics);

    H5Writer(H5Writer const&)            = delete;
    H5Writer(H5Writer&&)                 = delete;
    H5Writer& operator=(H5Writer const&) = delete;
    H5Writer& operator=(H5Writer&&)      = delete;


    //  State of this class is controlled via "dump()"
    //  block public access to internal state
    friend class FluidDiagnosticWriter<This>;
    friend class ElectromagDiagnosticWriter<This>;
    friend class MHDDiagnosticWriter<This>;
    friend class ParticlesDiagnosticWriter<This>;
    friend class MetaDiagnosticWriter<This>;
    friend class InfoDiagnosticWriter<This>;
    friend class H5TypeWriter<This>;

    // used by friends start
    std::string getPatchPathAddTimestamp(int iLevel, std::string globalCoords)
    {
        return getFullPatchPath(core::to_string_with_precision(timestamp_, timestamp_precision),
                                iLevel, globalCoords);
    }


    auto& patchPath() const { return patchPath_; }
    // used by friends end
};



template<typename ModelMapper_t>
void H5Writer<ModelMapper_t>::dump(std::vector<DiagnosticProperties*> const& diagnostics,
                                   double timestamp)
{
    timestamp_                   = timestamp;
    fileAttributes_["dimension"] = dimension;
    // fileAttributes_["layoutType"] = modelView_.getLayoutTypeString();
    // fileAttributes_["origin"]     = modelView_.origin();
    fileAttributes_["domain_box"]          = mapper_.hierarchy.domainBox();
    fileAttributes_["cell_width"]          = mapper_.cellWidth();
    fileAttributes_["boundary_conditions"] = mapper_.hierarchy.boundaryConditions();
    if constexpr (ModelMapper_t::has_hybrid_model)
        fileAttributes_["interpOrder"]
            = ModelMapper_t::HasModels_t::HybridModel_t::GridLayoutT::options.interp_order;

    for (auto* diagnostic : diagnostics)
        if (!file_flags.count(diagnostic->type + diagnostic->quantity))
            file_flags[diagnostic->type + diagnostic->quantity] = config_.flags;

    initializeDatasets_(diagnostics);
    writeDatasets_(diagnostics);

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
    std::size_t _minLevel = this->minLevel;
    std::size_t _maxLevel = this->maxLevel;

    this->minLevel = level;
    this->maxLevel = level;

    this->dump(diagnostics, timestamp);

    this->minLevel = _minLevel;
    this->maxLevel = _maxLevel;
}




template<typename ModelMapper_t>
void H5Writer<ModelMapper_t>::initializeDatasets_(
    std::vector<DiagnosticProperties*> const& diagnostics)
{
    std::size_t maxLocalLevel = 0;
    std::unordered_map<std::size_t, std::vector<std::string>> lvlPatchIDs;
    Attributes patchAttributes; // stores dataset info/size for synced MPI creation

    for (auto* diag : diagnostics)
        typeWriters_.at(diag->type)->createFiles(*diag);

    auto collectPatchAttributes
        = [&](auto& /*layout*/, std::string const& patchID, std::size_t iLevel, auto& modelView) {
              if (!lvlPatchIDs.count(iLevel))
                  lvlPatchIDs.emplace(iLevel, std::vector<std::string>());

              lvlPatchIDs.at(iLevel).emplace_back(patchID);

              for (auto* diag : diagnostics)
                  typeWriters_.at(diag->type)
                      ->getDataSetInfo(*diag, iLevel, patchID, patchAttributes, modelView);

              maxLocalLevel = iLevel;
          };

    mapper_.visitHierarchy(collectPatchAttributes, minLevel, maxLevel);

    // sets empty vectors in case current process lacks patch on a level
    std::size_t maxMPILevel = mpi::max(maxLocalLevel);
    for (std::size_t lvl = minLevel; lvl <= maxMPILevel; lvl++)
        if (!lvlPatchIDs.count(lvl))
            lvlPatchIDs.emplace(lvl, std::vector<std::string>());

    for (auto* diagnostic : diagnostics)
    {
        typeWriters_.at(diagnostic->type)
            ->initDataSets(*diagnostic, lvlPatchIDs, patchAttributes, maxMPILevel);
    }
}



template<typename ModelMapper_t>
void H5Writer<ModelMapper_t>::writeDatasets_(std::vector<DiagnosticProperties*> const& diagnostics)
{
    std::unordered_map<std::size_t, std::vector<std::pair<std::string, Attributes>>>
        patchAttributes;

    std::size_t maxLocalLevel = 0;

    auto collectPatch = [&](auto& gridLayout, std::string const& patchID, std::size_t iLevel,
                            auto& /*modelView*/) {
        if (!patchAttributes.count(iLevel))
            patchAttributes.emplace(iLevel, std::vector<std::pair<std::string, Attributes>>{});
        patchAttributes[iLevel].emplace_back(patchID, getPatchProperties(patchID, gridLayout));
        maxLocalLevel = iLevel;
    };
    mapper_.visitHierarchy(collectPatch, minLevel, maxLevel);

    for (auto* diagnostic : diagnostics)
    {
        auto& typeWriter = *typeWriters_.at(diagnostic->type);
        typeWriter.compute(*diagnostic); // compute to temporaries then write immediately!
        mapper_.visitHierarchy(
            [&](auto& /*gridLayout*/, std::string const& patchID, std::size_t iLevel,
                auto& modelView) {
                patchPath_        = getPatchPathAddTimestamp(iLevel, patchID);
                currentModelView_ = &modelView;
                typeWriter.write(*diagnostic);
            },
            minLevel, maxLevel);
    }

    std::size_t maxMPILevel = mpi::max(maxLocalLevel);
    // sets empty vectors in case current process lacks patch on a level
    for (std::size_t lvl = minLevel; lvl <= maxMPILevel; lvl++)
        if (!patchAttributes.count(lvl))
            patchAttributes.emplace(lvl, std::vector<std::pair<std::string, Attributes>>{});

    for (auto* diagnostic : diagnostics)
        typeWriters_.at(diagnostic->type)
            ->writeAttributes(*diagnostic, fileAttributes_, patchAttributes, maxMPILevel);
}



} /* namespace PHARE::diagnostic::h5 */

#endif /* PHARE_DETAIL_DIAGNOSTIC_HIGHFIVE_HPP */
