#ifndef PHARE_DETAIL_DIAGNOSTIC_HIGHFIVE_HPP
#define PHARE_DETAIL_DIAGNOSTIC_HIGHFIVE_HPP


#include "core/data/vecfield/vecfield_component.hpp"
#include "core/utilities/mpi_utils.hpp"
#include "core/utilities/types.hpp"
#include "core/utilities/meta/meta_utilities.hpp"

#include "hdf5/detail/h5/h5_file.hpp"

#include "diagnostic/detail/h5typewriter.hpp"
#include "diagnostic/diagnostic_manager.hpp"
#include "diagnostic/diagnostic_props.hpp"


#if !defined(PHARE_DIAG_DOUBLES)
#error // PHARE_DIAG_DOUBLES not defined
#endif


namespace PHARE::diagnostic::h5
{
using namespace hdf5::h5;

template<typename Writer>
class ElectromagDiagnosticWriter;
template<typename Writer>
class FluidDiagnosticWriter;
template<typename Writer>
class ParticlesDiagnosticWriter;
template<typename Writer>
class MetaDiagnosticWriter;
template<typename Writer>
class InfoDiagnosticWriter;



template<typename ModelView>
class H5Writer
{
    using FloatType = std::conditional_t<PHARE_DIAG_DOUBLES, double, float>;

    static constexpr std::size_t timestamp_precision = 10;

public:
    using This       = H5Writer<ModelView>;
    using GridLayout = typename ModelView::GridLayout;
    using Attributes = typename ModelView::PatchProperties;

    static constexpr auto dimension   = GridLayout::dimension;
    static constexpr auto interpOrder = GridLayout::interp_order;
    static constexpr auto READ_WRITE  = HiFile::ReadWrite | HiFile::Create;

    // flush_never: disables manual file closing, but still occurrs via RAII
    static constexpr std::size_t flush_never = 0;

    template<typename Hierarchy, typename Model>
    H5Writer(Hierarchy& hier, Model& model, std::string const hifivePath,
             unsigned _flags /* = HiFile::ReadWrite | HiFile::Create | HiFile::Truncate */)
        : flags{_flags}
        , filePath_{hifivePath}
        , modelView_{hier, model}
    {
    }

    ~H5Writer() {}

    template<typename Hierarchy, typename Model>
    static auto make_unique(Hierarchy& hier, Model& model, initializer::PHAREDict const& dict)
    {
        std::string filePath = dict["filePath"].template to<std::string>();
        unsigned flags       = READ_WRITE;
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
        return fileStr + ".h5";
    }

    auto makeFile(std::string const filename, unsigned file_flag)
    {
        return std::make_unique<HighFiveFile>(filePath_ + "/" + filename, file_flag);
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

    auto& modelView() { return modelView_; }

    std::size_t minLevel = 0, maxLevel = 10; // TODO hard-coded to be parametrized somehow
    unsigned flags;


private:
    double timestamp_ = 0;
    std::string filePath_;
    std::string patchPath_; // is passed around as "virtual write()" has no parameters
    ModelView modelView_;
    Attributes fileAttributes_;

    std::unordered_map<std::string, unsigned> file_flags;

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


    void initializeDatasets_(std::vector<DiagnosticProperties*> const& diagnotics);
    void writeDatasets_(std::vector<DiagnosticProperties*> const& diagnotics);

    H5Writer(H5Writer const&)            = delete;
    H5Writer(H5Writer&&)                 = delete;
    H5Writer& operator&(H5Writer const&) = delete;
    H5Writer& operator&(H5Writer&&)      = delete;


    //  State of this class is controlled via "dump()"
    //  block public access to internal state
    friend class FluidDiagnosticWriter<This>;
    friend class ElectromagDiagnosticWriter<This>;
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



template<typename ModelView>
void H5Writer<ModelView>::dump(std::vector<DiagnosticProperties*> const& diagnostics,
                               double timestamp)
{
    timestamp_                     = timestamp;
    fileAttributes_["dimension"]   = dimension;
    fileAttributes_["interpOrder"] = interpOrder;
    fileAttributes_["layoutType"]  = modelView_.getLayoutTypeString();
    fileAttributes_["domain_box"]  = modelView_.domainBox();
    fileAttributes_["cell_width"]  = modelView_.cellWidth();
    fileAttributes_["origin"]      = modelView_.origin();

    fileAttributes_["boundary_conditions"] = modelView_.boundaryConditions();

    for (auto* diagnostic : diagnostics)
        if (!file_flags.count(diagnostic->type + diagnostic->quantity))
            file_flags[diagnostic->type + diagnostic->quantity] = this->flags;

    initializeDatasets_(diagnostics);
    writeDatasets_(diagnostics);

    for (auto* diagnostic : diagnostics)
    {
        typeWriters_.at(diagnostic->type)->finalize(*diagnostic);
        // don't truncate past first dump
        file_flags[diagnostic->type + diagnostic->quantity] = READ_WRITE;
    }
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




template<typename ModelView>
void H5Writer<ModelView>::initializeDatasets_(std::vector<DiagnosticProperties*> const& diagnostics)
{
    std::size_t maxLocalLevel = 0;
    std::unordered_map<std::size_t, std::vector<std::string>> lvlPatchIDs;
    Attributes patchAttributes; // stores dataset info/size for synced MPI creation

    for (auto* diag : diagnostics)
        typeWriters_.at(diag->type)->createFiles(*diag);

    auto collectPatchAttributes = [&](GridLayout&, std::string patchID, std::size_t iLevel) {
        if (!lvlPatchIDs.count(iLevel))
            lvlPatchIDs.emplace(iLevel, std::vector<std::string>());

        lvlPatchIDs.at(iLevel).emplace_back(patchID);

        for (auto* diag : diagnostics)
        {
            typeWriters_.at(diag->type)->getDataSetInfo(*diag, iLevel, patchID, patchAttributes);
        }
        maxLocalLevel = iLevel;
    };

    modelView_.visitHierarchy(collectPatchAttributes, minLevel, maxLevel);

    // sets empty vectors in case current process lacks patch on a level
    std::size_t maxMPILevel = core::mpi::max(maxLocalLevel);
    for (std::size_t lvl = minLevel; lvl <= maxMPILevel; lvl++)
        if (!lvlPatchIDs.count(lvl))
            lvlPatchIDs.emplace(lvl, std::vector<std::string>());

    for (auto* diagnostic : diagnostics)
    {
        typeWriters_.at(diagnostic->type)
            ->initDataSets(*diagnostic, lvlPatchIDs, patchAttributes, maxMPILevel);
    }
}



template<typename ModelView>
void H5Writer<ModelView>::writeDatasets_(std::vector<DiagnosticProperties*> const& diagnostics)
{
    std::unordered_map<std::size_t, std::vector<std::pair<std::string, Attributes>>>
        patchAttributes;

    std::size_t maxLocalLevel = 0;
    auto writePatch = [&](GridLayout& gridLayout, std::string patchID, std::size_t iLevel) {
        if (!patchAttributes.count(iLevel))
            patchAttributes.emplace(iLevel, std::vector<std::pair<std::string, Attributes>>{});
        patchPath_ = getPatchPathAddTimestamp(iLevel, patchID);
        patchAttributes[iLevel].emplace_back(patchID,
                                             modelView_.getPatchProperties(patchID, gridLayout));
        for (auto* diagnostic : diagnostics)
            typeWriters_.at(diagnostic->type)->write(*diagnostic);
        maxLocalLevel = iLevel;
    };

    modelView_.visitHierarchy(writePatch, minLevel, maxLevel);

    std::size_t maxMPILevel = core::mpi::max(maxLocalLevel);
    // sets empty vectors in case current process lacks patch on a level
    for (std::size_t lvl = minLevel; lvl <= maxMPILevel; lvl++)
        if (!patchAttributes.count(lvl))
            patchAttributes.emplace(lvl, std::vector<std::pair<std::string, Attributes>>{});

    for (auto* diagnostic : diagnostics)
        typeWriters_.at(diagnostic->type)
            ->writeAttributes(*diagnostic, fileAttributes_, patchAttributes, maxMPILevel);
}



} /* namespace PHARE::diagnostic::h5 */

#endif /* PHARE_DETAIL_DIAGNOSTIC_HIGHFIVE_H */
