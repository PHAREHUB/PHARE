
#ifndef PHARE_DETAIL_DIAGNOSTIC_HIGHFIVE_H
#define PHARE_DETAIL_DIAGNOSTIC_HIGHFIVE_H

#include "highfive/H5DataSet.hpp"
#include "highfive/H5DataSpace.hpp"
#include "highfive/H5Easy.hpp"

#include "h5typewriter.h"
#include "h5file.h"

#include "diagnostic/diagnostic_manager.h"
#include "diagnostic/diagnostic_props.h"

#include "core/data/vecfield/vecfield_component.h"
#include "core/utilities/mpi_utils.h"
#include "core/utilities/types.h"


#if !defined(PHARE_DIAG_DOUBLES)
#error // PHARE_DIAG_DOUBLES not defined
#endif


namespace PHARE::diagnostic::h5
{
template<typename H5Writer>
class ElectromagDiagnosticWriter;
template<typename H5Writer>
class FluidDiagnosticWriter;
template<typename H5Writer>
class ParticlesDiagnosticWriter;



template<typename ModelView>
class Writer : public PHARE::diagnostic::IWriter
{
    using FloatType = std::conditional_t<PHARE_DIAG_DOUBLES, double, float>;

    static constexpr std::size_t timestamp_precision = 10;

public:
    using This       = Writer<ModelView>;
    using GridLayout = typename ModelView::GridLayout;
    using Attributes = typename ModelView::PatchProperties;

    static constexpr auto dimension   = GridLayout::dimension;
    static constexpr auto interpOrder = GridLayout::interp_order;
    static constexpr auto READ_WRITE  = HiFile::ReadWrite | HiFile::Create;

    // flush_never: disables manual file closing, but still occurrs via RAII
    static constexpr std::size_t flush_never = 0;

    template<typename Hierarchy, typename Model>
    Writer(Hierarchy& hier, Model& model, std::string const hifivePath,
           unsigned _flags /* = HiFile::ReadWrite | HiFile::Create | HiFile::Truncate */)
        : flags{_flags}
        , filePath_{hifivePath}
        , modelView_{hier, model}
    {
    }

    ~Writer() {}

    template<typename Hierarchy, typename Model>
    static decltype(auto) make_unique(Hierarchy& hier, Model& model,
                                      initializer::PHAREDict const& dict)
    {
        std::string filePath = dict["filePath"].template to<std::string>();
        unsigned flags       = READ_WRITE;
        if (dict.contains("mode") and dict["mode"].template to<std::string>() == "overwrite")
            flags |= HiFile::Truncate;
        return std::make_unique<This>(hier, model, filePath, flags);
    }


    void dump(std::vector<DiagnosticProperties*> const&, double current_timestamp) override;
    void dump_level(std::size_t level, std::vector<DiagnosticProperties*> const& diagnostics,
                    double timestamp) override;

    template<typename String>
    auto getDiagnosticWriterForType(String& type)
    {
        return writers.at(type);
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

    template<typename Type>
    static void createDataSet(HiFile& h5, std::string const& path, std::size_t size)
    {
        if constexpr (std::is_same_v<Type, double>) // force doubles for floats for storage
            This::createDatasetsPerMPI<FloatType>(h5, path, size);
        else
            This::createDatasetsPerMPI<Type>(h5, path, size);
    }



    template<typename Array, typename String>
    static void writeDataSet(HiFile& h5, String path, Array const* const array)
    {
        h5.getDataSet(path).write(array);
    }

    // global function when all path+key are the same
    template<typename Data>
    static void writeAttribute(HiFile& h5, std::string path, std::string key, Data const& value)
    {
        h5.getGroup(path)
            .template createAttribute<Data>(key, HighFive::DataSpace::From(value))
            .write(value);
    }

    template<typename Dict>
    static void writeAttributeDict(HiFile& h5, Dict dict, std::string path)
    {
        dict.visit([&](std::string const& key, const auto& val) {
            writeAttributesPerMPI(h5, path, key, val);
        });
    }

    // per MPI function where path differs per process
    template<typename Data>
    static void writeAttributesPerMPI(HiFile& h5, std::string path, std::string key,
                                      Data const& value);

    template<typename VecField>
    static void writeVecFieldAsDataset(HiFile& h5, std::string path, VecField& vecField)
    {
        for (auto& [id, type] : core::Components::componentMap)
            h5.getDataSet(path + "_" + id).write(vecField.getComponent(type).data());
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

    std::unordered_map<std::string, std::shared_ptr<H5TypeWriter<This>>> writers{
        {"fluid", make_writer<FluidDiagnosticWriter<This>>()},
        {"electromag", make_writer<ElectromagDiagnosticWriter<This>>()},
        {"particle", make_writer<ParticlesDiagnosticWriter<This>>()}};

    template<typename Writer>
    std::shared_ptr<H5TypeWriter<This>> make_writer()
    {
        return std::make_shared<Writer>(*this);
    }

    template<typename Type>
    static void createDatasetsPerMPI(HiFile& h5, std::string path, std::size_t dataSetSize);


    void initializeDatasets_(std::vector<DiagnosticProperties*> const& diagnotics);
    void writeDatasets_(std::vector<DiagnosticProperties*> const& diagnotics);

    Writer(const Writer&)             = delete;
    Writer(const Writer&&)            = delete;
    Writer& operator&(const Writer&)  = delete;
    Writer& operator&(const Writer&&) = delete;


    //  State of this class is controlled via "dump()"
    //  block public access to internal state
    friend class FluidDiagnosticWriter<This>;
    friend class ElectromagDiagnosticWriter<This>;
    friend class ParticlesDiagnosticWriter<This>;
    friend class H5TypeWriter<This>;

    // used by friends start
    std::string getPatchPathAddTimestamp(int iLevel, std::string globalCoords)
    {
        return getFullPatchPath(core::to_string_with_precision(timestamp_, timestamp_precision),
                                iLevel, globalCoords);
    }


    const auto& patchPath() const { return patchPath_; }
    // used by friends end
};



template<typename ModelView>
void Writer<ModelView>::dump(std::vector<DiagnosticProperties*> const& diagnostics,
                             double timestamp)
{
    timestamp_                     = timestamp;
    fileAttributes_["dimension"]   = dimension;
    fileAttributes_["interpOrder"] = interpOrder;
    fileAttributes_["layoutType"]  = modelView_.getLayoutTypeString();
    fileAttributes_["domain_box"]  = modelView_.domainBox();
    fileAttributes_["cell_width"]  = modelView_.cellWidth();
    fileAttributes_["origin"]      = modelView_.origin();

    for (auto* diagnostic : diagnostics)
        if (!file_flags.count(diagnostic->type + diagnostic->quantity))
            file_flags[diagnostic->type + diagnostic->quantity] = this->flags;

    initializeDatasets_(diagnostics);
    writeDatasets_(diagnostics);

    for (auto* diagnostic : diagnostics)
    {
        writers.at(diagnostic->type)->finalize(*diagnostic);
        // don't truncate past first dump
        file_flags[diagnostic->type + diagnostic->quantity] = READ_WRITE;
    }
}

template<typename ModelView>
void Writer<ModelView>::dump_level(std::size_t level,
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


namespace
{ // during attribute/dataset creation, we currently don't require the parents of the group to
  // exist, but create all that are missing - this used to exist in highfive
    inline std::string getParentName(const std::string& path)
    {
        std::size_t idx = path.find_last_of("/");
        if (idx == std::string::npos or idx == 0)
            return "/";
        return path.substr(0, idx);
    }

    inline void createGroupsToDataSet(HiFile& file, const std::string& path)
    {
        std::string group_name = getParentName(path);
        if (!file.exist(group_name))
            file.createGroup(group_name);
    }
}


/*
 * Communicate all dataset paths and sizes to all MPI process to allow each to create all
 * datasets independently. This is a requirement of HDF5.
 * in the case of a disparate number of datasets per MPI process, dataSetSize may be 0
 * such that the current process creates datasets for all other processes with non-zero
 * sizes. Recommended to use similar sized paths, if possible.
 */
template<typename ModelView>
template<typename Type>
void Writer<ModelView>::createDatasetsPerMPI(HiFile& h5file, std::string path,
                                             std::size_t dataSetSize)
{
    auto mpi_size = core::mpi::size();
    auto sizes    = core::mpi::collect(dataSetSize, mpi_size);
    auto paths    = core::mpi::collect(path, mpi_size);
    for (int i = 0; i < mpi_size; i++)
    {
        if (sizes[i] == 0)
            continue;
        createGroupsToDataSet(h5file, paths[i]);
        h5file.createDataSet<Type>(paths[i], HighFive::DataSpace(sizes[i]));
    }
}



/*
 * Communicate all attribute paths and values to all MPI process to allow each to create all
 * attributes independently. This is a requirement of HDF5.
 * in the case of a disparate number of attributes per MPI process, path may be an empty string
 * such that the current process creates attributes for all other processes with non-zero
 * sizes. Recommended to use similar sized paths, if possible. key is always assumed to the be
 * the same
 */
template<typename ModelView>
template<typename Data>
void Writer<ModelView>::writeAttributesPerMPI(HiFile& h5file, std::string path, std::string key,
                                              Data const& data)
{
    constexpr bool data_is_vector = core::is_std_vector_v<Data>;

    auto doAttribute = [&](auto node, auto const& _key, auto const& value) {
        if constexpr (data_is_vector)
            node.template createAttribute<typename Data::value_type>(
                    _key, HighFive::DataSpace(value.size()))
                .write(value.data());
        else
            node.template createAttribute<Data>(_key, HighFive::DataSpace::From(value))
                .write(value);
    };

    int mpi_size = core::mpi::size();
    auto values  = [&]() {
        if constexpr (data_is_vector)
            return core::mpi::collect_raw(data, mpi_size);
        else
            return core::mpi::collect(data, mpi_size);
    }();
    auto paths = core::mpi::collect(path, mpi_size);

    for (int i = 0; i < mpi_size; i++)
    {
        std::string const keyPath = paths[i] == "null" ? "" : paths[i];
        if (keyPath.empty())
            continue;

        if (h5file.exist(keyPath) && h5file.getObjectType(keyPath) == HighFive::ObjectType::Dataset)
        {
            if (!h5file.getDataSet(keyPath).hasAttribute(key))
                doAttribute(h5file.getDataSet(keyPath), key, values[i]);
        }
        else // group
        {
            createGroupsToDataSet(h5file, keyPath + "/dataset");
            if (!h5file.getGroup(keyPath).hasAttribute(key))
                doAttribute(h5file.getGroup(keyPath), key, values[i]);
        }
    }
}



template<typename ModelView>
void Writer<ModelView>::initializeDatasets_(std::vector<DiagnosticProperties*> const& diagnostics)
{
    std::size_t maxLocalLevel = 0;
    std::unordered_map<std::size_t, std::vector<std::string>> lvlPatchIDs;
    Attributes patchAttributes; // stores dataset info/size for synced MPI creation

    for (auto* diag : diagnostics)
        writers.at(diag->type)->createFiles(*diag);

    auto collectPatchAttributes = [&](GridLayout&, std::string patchID, std::size_t iLevel) {
        if (!lvlPatchIDs.count(iLevel))
            lvlPatchIDs.emplace(iLevel, std::vector<std::string>());

        lvlPatchIDs.at(iLevel).emplace_back(patchID);

        for (auto* diag : diagnostics)
        {
            writers.at(diag->type)->getDataSetInfo(*diag, iLevel, patchID, patchAttributes);
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
        writers.at(diagnostic->type)
            ->initDataSets(*diagnostic, lvlPatchIDs, patchAttributes, maxMPILevel);
    }
}



template<typename ModelView>
void Writer<ModelView>::writeDatasets_(std::vector<DiagnosticProperties*> const& diagnostics)
{
    std::unordered_map<std::size_t, std::vector<std::pair<std::string, Attributes>>>
        patchAttributes;

    std::size_t maxLocalLevel = 0;
    auto writePatch = [&](GridLayout& gridLayout, std::string patchID, std::size_t iLevel) {
        if (!patchAttributes.count(iLevel))
            patchAttributes.emplace(iLevel, std::vector<std::pair<std::string, Attributes>>{});
        patchPath_ = getPatchPathAddTimestamp(iLevel, patchID);
        patchAttributes[iLevel].emplace_back(patchID, modelView_.getPatchProperties(gridLayout));
        for (auto* diagnostic : diagnostics)
            writers.at(diagnostic->type)->write(*diagnostic);
        maxLocalLevel = iLevel;
    };

    modelView_.visitHierarchy(writePatch, minLevel, maxLevel);

    std::size_t maxMPILevel = core::mpi::max(maxLocalLevel);
    // sets empty vectors in case current process lacks patch on a level
    for (std::size_t lvl = minLevel; lvl <= maxMPILevel; lvl++)
        if (!patchAttributes.count(lvl))
            patchAttributes.emplace(lvl, std::vector<std::pair<std::string, Attributes>>{});

    for (auto* diagnostic : diagnostics)
        writers.at(diagnostic->type)
            ->writeAttributes(*diagnostic, fileAttributes_, patchAttributes, maxMPILevel);
}



} /* namespace PHARE::diagnostic::h5 */

#endif /* PHARE_DETAIL_DIAGNOSTIC_HIGHFIVE_H */
