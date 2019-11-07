#ifndef PHARE_AMR_DIAGNOSTIC_SAMRAI_HIGHFIVE_H
#define PHARE_AMR_DIAGNOSTIC_SAMRAI_HIGHFIVE_H

#include "kul/dbg.hpp"

#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5Easy.hpp>

#include "diagnostic_manager.h"

namespace PHARE
{
struct HighFiveFile
{
    HighFive::File file_;
    PHARE::diagnostic::Mode mode_ = PHARE::diagnostic::Mode::LIGHT;

    static auto createHighFiveFile(std::string const path)
    {
        return HighFive::File
        {
            path, HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Truncate
#if defined(H5_HAVE_PARALLEL)
                ,
                HighFive::MPIOFileDriver(MPI_COMM_WORLD, MPI_INFO_NULL)
#endif
        };
    }

    HighFiveFile(std::string const path,
                 PHARE::diagnostic::Mode mode = PHARE::diagnostic::Mode::LIGHT)
        : file_{createHighFiveFile(path)}
        , mode_{mode}
    {
    }
    ~HighFiveFile() {}

    HighFive::File& file() { return file_; }

    HighFiveFile(const HighFiveFile&)             = delete;
    HighFiveFile(const HighFiveFile&&)            = delete;
    HighFiveFile& operator&(const HighFiveFile&)  = delete;
    HighFiveFile& operator&(const HighFiveFile&&) = delete;
};


template<typename ModelView>
class HighFiveDiagnostic
{
public:
    using GridLayout                  = typename ModelView::GridLayout;
    static constexpr auto dimension   = ModelView::dimension;
    static constexpr auto interpOrder = GridLayout::interp_order;

    HighFiveDiagnostic(ModelView& modelView, std::string const hifivePath)
        : hi5_{hifivePath}
        , modelView_{modelView}
    {
    }

    ~HighFiveDiagnostic() {}

    void dump(std::vector<DiagnosticDAO*> const&);

    HighFive::File& file() { return hi5_.file(); }

    template<typename String>
    auto getDiagnosticWriterForType(String& writer)
    {
        return writers.at(writer);
    }

    std::string getPatchPath([[maybe_unused]] std::string time, int iLevel,
                             std::string globalCoords)
    {
        return "/t#/pl" + std::to_string(iLevel) + "/p" + globalCoords;
    }

private:
    HighFiveFile hi5_;
    ModelView& modelView_;
    std::string patchPath_; // is passed around as "virtual write()" has no parameters

    std::unordered_map<std::string, std::shared_ptr<PHARE::DiagnosticWriter>> writers{
        {"fluid", make_writer<FluidDiagnosticWriter>()},
        {"electromag", make_writer<ElectromagDiagnosticWriter>()},
        {"particles", make_writer<ParticlesDiagnosticWriter>()}};

    template<typename Writer>
    std::shared_ptr<PHARE::DiagnosticWriter> make_writer()
    {
        return std::make_shared<Writer>(*this);
    }

    template<typename Type>
    auto createDataSet(std::string const& path, size_t size)
    {
        this->template createDatasetsPerMPI<Type>(path, size);
        return hi5_.file_.getDataSet(path);
    }

    template<typename Type>
    void createDatasetsPerMPI(std::string path, size_t dataSetSize);

    template<typename Dataset, typename Array>
    void writeDataSetPart(Dataset dataSet, size_t start, size_t size, Array const& array)
    {
        dataSet.select({start}, {size}).write(array);
    }

    template<typename Array, typename String>
    void writeNewDataSet(String path, Array const* const array, size_t size)
    {
        createDataSet<Array>(path, size).write(array);
    }

    template<typename String, typename Data>
    void writeAttribute(String path, std::string const& key, Data const& value)
    {
        hi5_.file_.getGroup(path)
            .template createAttribute<Data>(key, HighFive::DataSpace::From(value))
            .write(value);
    }


    template<typename VecField>
    void writeVecFieldAsDataset(std::string path, VecField& vecField)
    {
        for (auto& [id, type] : Components::componentMap)
        {
            auto& field = vecField.getComponent(type);
            std::string dataSetPath{path + "/" + id};
            writeNewDataSet(dataSetPath, field.data(), field.size());
        }
    }

    template<typename Dict>
    void writeDict(Dict, std::string);

    HighFiveDiagnostic(const HighFiveDiagnostic&)             = delete;
    HighFiveDiagnostic(const HighFiveDiagnostic&&)            = delete;
    HighFiveDiagnostic& operator&(const HighFiveDiagnostic&)  = delete;
    HighFiveDiagnostic& operator&(const HighFiveDiagnostic&&) = delete;

    class Hi5DiagnosticWriter;
    class ElectromagDiagnosticWriter; /* : public Hi5DiagnosticWriter*/
    class ParticlesDiagnosticWriter;  /* : public Hi5DiagnosticWriter*/
    class FluidDiagnosticWriter;      /* : public Hi5DiagnosticWriter*/
};

/*TO DO
 * investigate level > 0 for MPI
 *  finalise HDF5 Path format
 */
template<typename ModelView>
void HighFiveDiagnostic<ModelView>::dump(std::vector<DiagnosticDAO*> const& diagnostics)
{
    writeAttribute("/", "dim", dimension);
    writeAttribute("/", "interpOrder", interpOrder);

    /*TODO
      add time/iterations
    */
    auto visitPatch = [&](GridLayout& gridLayout, std::string patchID, int iLevel) {
        patchPath_ = getPatchPath("time", 0, patchID);
        for (auto* diagnostic : diagnostics)
        {
            writers.at(diagnostic->type)->write(*diagnostic);
        }
        writeDict(modelView_.getPatchAttributes(gridLayout), patchPath_);
    };

    modelView_.visitHierarchy(visitPatch);
}

template<typename ModelView>
class HighFiveDiagnostic<ModelView>::Hi5DiagnosticWriter : public PHARE::DiagnosticWriter
{
public:
    Hi5DiagnosticWriter(HighFiveDiagnostic& _outer)
        : outer_(_outer)
    {
    }

protected:
    HighFiveDiagnostic& outer_;
};

template<typename ModelView>
class HighFiveDiagnostic<ModelView>::ParticlesDiagnosticWriter : public Hi5DiagnosticWriter
{
public:
    ParticlesDiagnosticWriter(HighFiveDiagnostic& _outer)
        : Hi5DiagnosticWriter(_outer)
    {
    }
    void write(DiagnosticDAO&) override;
    void compute(DiagnosticDAO&) override{};
};

template<typename T, std::size_t dimension>
inline constexpr auto is_array_dataset
    = (core::is_std_array_v<T, dimension> || core::is_std_array_v<T, 3>);

template<typename ModelView>
void HighFiveDiagnostic<ModelView>::ParticlesDiagnosticWriter::write([
    [maybe_unused]] DiagnosticDAO& diagnostic)
{
    auto& outer = this->outer_;

    auto createDataSet = [&outer](auto&& path, auto size, auto const& value) {
        using ValueType = std::decay_t<decltype(value)>;
        if constexpr (is_array_dataset<ValueType, dimension>)
            return outer.template createDataSet<typename ValueType::value_type>(
                path, size * value.size());
        else
            return outer.template createDataSet<ValueType>(path, size);
    };

    auto writeDatSet = [&outer](auto& dataset, auto& start, auto const& value) {
        using ValueType = std::decay_t<decltype(value)>;
        if constexpr (is_array_dataset<ValueType, dimension>)
            outer.writeDataSetPart(dataset, start * value.size(), value.size(), value);
        else
            outer.writeDataSetPart(dataset, start, 1, value); /*not array, write 1 value*/
    };

    auto writeParticles = [&](auto path, auto& particles) {
        if (!particles.size())
            return;

        size_t part_idx = 0;
        std::vector<HighFive::DataSet> datasets;
        auto packer = outer.modelView_.getParticlePacker(particles);

        std::apply(
            [&](auto&... args) {
                (
                    [&]() {
                        datasets.emplace_back(
                            createDataSet(path + packer.keys()[part_idx], particles.size(), args));
                        part_idx++;
                    }(),
                    ...);
            },
            packer.first());

        size_t idx = 0;
        while (packer.hasNext())
        {
            part_idx = 0;
            std::apply(
                [&](auto&... args) {
                    (
                        [&]() {
                            writeDatSet(datasets[part_idx], idx, args);
                            part_idx++;
                        }(),
                        ...);
                },
                packer.next());
            idx++;
        }
    };

    for (auto& pop : outer.modelView_.getIons())
    {
        std::string path(outer.patchPath_ + "/ions/pop/" + pop.name() + "/");
        writeParticles(path + "domain/", pop.domainParticles());
        writeParticles(path + "lvlGhost/", pop.levelGhostParticles());
        writeParticles(path + "patchGhost/", pop.patchGhostParticles());
    }
}

template<typename ModelView>
class HighFiveDiagnostic<ModelView>::ElectromagDiagnosticWriter : public Hi5DiagnosticWriter
{
public:
    ElectromagDiagnosticWriter(HighFiveDiagnostic& _outer)
        : Hi5DiagnosticWriter(_outer)
    {
    }
    void write(DiagnosticDAO&) override;
    void compute(DiagnosticDAO&) override{};
};


template<typename ModelView>
void HighFiveDiagnostic<ModelView>::ElectromagDiagnosticWriter::write([
    [maybe_unused]] DiagnosticDAO& diagnostic)
{
    auto& outer = this->outer_;

    for (auto* vecField : outer.modelView_.getElectromagFields())
    {
        outer.writeVecFieldAsDataset(outer.patchPath_ + "/" + vecField->name(), *vecField);
    }
}

template<typename ModelView>
class HighFiveDiagnostic<ModelView>::FluidDiagnosticWriter : public Hi5DiagnosticWriter
{
public:
    FluidDiagnosticWriter(HighFiveDiagnostic& _outer)
        : Hi5DiagnosticWriter(_outer)
    {
    }
    void write(DiagnosticDAO&) override;
    void compute(DiagnosticDAO&) override{};
};

template<typename ModelView>
void HighFiveDiagnostic<ModelView>::FluidDiagnosticWriter::write([
    [maybe_unused]] DiagnosticDAO& diagnostic)
{
    auto& outer = this->outer_;
    auto& ions  = outer.modelView_.getIons();
    std::string path(outer.patchPath_ + "/ions/");

    for (auto& pop : outer.modelView_.getIons())
    {
        std::string popPath(path + "pop/" + pop.name() + "/");
        auto& density = pop.density();
        outer.writeNewDataSet(popPath + "density", density.data(), density.size());
        outer.writeVecFieldAsDataset(popPath + "flux", pop.flux());
    }

    auto& density = ions.density();
    outer.writeNewDataSet(path + "density", density.data(), density.size());
    outer.writeVecFieldAsDataset(path + "bulkVelocity", ions.velocity());
}

/*
 * turns a dict of std::map<std::string, T> to hdf5 attributes
 */
template<typename ModelView>
template<typename Dict>
void HighFiveDiagnostic<ModelView>::writeDict(Dict dict, std::string path)
{
    using dict_map_t = typename Dict::map_t;
    auto visitor     = [&](auto&& map) {
        using Map_t = std::decay_t<decltype(map)>;
        if constexpr (std::is_same_v<Map_t, dict_map_t>)
            for (auto& pair : map)
                std::visit(
                    [&](auto&& val) {
                        using Val = std::decay_t<decltype(val)>;
                        if constexpr (core::is_dict_leaf<Val, Dict>::value)
                            writeAttribute(path, pair.first, val);
                        else
                            throw std::runtime_error(std::string("Expecting Writable value got ")
                                                     + typeid(Map_t).name());
                    },
                    pair.second.get()->data);
        else
            /* static_assert fails without if ! constexpr all possible types
             * regardless of what it actually is.
             */
            throw std::runtime_error(std::string("Expecting map<string, T> got ")
                                     + typeid(Map_t).name());
    };
    std::visit(visitor, dict.data);
}


/*
 * Communicate all dataset paths and sizes to all MPI process to allow each to create all datasets
 * independently.  This is a requirement of HDF5. all paths must be the same length
 */
template<typename ModelView>
template<typename Type>
void HighFiveDiagnostic<ModelView>::createDatasetsPerMPI(std::string path, size_t dataSetSize)
{
    int mpi_size;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    uint64_t pathSize = path.size();
    uint64_t datasetSizes[mpi_size];
    MPI_Allgather(&dataSetSize, 1, MPI_UINT64_T, datasetSizes, 1, MPI_UINT64_T, MPI_COMM_WORLD);

    char chars[pathSize * mpi_size];
    MPI_Allgather(path.c_str(), path.size(), MPI_CHAR, &chars, path.size(), MPI_CHAR,
                  MPI_COMM_WORLD);

    for (uint64_t i = 0; i < mpi_size; i++)
    {
        std::string datasetPath{&chars[pathSize * i], pathSize};
        H5Easy::detail::createGroupsToDataSet(hi5_.file_, datasetPath);
        hi5_.file_.createDataSet<Type>(datasetPath, HighFive::DataSpace(datasetSizes[i]));
    }
}

} /* namespace PHARE */

#endif /* PHARE_AMR_DIAGNOSTIC_SAMRAI_HIGHFIVE_H */
