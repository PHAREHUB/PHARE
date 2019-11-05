#ifndef PHARE_AMR_DIAGNOSTIC_SAMRAI_HIGHFIVE_H
#define PHARE_AMR_DIAGNOSTIC_SAMRAI_HIGHFIVE_H

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
    using GridLayout                   = typename ModelView::GridLayout;
    static constexpr auto dimension    = ModelView::dimension;
    static constexpr auto interp_order = GridLayout::interp_order;

    HighFiveDiagnostic(ModelView& modelView, std::string const hifivePath)
        : hi5_{hifivePath}
        , modelView_{modelView}
    {
    }

    void dump(std::vector<DiagnosticDAO*> const&);

    HighFive::File& file() { return hi5_.file(); }

    template<typename String>
    auto getDiagnosticWriterForType(String& writer)
    {
        return writers.at(writer);
    }

private:
    HighFiveFile hi5_;
    ModelView& modelView_;
    std::string patchPath_; // is passed around as "virtual write()" has no parameters

    std::unordered_map<std::string, std::shared_ptr<PHARE::DiagnosticWriter>> writers{
        {"electromag", make_writer<ElectromagDiagnosticWriter>()},
        {"particles", make_writer<ParticlesDiagnosticWriter>()}};

    template<typename Writer>
    std::shared_ptr<PHARE::DiagnosticWriter> make_writer()
    {
        return std::make_shared<Writer>(*this);
    }

    template<typename String>
    auto getOrCreateGroup(String const& path)
    {
        if (path != std::string("/") && !hi5_.file_.exist(path))
            hi5_.file_.createGroup(path);
        return hi5_.file_.getGroup(path);
    }

    template<typename Type, typename Size>
    auto getOrCreateDataSet(std::string const& path, Size size)
    {
        if (!hi5_.file_.exist(path))
            return createDataSet<Type>(path, size);
        return hi5_.file_.getDataSet(path);
    }

    auto getDataSet(std::string const& path) { return hi5_.file_.getDataSet(path); }

    template<typename Type, typename Size>
    auto createDataSet(std::string const& path, Size size)
    {
        H5Easy::detail::createGroupsToDataSet(hi5_.file_, path);
        if constexpr (std::is_same_v<Size, size_t>)
            hi5_.file_.template createDataSet<Type>(path, HighFive::DataSpace(size));
        else
            hi5_.file_.template createDataSet<Type>(path, HighFive::DataSpace::From(size));
        return hi5_.file_.getDataSet(path);
    }

    template<typename Dataset, typename Array>
    void writeDataSetPart(Dataset dataSet, size_t start, size_t size, Array const& array)
    {
        dataSet.select({start}, {size}).write(array);
    }

    template<typename Array, typename String>
    void writeDataSet(String path, Array const* const array, size_t size)
    {
        createDataSet<Array>(path, size).write(array);
    }

    template<typename String, typename Data>
    void writeAttribute(String path, std::string const& key, Data const& value)
    {
        getOrCreateGroup(path)
            .template createAttribute<Data>(key, HighFive::DataSpace::From(value))
            .write(value);
    }

    template<typename Dict>
    void writeDict(Dict, std::string const&);
    template<typename Dict> // template String causes internal compiler error in GCC 8.2
    void writeDict(Dict dict, std::string const&& str)
    {
        writeDict(dict, str);
    }

    HighFiveDiagnostic(const HighFiveDiagnostic&)             = delete;
    HighFiveDiagnostic(const HighFiveDiagnostic&&)            = delete;
    HighFiveDiagnostic& operator&(const HighFiveDiagnostic&)  = delete;
    HighFiveDiagnostic& operator&(const HighFiveDiagnostic&&) = delete;

    class Hi5DiagnosticWriter;
    class ElectromagDiagnosticWriter; // : public Hi5DiagnosticWriter
    class ParticlesDiagnosticWriter;  // : public Hi5DiagnosticWriter
    class FluidDiagnosticWriter;      // : public Hi5DiagnosticWriter
};                                    // namespace PHARE

/*TO DO
  investigate level > 0 for MPI
  finalise HDF5 Path format
*/
template<typename ModelView>
void HighFiveDiagnostic<ModelView>::dump(std::vector<DiagnosticDAO*> const& diagnostics)
{
    auto levelPath = [](auto idx) { return "/t#/pl" + std::to_string(idx); };

    writeAttribute("/", "dim", dimension);
    writeAttribute("/", "interpOrder", interp_order);

    /*TODO
      add time/iterations
    */
    modelView_.visitHierarchy([&](GridLayout& gridLayout, std::string patchID, int iLevel) {
        patchPath_ = levelPath(iLevel) + "/p" + patchID;
        getOrCreateGroup(patchPath_);
        for (auto* diagnostic : diagnostics)
        {
            writers.at(diagnostic->type)->write(*diagnostic);
        }
        writeDict(modelView_.getPatchAttributes(gridLayout), patchPath_);
    });
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
            outer.writeDataSetPart(dataset, start, 1, value); // not array, write 1 value
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

    for (auto& pop : outer.modelView_.getIons()) // bulkV
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


/*TODO
  finish
*/
template<typename ModelView>
void HighFiveDiagnostic<ModelView>::ElectromagDiagnosticWriter::write([
    [maybe_unused]] DiagnosticDAO& diagnostic)
{
    auto& outer = this->outer_;
    for (auto& fields : outer.modelView_.getElectromagFields())
    {
        for (auto& field : fields)
        {
            std::string fieldPath(outer.patchPath_ + "/" + field->id);
            outer.writeDataSet(fieldPath, field->data, field->size);
        }
    };
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
}

// turns a dict of std::map<std::string, T> to hdf5 attributes
template<typename ModelView>
template<typename Dict>
void HighFiveDiagnostic<ModelView>::writeDict(Dict dict, std::string const& path)
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
            // static_assert fails without if ! constexpr all possible types
            //   regardless of what it actually is.
            throw std::runtime_error(std::string("Expecting map<string, T> got ")
                                     + typeid(Map_t).name());
    };
    std::visit(visitor, dict.data);
}

} /*namespace PHARE*/

#endif /*PHARE_AMR_DIAGNOSTIC_SAMRAI_HIGHFIVE_H*/
