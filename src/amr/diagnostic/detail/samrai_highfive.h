#ifndef PHARE_AMR_DIAGNOSTIC_SAMRAI_HIGHFIVE_H
#define PHARE_AMR_DIAGNOSTIC_SAMRAI_HIGHFIVE_H

#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5Easy.hpp>

#include "diagnostic/samrai_diagnostic.h"

namespace PHARE
{
namespace hi5
{
    struct Diagnostic
    {
        HighFive::File& file_;
        PHARE::diagnostic::Mode mode_ = PHARE::diagnostic::Mode::LIGHT;

        Diagnostic(const Diagnostic&)             = delete;
        Diagnostic(const Diagnostic&&)            = delete;
        Diagnostic& operator&(const Diagnostic&)  = delete;
        Diagnostic& operator&(const Diagnostic&&) = delete;

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
    };
} // namespace hi5
template<typename Model>
class SamraiHighFiveDiagnostic : public SamraiDiagnostic<Model>
{
public:
    using SamraiDiagnostic<Model>::modelView_;

    using Hierarchy   = SAMRAI::hier::PatchHierarchy;
    using PatchLevel  = SAMRAI::hier::PatchLevel;
    using Diagnostics = typename DiagnosticsManager<Model>::DiagnosticWritingList;

    static constexpr auto dimensions = Model::dimension;

    SamraiHighFiveDiagnostic(Hierarchy& h, Model& r, hi5::Diagnostic& hifive)
        : PHARE::SamraiDiagnostic<Model>(h, r)
        , hi5_(hifive)
    {
    }

    void dump(Diagnostics const&);

    template<typename String>
    auto getWriter(String& writer)
    {
        return writers.at(writer);
    }

private:
    hi5::Diagnostic& hi5_;
    std::string patchPath_; // is passed around as "virtual write()" has no parameters

    std::unordered_map<std::string, std::shared_ptr<PHARE::DiagnosticWriter>> writers{
        {"electromag", make_writer<ElectromagDiagnosticWriter>()},
        {"particles", make_writer<ParticlesDiagnosticWriter>()}};

    template<typename Writer>
    std::shared_ptr<PHARE::DiagnosticWriter> make_writer()
    {
        return std::make_shared<Writer>(*this);
    }

    void write(PatchLevel&, std::string&&, Diagnostics const&);

    template<typename String>
    auto getOrCreateGroup(String const& path)
    {
        if (!hi5_.file_.exist(path))
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

    auto getDataSet(std::string const& path)
    {
        return hi5_.file_.getDataSet(path);
    }

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

    template<typename Data>
    void writeData(std::string const&& path, Data const& value)
    {
        createDataSet<Data>(path, value).write(value);
    }

    template<typename Dataset, typename Array>
    void writeDataSetPart(Dataset dataSet, size_t start, size_t end, Array const* const array)
    {
        dataSet.select({start, end}).write(array);
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

    SamraiHighFiveDiagnostic(const SamraiHighFiveDiagnostic&)             = delete;
    SamraiHighFiveDiagnostic(const SamraiHighFiveDiagnostic&&)            = delete;
    SamraiHighFiveDiagnostic& operator&(const SamraiHighFiveDiagnostic&)  = delete;
    SamraiHighFiveDiagnostic& operator&(const SamraiHighFiveDiagnostic&&) = delete;

    class Hi5DiagnosticWriter;
    class ElectromagDiagnosticWriter; // : public Hi5DiagnosticWriter
    class ParticlesDiagnosticWriter;  // : public Hi5DiagnosticWriter
    class FluidDiagnosticWriter;      // : public Hi5DiagnosticWriter
};

/*TO DO
  investigate level > 0 for MPI
  finalise HDF5 Path format
*/
template<typename Model>
void SamraiHighFiveDiagnostic<Model>::dump(Diagnostics const& diagnostics)
{
    auto levelPath = [](auto idx) { return "/t#/pl" + std::to_string(idx); };

    writeAttribute("/", "dim", dim);
    writeAttribute("/", "interpOrder", interpOrder);

    write(*this->hierarchy_.getPatchLevel(0), levelPath(0), diagnostics);

    /* CAUSES APP TO HANG ON EXIT WITH MPI */
    /*if (hi5_.mode_ == PHARE::diagnostic::Mode::FULL)
    {
        for (int iLevel = 1; iLevel < this->hierarchy_.getNumberOfLevels(); iLevel++)
        {
            write(this->hierarchy_.getPatchLevel(iLevel), levelPath(iLevel),
                diagnostics);
        }
    }*/
}

/*TODO
  add time/iterations
*/
template<typename Model>
void SamraiHighFiveDiagnostic<Model>::write(PatchLevel& level, std::string&& path,
                                            Diagnostics const& diagnostics)
{
    size_t patch_idx = 0;
    for (auto& patch : level)
    {
        auto guardedGrid = modelView_.guardedGrid(*patch);
        std::stringstream patchPath;
        patchPath << path << "/p" << patch_idx;
        patchPath_ = patchPath.str();
        getOrCreateGroup(patchPath_);

        for (auto& [diagnostic, writer] : diagnostics)
        {
            writer->write(diagnostic.get());
        }

        writeDict(modelView_.getPatchAttributes(guardedGrid), patchPath_);
        patch_idx++;
    }
}


template<typename Model>
class SamraiHighFiveDiagnostic<Model>::Hi5DiagnosticWriter : public PHARE::DiagnosticWriter
{
public:
    Hi5DiagnosticWriter(SamraiHighFiveDiagnostic& _outer)
        : outer_(_outer)
    {
    }

protected:
    SamraiHighFiveDiagnostic& outer_;
};

template<typename Model>
class SamraiHighFiveDiagnostic<Model>::ParticlesDiagnosticWriter : public Hi5DiagnosticWriter
{
public:
    ParticlesDiagnosticWriter(SamraiHighFiveDiagnostic& _outer)
        : Hi5DiagnosticWriter(_outer)
    {
    }
    void write(Diagnostic&) override;
    void compute(Diagnostic&) override{};
};
/*TODO
  finish
*/
template<typename Model>
void SamraiHighFiveDiagnostic<Model>::ParticlesDiagnosticWriter::write([
    [maybe_unused]] Diagnostic& diagnostic)
{
    auto& outer = this->outer_;
    auto inc    = [](auto& in, auto& idx) {
        idx++;
        return in;
    };
    auto createDataSetLam = [&outer](auto&& path, auto size, auto const& value) {
        using Val = std::decay_t<decltype(value)>;
        if constexpr (is_std_array<Val, dimensions>::value || is_std_array<Val, 3>::value)
            return outer.template createDataSet<typename Val::value_type>(path, size);
        else
            return outer.template createDataSet<Val>(path, size);
    };
    auto writeDatSetLam = [&outer](auto& datasets, auto& idx, auto& part_idx, auto& sizes,
                                   auto const& value) {
        using Val  = std::decay_t<decltype(value)>;
        auto start = idx * sizes[part_idx];
        if constexpr (is_std_array<Val, dimensions>::value || is_std_array<Val, 3>::value)
            outer.writeDataSetPart(datasets[part_idx], start, start + sizes[part_idx] - 1,
                                   &value[0]);
        else
            outer.writeDataSetPart(datasets[part_idx], start, start + sizes[part_idx] - 1, &value);
        part_idx++;
    };

    auto writeParticles = [&](auto& particles, auto&& path) {
        if (!particles.size())
            return;
        auto packer = outer.modelView_.getParticlePacker(particles);

        std::vector<HighFive::DataSet> datasets;
        auto sizes = packer.sizes();
        size_t idx = 0;
        std::apply(
            [&](auto&... args) {
                (inc(datasets.emplace_back(createDataSetLam(path + packer.keys()[idx],
                                                            sizes[idx] * particles.size(), args)),
                     idx),
                 ...);
            },
            packer.first());

        idx = 0;
        while (packer.hasNext())
        {
            size_t part_idx = 0;
            std::apply(
                [&](auto&... args) { (writeDatSetLam(datasets, idx, part_idx, sizes, args), ...); },
                packer.next());
            idx++;
        }
    };

    size_t pop_idx = 0;
    for (auto& pop : outer.modelView_.getParticles())
    {
        std::stringstream particlePath;
        particlePath << outer.patchPath_ << "/ions/pop/" << pop_idx++ << "/"; // bulkV

        std::stringstream domainPath;
        domainPath << particlePath.str() << "domain/";
        writeParticles(pop.domainParticles(), domainPath.str());
    }
}

template<typename Model>
class SamraiHighFiveDiagnostic<Model>::ElectromagDiagnosticWriter : public Hi5DiagnosticWriter
{
public:
    ElectromagDiagnosticWriter(SamraiHighFiveDiagnostic& _outer)
        : Hi5DiagnosticWriter(_outer)
    {
    }
    void write(Diagnostic&) override;
    void compute(Diagnostic&) override{};
};
/*TODO
  finish
*/
template<typename Model>
void SamraiHighFiveDiagnostic<Model>::ElectromagDiagnosticWriter::write([
    [maybe_unused]] Diagnostic& diagnostic)
{
    auto& outer = this->outer_;
    for (auto& fields : outer.modelView_.getElectromagFields())
    {
        for (auto& field : fields)
        {
            std::stringstream fieldPath;
            fieldPath << outer.patchPath_ << "/" << field->id;
            outer.writeDataSet(fieldPath.str(), field->data, field->size);
        }
    };
}

template<typename Model>
class SamraiHighFiveDiagnostic<Model>::FluidDiagnosticWriter : public Hi5DiagnosticWriter
{
public:
    FluidDiagnosticWriter(SamraiHighFiveDiagnostic& _outer)
        : Hi5DiagnosticWriter(_outer)
    {
    }
    void write(Diagnostic&) override;
    void compute(Diagnostic&) override{};
};

template<typename Model>
void SamraiHighFiveDiagnostic<Model>::FluidDiagnosticWriter::write([
    [maybe_unused]] Diagnostic& diagnostic)
{
}

// turns a dict of std::map<std::string, T> to hdf5 attributes
template<typename Model>
template<typename Dict>
void SamraiHighFiveDiagnostic<Model>::writeDict(Dict dict, std::string const& path)
{
    using dict_map_t = typename Dict::map_t;
    auto visitor     = [&](auto&& map) {
        using Map_t = std::decay_t<decltype(map)>;
        if constexpr (std::is_same_v<Map_t, dict_map_t>)
        {
            for (auto& pair : map)
            {
                std::visit(
                    [&](auto&& val) {
                        using Val = std::decay_t<decltype(val)>;
                        if constexpr (is_std_array<Val, dimensions>::value
                                      || is_std_array<Val, 3>::value)
                        {
                            writeDataSet(path + pair.first, &val[0], val.size());
                        }
                        else if constexpr (
                            !std::is_same_v<Val,
                                            cppdict::NoValue> && !std::is_same_v<Val, dict_map_t>)
                        {
                            writeAttribute(path, pair.first, val);
                        }
                        else
                            throw std::runtime_error(std::string("Expecting Writable value got ")
                                                     + typeid(Map_t).name());
                    },
                    pair.second.get()->data);
            }
        }
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
