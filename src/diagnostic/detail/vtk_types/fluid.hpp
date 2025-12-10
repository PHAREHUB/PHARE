#ifndef PHARE_DIAGNOSTIC_DETAIL_VTK_TYPES_FLUID_HPP
#define PHARE_DIAGNOSTIC_DETAIL_VTK_TYPES_FLUID_HPP

#include "core/logger.hpp"

#include "amr/physical_models/mhd_model.hpp"
#include "amr/physical_models/hybrid_model.hpp"
#include "diagnostic/detail/vtkh5_type_writer.hpp"

#include <stdexcept>
#include <string>
#include <optional>
#include <unordered_map>

namespace PHARE::diagnostic::vtkh5
{

template<typename H5Writer>
class FluidDiagnosticWriter : public H5TypeWriter<H5Writer>
{
    using Super              = H5TypeWriter<H5Writer>;
    using VTKFileWriter      = Super::VTKFileWriter;
    using VTKFileInitializer = Super::VTKFileInitializer;

public:
    FluidDiagnosticWriter(H5Writer& h5Writer)
        : Super{h5Writer}
    {
    }

    void setup(DiagnosticProperties&) override;
    void write(DiagnosticProperties&) override;
    void compute(DiagnosticProperties&) override {};

private:
    struct Info
    {
        std::vector<std::size_t> offset_per_level = std::vector<std::size_t>(amr::MAX_LEVEL);
    };

    struct HybridFluidInitializer
    {
        std::optional<std::size_t> operator()(auto const& level);

        FluidDiagnosticWriter* writer;
        DiagnosticProperties& diagnostic;
        VTKFileInitializer& file_initializer;
    };

    struct MhdFluidInitializer
    {
        std::optional<std::size_t> operator()(auto const& level);

        FluidDiagnosticWriter* writer;
        DiagnosticProperties& diagnostic;
        VTKFileInitializer& file_initializer;
    };

    struct HybridFluidWriter
    {
        void operator()(auto const& layout);

        FluidDiagnosticWriter* writer;
        DiagnosticProperties& diagnostic;
        VTKFileWriter& file_writer;
    };

    struct MhdFluidWriter
    {
        void operator()(auto const& layout);

        FluidDiagnosticWriter* writer;
        DiagnosticProperties& diagnostic;
        VTKFileWriter& file_writer;
    };

    auto static isActiveDiag(DiagnosticProperties const& diagnostic, std::string const& tree,
                             std::string const& var)
    {
        return diagnostic.quantity == tree + var;
    };

    std::unordered_map<std::string, Info> mem;
};



template<typename H5Writer>
std::optional<std::size_t>
FluidDiagnosticWriter<H5Writer>::HybridFluidInitializer::operator()(auto const& level)
{
    auto& modelView = writer->h5Writer_.modelView();
    auto& ions      = modelView.getIons();
    std::string const tree{"/ions/"};

    if (isActiveDiag(diagnostic, tree, "charge_density"))
        return file_initializer.initFieldFileLevel(level);
    if (isActiveDiag(diagnostic, tree, "mass_density"))
        return file_initializer.initFieldFileLevel(level);
    if (isActiveDiag(diagnostic, tree, "bulkVelocity"))
        return file_initializer.template initTensorFieldFileLevel<1>(level);

    for (auto& pop : ions)
    {
        auto const pop_tree = tree + "pop/" + pop.name() + "/";
        if (isActiveDiag(diagnostic, pop_tree, "density"))
            return file_initializer.initFieldFileLevel(level);
        else if (isActiveDiag(diagnostic, pop_tree, "charge_density"))
            return file_initializer.initFieldFileLevel(level);
        else if (isActiveDiag(diagnostic, pop_tree, "flux"))
            return file_initializer.template initTensorFieldFileLevel<1>(level);
    }

    return std::nullopt;
}



template<typename H5Writer>
std::optional<std::size_t>
FluidDiagnosticWriter<H5Writer>::MhdFluidInitializer::operator()(auto const& level)
{
    return std::nullopt;
}


template<typename H5Writer>
void FluidDiagnosticWriter<H5Writer>::setup(DiagnosticProperties& diagnostic)
{
    PHARE_LOG_SCOPE(3, "FluidDiagnosticWriter<H5Writer>::setup");

    using Model_t   = H5Writer::ModelView::Model_t;
    auto& modelView = this->h5Writer_.modelView();

    VTKFileInitializer initializer{diagnostic, this};

    if (mem.count(diagnostic.quantity) == 0)
        mem.try_emplace(diagnostic.quantity);
    auto& info = mem[diagnostic.quantity];

    auto const init = [&](auto const& level) -> std::optional<std::size_t> {
        //

        if constexpr (solver::is_hybrid_model_v<Model_t>)
            if (auto ret = HybridFluidInitializer{this, diagnostic, initializer}(level))
                return ret;

        if constexpr (solver::is_mhd_model_v<Model_t>)
            if (auto ret = MhdFluidInitializer{this, diagnostic, initializer}(level))
                return ret;

        return std::nullopt;
    };

    modelView.onLevels(
        [&](auto const& level) {
            PHARE_LOG_SCOPE(3, "FluidDiagnosticWriter<H5Writer>::setup");

            auto const ilvl = level.getLevelNumber();
            initializer.initFileLevel(ilvl);
            if (auto const offset = init(level))
                info.offset_per_level[ilvl] = *offset;
        },
        this->h5Writer_.minLevel, this->h5Writer_.maxLevel);
}



template<typename H5Writer>
void FluidDiagnosticWriter<H5Writer>::HybridFluidWriter::operator()(auto const& layout)
{
    auto& modelView = writer->h5Writer_.modelView();
    auto& ions      = modelView.getIons();
    std::string const tree{"/ions/"};

    if (isActiveDiag(diagnostic, tree, "charge_density"))
        file_writer.writeField(ions.chargeDensity(), layout);
    else if (isActiveDiag(diagnostic, tree, "mass_density"))
        file_writer.writeField(ions.massDensity(), layout);
    else if (isActiveDiag(diagnostic, tree, "bulkVelocity"))
        file_writer.template writeTensorField<1>(ions.velocity(), layout);
    else
    {
        for (auto& pop : ions)
        {
            auto const pop_tree = tree + "pop/" + pop.name() + "/";
            if (isActiveDiag(diagnostic, pop_tree, "density"))
                file_writer.writeField(pop.particleDensity(), layout);
            else if (isActiveDiag(diagnostic, pop_tree, "charge_density"))
                file_writer.writeField(pop.chargeDensity(), layout);
            else if (isActiveDiag(diagnostic, pop_tree, "flux"))
                file_writer.template writeTensorField<1>(pop.flux(), layout);
        }
    }
}



template<typename H5Writer>
void FluidDiagnosticWriter<H5Writer>::MhdFluidWriter::operator()(auto const& layout)
{
    throw std::runtime_error("not implemented");
}



template<typename H5Writer>
void FluidDiagnosticWriter<H5Writer>::write(DiagnosticProperties& diagnostic)
{
    PHARE_LOG_SCOPE(3, "FluidDiagnosticWriter<H5Writer>::write");

    using Model_t    = H5Writer::ModelView::Model_t;
    auto& modelView  = this->h5Writer_.modelView();
    auto const& info = mem[diagnostic.quantity];

    modelView.onLevels(
        [&](auto const& level) {
            auto const ilvl = level.getLevelNumber();

            VTKFileWriter writer{diagnostic, this, info.offset_per_level[ilvl]};

            auto const write_quantity = [&](auto& layout, auto const&, auto const) {
                if constexpr (solver::is_hybrid_model_v<Model_t>)
                    HybridFluidWriter{this, diagnostic, writer}(layout);

                if constexpr (solver::is_mhd_model_v<Model_t>)
                    MhdFluidWriter{this, diagnostic, writer}(layout);
            };

            modelView.visitHierarchy(write_quantity, ilvl, ilvl);
        },
        this->h5Writer_.minLevel, this->h5Writer_.maxLevel);
}


} // namespace PHARE::diagnostic::vtkh5



#endif /* PHARE_DIAGNOSTIC_DETAIL_VTK_TYPES_FLUID_H */
