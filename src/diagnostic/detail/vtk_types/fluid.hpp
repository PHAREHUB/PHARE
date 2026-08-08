#ifndef PHARE_DIAGNOSTIC_DETAIL_VTK_TYPES_FLUID_HPP
#define PHARE_DIAGNOSTIC_DETAIL_VTK_TYPES_FLUID_HPP

#include "core/logger.hpp"

#include "amr/physical_models/mhd_model.hpp"
#include "amr/physical_models/hybrid_model.hpp"
#include "diagnostic/detail/vtkh5_type_writer.hpp"

#include <string>
#include <vector>
#include <optional>
#include <stdexcept>
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
    void compute(DiagnosticProperties&) override;

    // resolves the ModelView for diagnostic.type ("fluid"->hybrid, "mhd"->mhd) by checking
    // is_hybrid_model_v/is_mhd_model_v on the alternatives actually present in modelViews -
    // unlike mapper().hyridModelView()/mhdModelView(), never names ModelView<Hierarchy, void>,
    // so it stays safe to call even when the other model isn't in this build's Models... pack
    auto& ownModelView(DiagnosticProperties& diagnostic)
    {
        for (auto& mv : this->h5Writer_.mapper().modelViews)
        {
            bool const isMine = std::visit(
                [&](auto& modelView) {
                    using Model_t = typename std::decay_t<decltype(modelView)>::Model_t;
                    if (diagnostic.type == "fluid")
                        return solver::is_hybrid_model_v<Model_t>;
                    if (diagnostic.type == "mhd")
                        return solver::is_mhd_model_v<Model_t>;
                    return false;
                },
                mv);
            if (isMine)
                return mv;
        }
        throw std::runtime_error("FluidDiagnosticWriter: no model registered for diagnostic type "
                                 + diagnostic.type);
    }

private:
    struct Info
    {
        std::vector<std::size_t> offset_per_level
            = std::vector<std::size_t>(amr::MAX_LEVEL_IDX + 1);
    };

    auto static isActiveDiag(DiagnosticProperties const& diagnostic, std::string const& tree,
                             std::string const& var)
    {
        return diagnostic.quantity == tree + var;
    };

    // returns nullopt for the model that doesn't match Model_t - the caller (setup's init
    // lambda) already resolved to the right alternative via ownModelView, so in practice only
    // one branch is ever meaningfully exercised, but std::visit still instantiates this body
    // once per alternative in ModelViewVariant, hence the if constexpr guards
    template<typename ModelView>
    std::optional<std::size_t>
    initHybridFluidLevel(ModelView& modelView, DiagnosticProperties& diagnostic,
                         VTKFileInitializer& file_initializer, auto const ilvl)
    {
        using Model_t = ModelView::Model_t;
        if constexpr (solver::is_hybrid_model_v<Model_t>)
        {
            using Accessors = ModelView::FluidAccessors;

            auto& accessors = Accessors::getOrCreateFor(modelView);
            auto& model     = modelView.model();

            std::optional<std::size_t> offset;
            accessors.dispatch(
                diagnostic.quantity, model,
                [&](auto& field, std::string const&, std::string const&) {
                    using Quantity = std::decay_t<decltype(field)>;
                    if constexpr (std::is_same_v<Quantity, typename Accessors::Field_t>)
                        offset = file_initializer.initFieldFileLevel(ilvl);
                    else
                        offset = file_initializer.template initTensorFieldFileLevel<Quantity::rank>(
                            ilvl);
                });
            return offset;
        }
        return std::nullopt;
    }

    template<typename ModelView>
    std::optional<std::size_t>
    initMhdFluidLevel(ModelView& modelView, DiagnosticProperties& diagnostic,
                      VTKFileInitializer& file_initializer, auto const ilvl)
    {
        using Model_t = ModelView::Model_t;
        if constexpr (solver::is_mhd_model_v<Model_t>)
        {
            using Accessors = ModelView::MHDAccessors;

            auto& accessors = Accessors::getOrCreateFor(modelView);
            auto& model     = modelView.model();

            std::optional<std::size_t> offset;
            accessors.dispatch(
                diagnostic.quantity, model,
                [&](auto& field, std::string const&, std::string const&) {
                    using Quantity = std::decay_t<decltype(field)>;
                    if constexpr (std::is_same_v<Quantity, typename Accessors::Field_t>)
                        offset = file_initializer.initFieldFileLevel(ilvl);
                    else
                        offset = file_initializer.template initTensorFieldFileLevel<1>(ilvl);
                });
            return offset;
        }
        return std::nullopt;
    }

    template<typename ModelView>
    void writeHybridFluid(ModelView& modelView, DiagnosticProperties& diagnostic,
                          VTKFileWriter& file_writer, auto const& layout)
    {
        using Model_t = ModelView::Model_t;
        if constexpr (solver::is_hybrid_model_v<Model_t>)
        {
            using Accessors = ModelView::FluidAccessors;

            auto& accessors = Accessors::getOrCreateFor(modelView);
            auto& model     = modelView.model();

            accessors.dispatch(diagnostic.quantity, model,
                               [&](auto& quantity, std::string const&, std::string const&) {
                                   file_writer.write(quantity, layout);
                               });
        }
    }

    template<typename ModelView>
    void writeMhdFluid(ModelView& modelView, DiagnosticProperties& diagnostic,
                       VTKFileWriter& file_writer, auto const& layout)
    {
        using Model_t = ModelView::Model_t;
        if constexpr (solver::is_mhd_model_v<Model_t>)
        {
            using Accessors = ModelView::MHDAccessors;

            auto& accessors = Accessors::getOrCreateFor(modelView);
            auto& model     = modelView.model();

            accessors.dispatch(diagnostic.quantity, model,
                               [&](auto& quantity, std::string const&, std::string const&) {
                                   file_writer.write(quantity, layout);
                               });
        }
    }

    template<typename ModelView>
    void computeMhdFluid(ModelView& modelView, DiagnosticProperties& diagnostic)
    {
        using Model_t = ModelView::Model_t;
        if constexpr (solver::is_mhd_model_v<Model_t>)
        {
            using Computers = ModelView::MHDComputers;

            auto& computers    = Computers::getOrCreateFor(modelView);
            double const gamma = diagnostic.fileAttributes["heat_capacity_ratio"];

            computers.compute(diagnostic.quantity, modelView, gamma, this->h5Writer_.minLevel,
                              this->h5Writer_.maxLevel, this->h5Writer_.timestamp());
        }
    }

    std::unordered_map<std::string, Info> mem;
};



template<typename H5Writer>
void FluidDiagnosticWriter<H5Writer>::setup(DiagnosticProperties& diagnostic)
{
    PHARE_LOG_SCOPE(3, "FluidDiagnosticWriter<H5Writer>::setup");

    VTKFileInitializer initializer{diagnostic, this};

    if (mem.count(diagnostic.quantity) == 0)
        mem.try_emplace(diagnostic.quantity);
    auto& info = mem[diagnostic.quantity];

    // this writer instance is registered once for "fluid" (hybrid) diagnostics and once for
    // "mhd" ones - each instance only ever sees diagnostics of its own type. population names
    // (used below to match per-pop quantities) are global to the simulation, not per-level,
    // so resolving the model once here (rather than per-level) is correct even for levels
    // that don't exist yet.
    auto& mv = ownModelView(diagnostic);

    auto const init = [&](auto const ilvl) -> std::optional<std::size_t> {
        return std::visit(
            [&](auto& modelView) -> std::optional<std::size_t> {
                if (auto ret = initHybridFluidLevel(modelView, diagnostic, initializer, ilvl))
                    return ret;
                if (auto ret = initMhdFluidLevel(modelView, diagnostic, initializer, ilvl))
                    return ret;
                return std::nullopt;
            },
            mv);
    };

    this->h5Writer_.mapper().onLevels(
        [&](auto const& level) {
            PHARE_LOG_SCOPE(3, "FluidDiagnosticWriter<H5Writer>::setup_level");

            auto const ilvl = level.getLevelNumber();
            if (auto const offset = init(ilvl))
                info.offset_per_level[ilvl] = *offset;
        },
        [&](int const ilvl) {
            PHARE_LOG_SCOPE(3, "FluidDiagnosticWriter<H5Writer>::setup_missing_level");

            init(ilvl);
        },
        this->h5Writer_.minLevel, this->h5Writer_.maxLevel);
}



template<typename H5Writer>
void FluidDiagnosticWriter<H5Writer>::write(DiagnosticProperties& diagnostic)
{
    PHARE_LOG_SCOPE(3, "FluidDiagnosticWriter<H5Writer>::write");

    auto& mapper     = this->h5Writer_.mapper();
    auto const& info = mem[diagnostic.quantity];
    auto& mv         = ownModelView(diagnostic);

    mapper.onLevels(
        [&](auto const& level) {
            auto const ilvl = level.getLevelNumber();

            VTKFileWriter fileWriter{diagnostic, this, info.offset_per_level[ilvl]};

            std::visit(
                [&](auto& modelView) {
                    auto const write_quantity = [&](auto& layout, auto const&, auto const) {
                        if (diagnostic.type == "fluid")
                            writeHybridFluid(modelView, diagnostic, fileWriter, layout);
                        else if (diagnostic.type == "mhd")
                            writeMhdFluid(modelView, diagnostic, fileWriter, layout);
                    };

                    modelView.visitHierarchy(write_quantity, ilvl, ilvl);
                },
                mv);
        },
        this->h5Writer_.minLevel, this->h5Writer_.maxLevel);
}



template<typename H5Writer>
void FluidDiagnosticWriter<H5Writer>::compute(DiagnosticProperties& diagnostic)
{
    auto& mv = ownModelView(diagnostic);
    std::visit([&](auto& modelView) { computeMhdFluid(modelView, diagnostic); }, mv);
}

} // namespace PHARE::diagnostic::vtkh5



#endif /* PHARE_DIAGNOSTIC_DETAIL_VTK_TYPES_FLUID_HPP */
