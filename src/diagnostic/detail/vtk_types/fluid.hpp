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
    using ModelView          = H5Writer::ModelView;
    using Model_t            = H5Writer::ModelView::Model_t;
    using Field_t            = Model_t::field_type;

public:
    FluidDiagnosticWriter(H5Writer& h5Writer)
        : Super{h5Writer}
    {
    }

    void setup(DiagnosticProperties&) override;
    void write(DiagnosticProperties&) override;
    void compute_as_needed(DiagnosticProperties&) override;

private:
    struct Info
    {
        std::vector<std::size_t> offset_per_level
            = std::vector<std::size_t>(amr::MAX_LEVEL_IDX + 1);
    };

    // action(field, name, ownerKey) on the quantity of the diagnostic, if this model has it
    void dispatch_(DiagnosticProperties const& diagnostic, auto&& action)
        requires(solver::is_hybrid_model_v<Model_t>);
    void dispatch_(DiagnosticProperties const& diagnostic, auto&& action)
        requires(solver::is_mhd_model_v<Model_t>);

    void compute_as_needed_(DiagnosticProperties& diagnostic)
        requires(solver::is_hybrid_model_v<Model_t>);
    void compute_as_needed_(DiagnosticProperties& diagnostic)
        requires(solver::is_mhd_model_v<Model_t>);

    std::unordered_map<std::string, Info> mem;
};



template<typename H5Writer>
void FluidDiagnosticWriter<H5Writer>::setup(DiagnosticProperties& diagnostic)
{
    PHARE_LOG_SCOPE(3, "FluidDiagnosticWriter<H5Writer>::setup");

    auto& modelView = this->h5Writer_.modelView();

    VTKFileInitializer initializer{diagnostic, this};

    if (mem.count(diagnostic.quantity) == 0)
        mem.try_emplace(diagnostic.quantity);
    auto& info = mem[diagnostic.quantity];

    auto const init = [&](auto const ilvl) -> std::optional<std::size_t> {
        std::optional<std::size_t> offset;
        dispatch_(diagnostic, [&](auto& field, auto const&, auto const&) {
            using Quantity = std::decay_t<decltype(field)>;
            if constexpr (std::is_same_v<Quantity, Field_t>)
                offset = initializer.initFieldFileLevel(ilvl);
            else
                offset = initializer.template initTensorFieldFileLevel<Quantity::rank>(ilvl);
        });
        return offset;
    };

    modelView.onLevels(
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

    if constexpr (solver::is_hybrid_model_v<Model_t>)
        for (auto const& pop : modelView.getIons())
            if (diagnostic.quantity.starts_with("/ions/pop/" + pop.name() + "/"))
            {
                typename H5Writer::Attributes popAttributes;
                popAttributes["pop_mass"]   = pop.mass();
                popAttributes["pop_charge"] = pop.charge();
                this->h5Writer_.writeGlobalAttributeDict(this->getOrCreateH5File(diagnostic),
                                                         popAttributes, "/");
            }
}



template<typename H5Writer>
void FluidDiagnosticWriter<H5Writer>::write(DiagnosticProperties& diagnostic)
{
    PHARE_LOG_SCOPE(3, "FluidDiagnosticWriter<H5Writer>::write");

    auto& modelView  = this->h5Writer_.modelView();
    auto const& info = mem[diagnostic.quantity];

    modelView.onLevels(
        [&](auto const& level) {
            auto const ilvl = level.getLevelNumber();

            VTKFileWriter writer{diagnostic, this, info.offset_per_level[ilvl]};

            auto const write_quantity = [&](auto& layout, auto const&, auto const) {
                dispatch_(diagnostic, [&](auto& field, auto const&, auto const&) {
                    using Quantity = std::decay_t<decltype(field)>;
                    if constexpr (std::is_same_v<Quantity, Field_t>)
                        writer.writeField(field, layout);
                    else
                        writer.template writeTensorField<Quantity::rank>(field, layout);
                });
            };

            modelView.visitHierarchy(write_quantity, ilvl, ilvl);
        },
        this->h5Writer_.minLevel, this->h5Writer_.maxLevel);
}



template<typename H5Writer>
void FluidDiagnosticWriter<H5Writer>::dispatch_(DiagnosticProperties const& diagnostic,
                                                auto&& action)
    requires(solver::is_hybrid_model_v<Model_t>)
{
    auto& modelView = this->h5Writer_.modelView();
    ModelView::FluidAccessors::getOrCreateFor(modelView).dispatch(diagnostic.quantity, modelView,
                                                                  action);
}

template<typename H5Writer>
void FluidDiagnosticWriter<H5Writer>::dispatch_(DiagnosticProperties const& diagnostic,
                                                auto&& action)
    requires(solver::is_mhd_model_v<Model_t>)
{
    auto& modelView = this->h5Writer_.modelView();
    ModelView::MHDAccessors::getOrCreateFor(modelView).dispatch(diagnostic.quantity, modelView,
                                                                action);
}



template<typename H5Writer>
void FluidDiagnosticWriter<H5Writer>::compute_as_needed(DiagnosticProperties& diagnostic)
{
    compute_as_needed_(diagnostic);
}

template<typename H5Writer>
void FluidDiagnosticWriter<H5Writer>::compute_as_needed_(DiagnosticProperties& diagnostic)
    requires(solver::is_hybrid_model_v<Model_t>)
{
    auto& h5Writer  = this->h5Writer_;
    auto& modelView = h5Writer.modelView();

    ModelView::FluidComputers::getOrCreateFor(modelView).compute(
        diagnostic.quantity, modelView, h5Writer.minLevel, h5Writer.maxLevel, h5Writer.timestamp());
}

template<typename H5Writer>
void FluidDiagnosticWriter<H5Writer>::compute_as_needed_(DiagnosticProperties& diagnostic)
    requires(solver::is_mhd_model_v<Model_t>)
{
    auto& h5Writer  = this->h5Writer_;
    auto& modelView = h5Writer.modelView();
    // or FloatType if we want to expose that to DiagnosticProperties
    double const gamma = diagnostic.fileAttributes["heat_capacity_ratio"].template to<double>();

    ModelView::MHDComputers::getOrCreateFor(modelView).compute(
        diagnostic.quantity, modelView, gamma, h5Writer.minLevel, h5Writer.maxLevel,
        h5Writer.timestamp());
}

} // namespace PHARE::diagnostic::vtkh5



#endif /* PHARE_DIAGNOSTIC_DETAIL_VTK_TYPES_FLUID_HPP */
