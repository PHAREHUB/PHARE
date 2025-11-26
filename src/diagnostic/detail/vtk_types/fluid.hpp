#ifndef PHARE_DIAGNOSTIC_DETAIL_VTK_TYPES_FLUID_HPP
#define PHARE_DIAGNOSTIC_DETAIL_VTK_TYPES_FLUID_HPP

#include "core/logger.hpp"
#include "diagnostic/detail/vtkh5_type_writer.hpp"

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
        Info()
            : offset_per_level(10)
        {
        }

        std::vector<std::size_t> offset_per_level;
    };

    auto static isActiveDiag(DiagnosticProperties const& diagnostic, std::string const& tree,
                             std::string const& var)
    {
        return diagnostic.quantity == tree + var;
    };

    std::unordered_map<std::string, Info> mem;
};


template<typename H5Writer>
void FluidDiagnosticWriter<H5Writer>::setup(DiagnosticProperties& diagnostic)
{
    auto& modelView = this->h5Writer_.modelView();
    auto& ions      = modelView.getIons();
    VTKFileInitializer initializer{diagnostic, this};
    std::string const tree{"/ions/"};

    if (mem.count(diagnostic.quantity) == 0)
        mem.try_emplace(diagnostic.quantity);
    auto& info = mem[diagnostic.quantity];

    auto const init = [&](auto const& level) -> std::optional<std::size_t> {
        if (isActiveDiag(diagnostic, tree, "charge_density"))
            return initializer.initFieldFileLevel(level);
        else if (isActiveDiag(diagnostic, tree, "mass_density"))
            return initializer.initFieldFileLevel(level);
        else if (isActiveDiag(diagnostic, tree, "bulkVelocity"))
            return initializer.template initTensorFieldFileLevel<1>(level);
        else
        {
            for (auto& pop : ions)
            {
                auto const pop_tree = tree + "pop/" + pop.name() + "/";
                if (isActiveDiag(diagnostic, pop_tree, "density"))
                    return initializer.initFieldFileLevel(level);
                else if (isActiveDiag(diagnostic, pop_tree, "charge_density"))
                    return initializer.initFieldFileLevel(level);
                else if (isActiveDiag(diagnostic, pop_tree, "flux"))
                    return initializer.template initTensorFieldFileLevel<1>(level);
            }
        }
        return std::nullopt;
    };

    modelView.onLevels(
        [&](auto const& level) {
            PHARE_LOG_SCOPE(3, "FluidDiagnosticWriter<H5Writer>::write_level");

            auto const ilvl = level.getLevelNumber();
            initializer.initFileLevel(ilvl);
            if (auto const offset = init(level))
                info.offset_per_level[ilvl] = *offset;
        },
        this->h5Writer_.minLevel, this->h5Writer_.maxLevel);
}

template<typename H5Writer>
void FluidDiagnosticWriter<H5Writer>::write(DiagnosticProperties& diagnostic)
{
    PHARE_LOG_SCOPE(3, "FluidDiagnosticWriter<H5Writer>::write");

    auto& modelView = this->h5Writer_.modelView();
    auto& ions      = modelView.getIons();
    auto& info      = mem[diagnostic.quantity];
    std::string const tree{"/ions/"};

    auto const write_quantity = [&](auto& layout, auto const&, auto const iLevel) {
        PHARE_LOG_SCOPE(3, "FluidDiagnosticWriter<H5Writer>::write_quantity");

        VTKFileWriter writer{diagnostic, this, info.offset_per_level[iLevel]};

        if (isActiveDiag(diagnostic, tree, "charge_density"))
            writer.writeField(ions.chargeDensity(), layout);
        else if (isActiveDiag(diagnostic, tree, "mass_density"))
            writer.writeField(ions.massDensity(), layout);
        else if (isActiveDiag(diagnostic, tree, "bulkVelocity"))
            writer.template writeTensorField<1>(ions.velocity(), layout);
        else
        {
            for (auto& pop : ions)
            {
                auto const pop_tree = tree + "pop/" + pop.name() + "/";
                if (isActiveDiag(diagnostic, pop_tree, "density"))
                    writer.writeField(pop.particleDensity(), layout);
                else if (isActiveDiag(diagnostic, pop_tree, "charge_density"))
                    writer.writeField(pop.chargeDensity(), layout);
                else if (isActiveDiag(diagnostic, pop_tree, "flux"))
                    writer.template writeTensorField<1>(pop.flux(), layout);
            }
        }
    };

    modelView.visitHierarchy(write_quantity, this->h5Writer_.minLevel, this->h5Writer_.maxLevel);
}


} // namespace PHARE::diagnostic::vtkh5



#endif /* PHARE_DIAGNOSTIC_DETAIL_VTK_TYPES_FLUID_H */
