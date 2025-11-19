#ifndef PHARE_DIAGNOSTIC_DETAIL_VTK_TYPES_FLUID_HPP
#define PHARE_DIAGNOSTIC_DETAIL_VTK_TYPES_FLUID_HPP


#include "diagnostic/detail/vtkh5_type_writer.hpp"

#include <string>

namespace PHARE::diagnostic::vtkh5
{


template<typename H5Writer>
class FluidDiagnosticWriter : public H5TypeWriter<H5Writer>
{
    using Super         = H5TypeWriter<H5Writer>;
    using VTKFileWriter = Super::VTKFileWriter;

public:
    FluidDiagnosticWriter(H5Writer& h5Writer)
        : Super{h5Writer}
    {
    }

    void write(DiagnosticProperties&) override;
    void compute(DiagnosticProperties&) override {};

private:
    auto static isActiveDiag(DiagnosticProperties const& diagnostic, std::string const& tree,
                             std::string const& var)
    {
        return diagnostic.quantity == tree + var;
    };
};


template<typename H5Writer>
void FluidDiagnosticWriter<H5Writer>::write(DiagnosticProperties& diagnostic)
{
    auto& modelView = this->h5Writer_.modelView();
    auto& ions      = modelView.getIons();
    VTKFileWriter writer{diagnostic, this};

    std::string const tree{"/ions/"};

    auto const write_quantity = [&](auto& layout, auto const&, auto const iLevel) {
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

    modelView.onLevels(
        [&](auto const& lvl) {
            auto const ilvl = lvl.getLevelNumber();
            writer.initFileLevel(ilvl);

            auto boxes = modelView.localLevelBoxes(ilvl);
            if (isActiveDiag(diagnostic, tree, "charge_density"))
                writer.initFieldFileLevel(ilvl, boxes);
            else if (isActiveDiag(diagnostic, tree, "mass_density"))
                writer.initFieldFileLevel(ilvl, boxes);
            else if (isActiveDiag(diagnostic, tree, "bulkVelocity"))
                writer.template initTensorFieldFileLevel<1>(ilvl, boxes);
            else
            {
                for (auto& pop : ions)
                {
                    auto const pop_tree = tree + "pop/" + pop.name() + "/";
                    if (isActiveDiag(diagnostic, pop_tree, "density"))
                        writer.initFieldFileLevel(ilvl, boxes);
                    else if (isActiveDiag(diagnostic, pop_tree, "charge_density"))
                        writer.initFieldFileLevel(ilvl, boxes);
                    else if (isActiveDiag(diagnostic, pop_tree, "flux"))
                        writer.template initTensorFieldFileLevel<1>(ilvl, boxes);
                }
            }

            modelView.visitHierarchy(write_quantity, ilvl, ilvl);
        },
        this->h5Writer_.minLevel, this->h5Writer_.maxLevel);
}


} // namespace PHARE::diagnostic::vtkh5



#endif /* PHARE_DIAGNOSTIC_DETAIL_VTK_TYPES_FLUID_H */
