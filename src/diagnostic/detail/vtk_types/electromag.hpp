#ifndef PHARE_DIAGNOSTIC_DETAIL_VTK_TYPES_ELECTROMAG_HPP
#define PHARE_DIAGNOSTIC_DETAIL_VTK_TYPES_ELECTROMAG_HPP

#include "diagnostic/detail/vtkh5_type_writer.hpp"



namespace PHARE::diagnostic::vtkh5
{


template<typename H5Writer>
class ElectromagDiagnosticWriter : public H5TypeWriter<H5Writer>
{
    using Super         = H5TypeWriter<H5Writer>;
    using VTKFileWriter = Super::VTKFileWriter;

public:
    ElectromagDiagnosticWriter(H5Writer& h5Writer)
        : Super{h5Writer}
    {
    }
    void write(DiagnosticProperties&) override;

    void compute(DiagnosticProperties&) override {}
};


template<typename H5Writer>
void ElectromagDiagnosticWriter<H5Writer>::write(DiagnosticProperties& diagnostic)
{
    auto& modelView = this->h5Writer_.modelView();

    VTKFileWriter writer{diagnostic, this};

    auto const write_quantity = [&](auto& layout, auto const&, auto const iLevel) {
        for (auto* vecField : this->h5Writer_.modelView().getElectromagFields())
            if (diagnostic.quantity == "/" + vecField->name())
                writer.template writeTensorField<1>(*vecField, layout);
    };

    modelView.onLevels(
        [&](auto const& lvl) {
            auto const ilvl = lvl.getLevelNumber();
            writer.initFileLevel(ilvl);

            auto boxes = modelView.localLevelBoxes(ilvl);
            for (auto* vecField : this->h5Writer_.modelView().getElectromagFields())
                if (diagnostic.quantity == "/" + vecField->name())
                    writer.template initTensorFieldFileLevel<1>(ilvl, boxes);

            modelView.visitHierarchy(write_quantity, ilvl, ilvl);
        },
        this->h5Writer_.minLevel, this->h5Writer_.maxLevel);
}




} // namespace PHARE::diagnostic::vtkh5

#endif /* PHARE_DIAGNOSTIC_DETAIL_VTK_TYPES_ELECTROMAG_HPP */
