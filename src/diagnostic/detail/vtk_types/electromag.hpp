#ifndef PHARE_DIAGNOSTIC_DETAIL_VTK_TYPES_ELECTROMAG_HPP
#define PHARE_DIAGNOSTIC_DETAIL_VTK_TYPES_ELECTROMAG_HPP

#include "diagnostic/detail/vtkh5_type_writer.hpp"

#include <string>
#include <optional>
#include <unordered_map>

namespace PHARE::diagnostic::vtkh5
{

template<typename H5Writer>
class ElectromagDiagnosticWriter : public H5TypeWriter<H5Writer>
{
    using Super              = H5TypeWriter<H5Writer>;
    using VTKFileWriter      = Super::VTKFileWriter;
    using VTKFileInitializer = Super::VTKFileInitializer;

public:
    ElectromagDiagnosticWriter(H5Writer& h5Writer)
        : Super{h5Writer}
    {
    }

    void setup(DiagnosticProperties&) override;
    void write(DiagnosticProperties&) override;
    void compute(DiagnosticProperties&) override {}

private:
    struct Info
    {
        std::vector<std::size_t> offset_per_level = std::vector<std::size_t>(H5Writer::MAX_LEVEL);
    };

    std::unordered_map<std::string, Info> mem;
};


template<typename H5Writer>
void ElectromagDiagnosticWriter<H5Writer>::setup(DiagnosticProperties& diagnostic)
{
    auto& modelView = this->h5Writer_.modelView();
    VTKFileInitializer initializer{diagnostic, this};

    if (mem.count(diagnostic.quantity) == 0)
        mem.try_emplace(diagnostic.quantity);
    auto& info = mem[diagnostic.quantity];

    auto const init = [&](auto const& level) -> std::optional<std::size_t> {
        for (auto* vecField : this->h5Writer_.modelView().getElectromagFields())
            if (diagnostic.quantity == "/" + vecField->name())
                return initializer.template initTensorFieldFileLevel<1>(level);
        return std::nullopt;
    };

    modelView.onLevels(
        [&](auto const& level) {
            auto const ilvl = level.getLevelNumber();
            initializer.initFileLevel(ilvl);
            if (auto const offset = init(level))
                info.offset_per_level[ilvl] = *offset;
        },
        this->h5Writer_.minLevel, this->h5Writer_.maxLevel);
}

template<typename H5Writer>
void ElectromagDiagnosticWriter<H5Writer>::write(DiagnosticProperties& diagnostic)
{
    auto& modelView = this->h5Writer_.modelView();
    auto& info      = mem[diagnostic.quantity];

    auto const write_quantity = [&](auto& layout, auto const&, auto const iLevel) {
        VTKFileWriter writer{diagnostic, this, info.offset_per_level[iLevel]};

        for (auto* vecField : this->h5Writer_.modelView().getElectromagFields())
            if (diagnostic.quantity == "/" + vecField->name())
                writer.template writeTensorField<1>(*vecField, layout);
    };

    modelView.visitHierarchy(write_quantity, this->h5Writer_.minLevel, this->h5Writer_.maxLevel);
}




} // namespace PHARE::diagnostic::vtkh5

#endif /* PHARE_DIAGNOSTIC_DETAIL_VTK_TYPES_ELECTROMAG_HPP */
