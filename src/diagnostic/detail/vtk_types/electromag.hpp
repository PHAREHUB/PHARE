#ifndef PHARE_DIAGNOSTIC_DETAIL_VTK_TYPES_ELECTROMAG_HPP
#define PHARE_DIAGNOSTIC_DETAIL_VTK_TYPES_ELECTROMAG_HPP

#include "diagnostic/detail/vtkh5_type_writer.hpp"

#include <string>
#include <vector>
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
        std::vector<std::size_t> offset_per_level = std::vector<std::size_t>(amr::MAX_LEVEL_IDX);
    };

    std::unordered_map<std::string, Info> mem;

    auto isActiveDiag(DiagnosticProperties const& diagnostic, std::string const& tree,
                      std::string var)
    {
        return diagnostic.quantity == tree + var;
    };
};


template<typename H5Writer>
void ElectromagDiagnosticWriter<H5Writer>::setup(DiagnosticProperties& diagnostic)
{
    VTKFileInitializer initializer{diagnostic, this};

    if (mem.count(diagnostic.quantity) == 0)
        mem.try_emplace(diagnostic.quantity);
    auto& info = mem[diagnostic.quantity];

    // EM_B/EM_E exist on every level regardless of which model owns it, so no model
    // resolution is needed here at all - init only touches file/box bookkeeping
    auto const init = [&](auto const& level) -> std::optional<std::size_t> {
        if (isActiveDiag(diagnostic, "/", "EM_B"))
            return initializer.template initTensorFieldFileLevel<1>(level);

        if (isActiveDiag(diagnostic, "/", "EM_E"))
            return initializer.template initTensorFieldFileLevel<1>(level);

        return std::nullopt;
    };

    this->h5Writer_.mapper().onLevels(
        [&](auto const& level) {
            auto const ilvl = level.getLevelNumber();
            if (auto const offset = init(ilvl))
                info.offset_per_level[ilvl] = *offset;
        },
        [&](int const ilvl) { // missing level
            init(ilvl);
        },
        this->h5Writer_.minLevel, this->h5Writer_.maxLevel);
}



template<typename H5Writer>
void ElectromagDiagnosticWriter<H5Writer>::write(DiagnosticProperties& diagnostic)
{
    auto& mapper = this->h5Writer_.mapper();
    auto& info   = mem[diagnostic.quantity];

    mapper.onLevels(
        [&](auto const& level) {
            auto const ilvl = level.getLevelNumber();

            VTKFileWriter fileWriter{diagnostic, this, info.offset_per_level[ilvl]};

            // EM_B/EM_E both exist on hybrid and mhd - resolve whichever model owns this
            // level and dispatch through it, rather than assuming a single fixed ModelView
            std::visit(
                [&](auto& modelView) {
                    auto const write_quantity = [&](auto& layout, auto const&, auto const) {
                        PHARE_LOG_SCOPE(3, "ElectromagDiagnosticWriter<H5Writer>::write_quantity");

                        if (isActiveDiag(diagnostic, "/", "EM_B"))
                            fileWriter.write(modelView.getB(), layout);
                        if (isActiveDiag(diagnostic, "/", "EM_E"))
                            fileWriter.write(modelView.getE(), layout);
                    };

                    modelView.visitHierarchy(write_quantity, ilvl, ilvl);
                },
                mapper.levelModelView(ilvl));
        },
        this->h5Writer_.minLevel, this->h5Writer_.maxLevel);
}


} // namespace PHARE::diagnostic::vtkh5

#endif /* PHARE_DIAGNOSTIC_DETAIL_VTK_TYPES_ELECTROMAG_HPP */
