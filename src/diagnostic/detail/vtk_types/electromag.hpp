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
    auto& modelView = this->h5Writer_.modelView();
    VTKFileInitializer initializer{diagnostic, this};

    if (mem.count(diagnostic.quantity) == 0)
        mem.try_emplace(diagnostic.quantity);
    auto& info = mem[diagnostic.quantity];

    // assumes exists for all models
    auto const init = [&](auto const& level) -> std::optional<std::size_t> {
        if (isActiveDiag(diagnostic, "/", "EM_B"))
        {
            return initializer.template initTensorFieldFileLevel<1>(level);
        }
        if (isActiveDiag(diagnostic, "/", "EM_E"))
        {
            return initializer.template initTensorFieldFileLevel<1>(level);
        }

        return std::nullopt;
    };

    modelView.onLevels(
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
    auto& modelView = this->h5Writer_.modelView();
    auto& info      = mem[diagnostic.quantity];

    modelView.onLevels(
        [&](auto const& level) {
            auto const ilvl = level.getLevelNumber();

            VTKFileWriter writer{diagnostic, this, info.offset_per_level[ilvl]};

            auto const write_quantity = [&](auto& layout, auto const&, auto const) {
                PHARE_LOG_SCOPE(3, "FluidDiagnosticWriter<H5Writer>::write_quantity");

                if (isActiveDiag(diagnostic, "/", "EM_B"))
                {
                    auto& B = this->h5Writer_.modelView().getB();
                    writer.template writeTensorField<1>(B, layout);
                }
                if (isActiveDiag(diagnostic, "/", "EM_E"))
                {
                    auto& E = this->h5Writer_.modelView().getE();
                    writer.template writeTensorField<1>(E, layout);
                }
            };

            modelView.visitHierarchy(write_quantity, ilvl, ilvl);
        },
        this->h5Writer_.minLevel, this->h5Writer_.maxLevel);
}


} // namespace PHARE::diagnostic::vtkh5

#endif /* PHARE_DIAGNOSTIC_DETAIL_VTK_TYPES_ELECTROMAG_HPP */
