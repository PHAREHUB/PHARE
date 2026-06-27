#ifndef PHARE_DIAGNOSTIC_DETAIL_VTK_TYPES_ELECTROMAG_HPP
#define PHARE_DIAGNOSTIC_DETAIL_VTK_TYPES_ELECTROMAG_HPP

#include "amr/utilities/box/amr_box.hpp"

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
    static constexpr auto dimension = H5Writer::dimension;

    struct Info
    {
        std::vector<std::size_t> offset_per_level
            = std::vector<std::size_t>(amr::MAX_LEVEL_IDX + 1);
        std::optional<core::Box<int, dimension>> slice_box;
    };

    void setup_normal(DiagnosticProperties&, Info&);
    void setup_sliced(DiagnosticProperties&, Info&, core::Box<int, H5Writer::dimension> const&);

    static auto get_slice_box(DiagnosticProperties const&)
        -> std::optional<core::Box<int, dimension>>;

    std::unordered_map<std::string, Info> mem;

    auto isActiveDiag(DiagnosticProperties const& diagnostic, std::string const& tree,
                      std::string var)
    {
        return diagnostic.quantity == tree + var;
    };
};


template<typename H5Writer>
auto ElectromagDiagnosticWriter<H5Writer>::get_slice_box(DiagnosticProperties const& diagnostic)
    -> std::optional<core::Box<int, dimension>>
{
    for (auto const& [key, val_node] : diagnostic.fileAttributes)
        if (key.find("_slice_") != std::string::npos)
            return core::Box<int, dimension>::fromString(val_node->template to<std::string>());
    return std::nullopt;
}


template<typename H5Writer>
void ElectromagDiagnosticWriter<H5Writer>::setup(DiagnosticProperties& diagnostic)
{
    if (mem.count(diagnostic.fileKey()) == 0)
        mem.try_emplace(diagnostic.fileKey());
    auto& info     = mem[diagnostic.fileKey()];
    info.slice_box = get_slice_box(diagnostic);

    if (info.slice_box)
        setup_sliced(diagnostic, info, *info.slice_box);
    else
        setup_normal(diagnostic, info);
}


template<typename H5Writer>
void ElectromagDiagnosticWriter<H5Writer>::setup_normal(DiagnosticProperties& diagnostic,
                                                        Info& info)
{
    auto& modelView = this->h5Writer_.modelView();
    VTKFileInitializer initializer{diagnostic, this};

    auto const init = [&](auto const ilvl) -> std::optional<std::size_t> {
        if (isActiveDiag(diagnostic, "/", "EM_B") || isActiveDiag(diagnostic, "/", "EM_E"))
            return initializer.template initTensorFieldFileLevel<1>(ilvl);
        return std::nullopt;
    };

    modelView.onLevels(
        [&](auto const& level) {
            auto const ilvl = level.getLevelNumber();
            if (auto const offset = init(ilvl))
                info.offset_per_level[ilvl] = *offset;
        },
        [&](int const ilvl) { init(ilvl); }, this->h5Writer_.minLevel, this->h5Writer_.maxLevel);
}


template<typename H5Writer>
void ElectromagDiagnosticWriter<H5Writer>::setup_sliced(
    DiagnosticProperties& diagnostic, Info& info,
    core::Box<int, H5Writer::dimension> const& slice_box)
{
    auto& modelView = this->h5Writer_.modelView();
    VTKFileInitializer initializer{diagnostic, this};

    auto const init = [&](auto const ilvl) -> std::optional<std::size_t> {
        if (isActiveDiag(diagnostic, "/", "EM_B") || isActiveDiag(diagnostic, "/", "EM_E"))
            return initializer.template initTensorFieldFileLevelWithSlice<1>(
                ilvl, amr::refine_box(slice_box, ilvl));
        return std::nullopt;
    };

    modelView.onLevels(
        [&](auto const& level) {
            auto const ilvl = level.getLevelNumber();
            if (auto const offset = init(ilvl))
                info.offset_per_level[ilvl] = *offset;
        },
        [&](int const ilvl) { init(ilvl); }, this->h5Writer_.minLevel, this->h5Writer_.maxLevel);
}


template<typename H5Writer>
void ElectromagDiagnosticWriter<H5Writer>::write(DiagnosticProperties& diagnostic)
{
    auto& modelView = this->h5Writer_.modelView();
    auto& info      = mem[diagnostic.fileKey()];

    modelView.onLevels(
        [&](auto const& level) {
            auto const ilvl = level.getLevelNumber();

            VTKFileWriter writer{diagnostic, this, info.offset_per_level[ilvl]};

            auto const write_quantity = [&](auto& layout, auto const&, auto const) {
                PHARE_LOG_SCOPE(3, "FluidDiagnosticWriter<H5Writer>::write_quantity");

                if (isActiveDiag(diagnostic, "/", "EM_B"))
                {
                    auto& B = this->h5Writer_.modelView().getB();
                    if (info.slice_box)
                        writer.template writeTensorFieldSlice<1>(
                            B, layout, amr::refine_box(*info.slice_box, ilvl));
                    else
                        writer.template writeTensorField<1>(B, layout);
                }
                if (isActiveDiag(diagnostic, "/", "EM_E"))
                {
                    auto& E = this->h5Writer_.modelView().getE();
                    if (info.slice_box)
                        writer.template writeTensorFieldSlice<1>(
                            E, layout, amr::refine_box(*info.slice_box, ilvl));
                    else
                        writer.template writeTensorField<1>(E, layout);
                }
            };

            modelView.visitHierarchy(write_quantity, ilvl, ilvl);
        },
        this->h5Writer_.minLevel, this->h5Writer_.maxLevel);
}


} // namespace PHARE::diagnostic::vtkh5

#endif /* PHARE_DIAGNOSTIC_DETAIL_VTK_TYPES_ELECTROMAG_HPP */
