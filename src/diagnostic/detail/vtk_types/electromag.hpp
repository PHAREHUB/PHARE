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
    using GridLayout         = H5Writer::GridLayout;

public:
    ElectromagDiagnosticWriter(H5Writer& h5Writer)
        : Super{h5Writer}
    {
    }

    void setup(DiagnosticProperties&) override;
    void write(DiagnosticProperties&) override;
    void compute(DiagnosticProperties&) override;

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
void ElectromagDiagnosticWriter<H5Writer>::compute(DiagnosticProperties& diagnostic)
{
    auto& modelView = this->h5Writer_.modelView();
    auto minLvl     = this->h5Writer_.minLevel;
    auto maxLvl     = this->h5Writer_.maxLevel;

    if constexpr (requires {
                      modelView.getB1();
                      modelView.getB0();
                  })
    {
        if (isActiveDiag(diagnostic, "/", "EM_B"))
        {
            auto& B        = modelView.getB();
            auto const& B1 = modelView.getB1();
            auto const& B0 = modelView.getB0();

            auto reconstructB = [&](GridLayout& layout, std::string const&, std::size_t) {
                auto& Bx        = B.getComponent(core::Component::X);
                auto& By        = B.getComponent(core::Component::Y);
                auto& Bz        = B.getComponent(core::Component::Z);
                auto const& B1x = B1.getComponent(core::Component::X);
                auto const& B1y = B1.getComponent(core::Component::Y);
                auto const& B1z = B1.getComponent(core::Component::Z);
                auto const& B0x = B0.getComponent(core::Component::X);
                auto const& B0y = B0.getComponent(core::Component::Y);
                auto const& B0z = B0.getComponent(core::Component::Z);

                auto const rebuildComponent
                    = [&](auto& dst, auto const& perturbed, auto const& background) {
                          layout.evalOnGhostBox(dst, [&](auto&... args) mutable {
                              dst(args...) = perturbed(args...) + background(args...);
                          });
                      };

                rebuildComponent(Bx, B1x, B0x);
                rebuildComponent(By, B1y, B0y);
                rebuildComponent(Bz, B1z, B0z);
            };

            modelView.visitHierarchy(reconstructB, minLvl, maxLvl);
            return;
        }
    }

    if constexpr (requires {
                      modelView.getDivB();
                      modelView.getB1();
                      modelView.getB0();
                  })
    {
        if (isActiveDiag(diagnostic, "/", "EM_divB"))
        {
            // Cell-centered divergence of the total field B = B0 + B1. deriv is linear, so sum the
            // B0 and B1 contributions per direction; each deriv<dir> on a face-centered component
            // lands on the cell centre where divB lives.
            auto& divB     = modelView.getDivB();
            auto const& B1 = modelView.getB1();
            auto const& B0 = modelView.getB0();

            modelView.visitHierarchy(
                [&](GridLayout& layout, std::string const&, std::size_t) {
                    auto const& B1x = B1.getComponent(core::Component::X);
                    auto const& B0x = B0.getComponent(core::Component::X);
                    // divB is pure dual (cell-centered); B1 and B0 are defined on the full ghost
                    // box, so divB is evaluated there too (not just the interior) — the dump reads
                    // divB ghost cells, which would otherwise be uninitialized (NaN).
                    layout.evalOnGhostBox(divB, [&](auto&... args) mutable {
                        auto d = layout.template deriv<core::Direction::X>(B1x, {args...})
                                 + layout.template deriv<core::Direction::X>(B0x, {args...});
                        if constexpr (GridLayout::dimension >= 2)
                        {
                            auto const& B1y = B1.getComponent(core::Component::Y);
                            auto const& B0y = B0.getComponent(core::Component::Y);
                            d += layout.template deriv<core::Direction::Y>(B1y, {args...})
                                 + layout.template deriv<core::Direction::Y>(B0y, {args...});
                        }
                        if constexpr (GridLayout::dimension == 3)
                        {
                            auto const& B1z = B1.getComponent(core::Component::Z);
                            auto const& B0z = B0.getComponent(core::Component::Z);
                            d += layout.template deriv<core::Direction::Z>(B1z, {args...})
                                 + layout.template deriv<core::Direction::Z>(B0z, {args...});
                        }
                        divB(args...) = d;
                    });
                },
                minLvl, maxLvl);
            return;
        }
    }
}


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
        if constexpr (requires { modelView.getE(); })
            if (isActiveDiag(diagnostic, "/", "EM_E"))
            {
                return initializer.template initTensorFieldFileLevel<1>(level);
            }
        if constexpr (requires { modelView.getB1(); })
            if (isActiveDiag(diagnostic, "/", "EM_B1"))
            {
                return initializer.template initTensorFieldFileLevel<1>(level);
            }
        if constexpr (requires { modelView.getDivB(); })
            if (isActiveDiag(diagnostic, "/", "EM_divB"))
            {
                return initializer.initFieldFileLevel(level);
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
                if constexpr (requires { modelView.getE(); })
                    if (isActiveDiag(diagnostic, "/", "EM_E"))
                    {
                        auto& E = this->h5Writer_.modelView().getE();
                        writer.template writeTensorField<1>(E, layout);
                    }
                if constexpr (requires { modelView.getB1(); })
                    if (isActiveDiag(diagnostic, "/", "EM_B1"))
                    {
                        auto const& B1 = modelView.getB1();
                        writer.template writeTensorField<1>(B1, layout);
                    }
                if constexpr (requires { modelView.getDivB(); })
                    if (isActiveDiag(diagnostic, "/", "EM_divB"))
                    {
                        auto& divB = modelView.getDivB();
                        writer.writeField(divB, layout);
                    }
            };

            modelView.visitHierarchy(write_quantity, ilvl, ilvl);
        },
        this->h5Writer_.minLevel, this->h5Writer_.maxLevel);
}


} // namespace PHARE::diagnostic::vtkh5

#endif /* PHARE_DIAGNOSTIC_DETAIL_VTK_TYPES_ELECTROMAG_HPP */
