#ifndef PHARE_CORE_NUMERICS_CONSTRAINED_TRANSPORT_POYNTING_CORRECTION_HPP
#define PHARE_CORE_NUMERICS_CONSTRAINED_TRANSPORT_POYNTING_CORRECTION_HPP

#include "core/utilities/types.hpp"
#include "core/utilities/index/index.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/data/vecfield/vecfield_component.hpp"
#include "core/numerics/constrained_transport/upwind_edge_b.hpp"

#include <type_traits>

namespace PHARE::core
{
template<typename GridLayout, template<typename> typename Reconstruction>
class PoyntingCorrection
{
    constexpr static auto dimension = GridLayout::dimension;
    using EdgeB_t                   = UpwindEdgeB<GridLayout, Reconstruction>;

public:
    explicit PoyntingCorrection(GridLayout const& layout)
        : layout_{layout}
    {
    }

    void operator()(auto const& ct_state, auto const& state, auto& fluxes) const
    {
        // Apply Poynting flux correction to energy: ∂E/∂t -= ∇·(E×B)
        // Because E from CT includes resistive (ηJ) and hyper-resistive terms, E×B
        // automatically captures the resistive Poynting flux.
        for_N<dimension>([&](auto i) {
            constexpr auto direction = static_cast<Direction>(i());

            layout_.evalOnBox(fluxes.template expose_centering<direction>(), [&](auto&... indices) {
                auto& F_Etot = fluxes.template get_dir<direction>({indices...}).Etot();
                poynting_energy_flux_<direction>(ct_state, state.E, state.B,
                                                 MeshIndex<dimension>{indices...}, F_Etot);
            });
        });
    }

private:
    template<auto direction>
    void poynting_energy_flux_(auto const& ct_state, auto const& E, auto const& B,
                               MeshIndex<dimension> const& index, auto& F_Etot) const
    {
        // Magnetic energy flux through the face: S.n = (E x B).n. Since S_i = eps_ijk E_j B_k,
        // a face never needs its own E component, which leaves exactly two products per face.

        // Same as GridLayout::project, but for an expression rather than a stored field: the
        // integrand is a product of two fields, and since both are point values the product
        // must be formed on the edge before being projected onto the face. Projecting the
        // factors separately and multiplying on the face would be a different (and, beyond
        // order 2, wrong) discretisation.
        auto project_on_face = [&](auto const& stencil, auto&& integrand) {
            std::decay_t<decltype(integrand(index))> result = 0.;

            for (auto const& wp : stencil)
                result += wp.coef * integrand(index + wp.indexes);

            return result;
        };

        auto const& Ex = E(Component::X);
        auto const& Ey = E(Component::Y);
        auto const& Ez = E(Component::Z);

        using implT = GridLayout::implT;

        if constexpr (direction == Direction::X)
        {
            constexpr auto zEdgeToXFace = implT::edgeZToFaceX();
            constexpr auto yEdgeToXFace = implT::edgeYToFaceX();

            auto const EzBy = project_on_face(zEdgeToXFace, [&](auto const& i) {
                return Ez(i) * EdgeB_t::template at<Component::Z, Component::Y>(ct_state, B, i);
            });
            auto const EyBz = project_on_face(yEdgeToXFace, [&](auto const& i) {
                return Ey(i) * EdgeB_t::template at<Component::Y, Component::Z>(ct_state, B, i);
            });

            F_Etot += EyBz - EzBy;
        }
        else if constexpr (direction == Direction::Y)
        {
            constexpr auto zEdgeToYFace = implT::edgeZToFaceY();
            constexpr auto xEdgeToYFace = implT::edgeXToFaceY();

            auto const EzBx = project_on_face(zEdgeToYFace, [&](auto const& i) {
                return Ez(i) * EdgeB_t::template at<Component::Z, Component::X>(ct_state, B, i);
            });
            auto const ExBz = project_on_face(xEdgeToYFace, [&](auto const& i) {
                return Ex(i) * EdgeB_t::template at<Component::X, Component::Z>(ct_state, B, i);
            });

            F_Etot += EzBx - ExBz;
        }
        else if constexpr (direction == Direction::Z)
        {
            constexpr auto xEdgeToZFace = implT::edgeXToFaceZ();
            constexpr auto yEdgeToZFace = implT::edgeYToFaceZ();

            auto const ExBy = project_on_face(xEdgeToZFace, [&](auto const& i) {
                return Ex(i) * EdgeB_t::template at<Component::X, Component::Y>(ct_state, B, i);
            });
            auto const EyBx = project_on_face(yEdgeToZFace, [&](auto const& i) {
                return Ey(i) * EdgeB_t::template at<Component::Y, Component::X>(ct_state, B, i);
            });

            F_Etot += ExBy - EyBx;
        }
    }

    GridLayout layout_;
};

} // namespace PHARE::core

#endif
