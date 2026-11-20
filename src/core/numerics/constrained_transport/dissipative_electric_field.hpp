#ifndef PHARE_CORE_NUMERICS_CONSTRAINED_TRANSPORT_DISSIPATIVE_ELECTRIC_FIELD_HPP
#define PHARE_CORE_NUMERICS_CONSTRAINED_TRANSPORT_DISSIPATIVE_ELECTRIC_FIELD_HPP

#include "core/def.hpp"
#include "core/numerics/ohm/ohm.hpp"
#include "core/numerics/differential_operators.hpp"
#include "core/utilities/index/index.hpp"
#include "core/utilities/meta/meta_utilities.hpp"
#include "core/models/quantities/mhd_quantities.hpp"
#include "core/data/vecfield/vecfield_component.hpp"

#include <cmath>
#include <vector>
#include <algorithm>

namespace PHARE::core
{

template<typename VecField>
class DissipativeElectricFieldState
{
public:
    DissipativeElectricFieldState() = default;
    explicit DissipativeElectricFieldState(bool const isDissipative)
    {
        if (isDissipative)
            E_.emplace_back("E_diss", MHDQuantity::Vector::E);
    }

    NO_DISCARD std::vector<VecField>& getRunTimeResourcesViewList() { return E_; }
    NO_DISCARD std::vector<VecField> const& getRunTimeResourcesViewList() const { return E_; }

    NO_DISCARD auto& E() { return E_[0]; }
    NO_DISCARD auto const& E() const { return E_[0]; }

private:
    std::vector<VecField> E_;
};


template<typename GridLayout>
class DissipativeElectricField : public OhmInfo
{
    using Super                     = OhmInfo;
    constexpr static auto dimension = GridLayout::dimension;

public:
    using Info_t = Super;

    DissipativeElectricField(OhmInfo const& info, GridLayout const& layout)
        : Super{info}
        , layout_{layout}
    {
    }

    void operator()(auto& dissipative_electric_state, auto const& mhd_state) const
    {
        auto& E         = dissipative_electric_state.E();
        auto const& J   = mhd_state.J;
        auto const& B   = mhd_state.B;
        auto const& rho = mhd_state.rho;

        auto const& meshSize = layout_.meshSize();
        auto const dx        = *std::min_element(meshSize.begin(), meshSize.end());
        auto const nu_dx2    = nu * dx * dx;

        Constexprifier{isResistive(), isHyperResistive(),
                       hyper_mode}([&]<bool isResistive, bool isHyperResistive, HyperMode hyper>() {
            for_N<3>([&](auto i) {
                constexpr auto component = static_cast<Component>(i());
                auto& Ec                 = E(component);
                auto const& Jc           = J(component);
                layout_.evalOnLines(Ec, [&](MeshIndex<dimension> const& first,
                                            std::uint32_t const n) {
                    double* ec       = &Ec(first);
                    double const* jc = &Jc(first);

                    if constexpr (isHyperResistive)
                    {
                        // laplacian(layout_, Jc, .) on the line: Jc one point before and after
                        // along each direction
                        std::array<double const*, dimension> jm, jp;
                        std::array<double, dimension> inv2;
                        for_N<dimension>([&](auto d) {
                            auto before = first, after = first;
                            before[d] -= 1;
                            after[d] += 1;
                            jm[d]          = &Jc(before);
                            jp[d]          = &Jc(after);
                            auto const inv = layout_.inverseMeshSize(static_cast<Direction>(d()));
                            inv2[d]        = inv * inv;
                        });
                        auto const coef = hyper_coef_<component, hyper>(B, rho, nu_dx2, first);

                        PHARE_VECTORIZE_IVDEP
                        for (std::uint32_t k = 0; k < n; ++k)
                        {
                            double e = 0.;
                            if constexpr (isResistive)
                                e += eta * jc[k];
                            double lap = 0.;
                            for (std::size_t d = 0; d < dimension; ++d)
                                lap += inv2[d] * (jp[d][k] - 2.0 * jc[k] + jm[d][k]);
                            e -= coef(k) * lap;
                            ec[k] = e;
                        }
                    }
                    else
                    {
                        PHARE_VECTORIZE_IVDEP
                        for (std::uint32_t k = 0; k < n; ++k)
                        {
                            double e = 0.;
                            if constexpr (isResistive)
                                e += eta * jc[k];
                            ec[k] = e;
                        }
                    }
                });
            });
        });
    }

private:
    // the hyper-resistive coefficient on the line from `first`, as coef(k);
    // nu_dx2 = nu dx^2, dx the smallest mesh size
    template<auto component, HyperMode hyper>
    auto hyper_coef_(auto const& B, auto const& rho, double const nu_dx2,
                     MeshIndex<dimension> const& first) const
    {
        if constexpr (hyper == HyperMode::constant)
            return [nu = nu](std::uint32_t) PHARE_LAMBDA_ALWAYS_INLINE { return nu; };
        else
        {
            auto coef = [&]<auto dst>() {
                return [bx = GridLayout::template projectLine<MHDQuantity::Scalar::Bx, dst>(
                            B[Component::X], first),
                        by = GridLayout::template projectLine<MHDQuantity::Scalar::By, dst>(
                            B[Component::Y], first),
                        bz = GridLayout::template projectLine<MHDQuantity::Scalar::Bz, dst>(
                            B[Component::Z], first),
                        n = GridLayout::template projectLine<MHDQuantity::Scalar::rho, dst>(rho,
                                                                                            first),
                        nu_dx2](std::uint32_t const k) PHARE_LAMBDA_ALWAYS_INLINE {
                    auto const bxk = bx(k), byk = by(k), bzk = bz(k);
                    auto const b = std::sqrt(bxk * bxk + byk * byk + bzk * bzk);
                    return nu_dx2 * (b / n(k) + 1);
                };
            };

            if constexpr (component == Component::X)
                return coef.template operator()<MHDQuantity::Scalar::Ex>();
            else if constexpr (component == Component::Y)
                return coef.template operator()<MHDQuantity::Scalar::Ey>();
            else
                return coef.template operator()<MHDQuantity::Scalar::Ez>();
        }
    }

    GridLayout layout_;
};

} // namespace PHARE::core

#endif
