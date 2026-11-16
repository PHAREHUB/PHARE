#ifndef PHARE_OHM_HPP
#define PHARE_OHM_HPP


#include "core/utilities/index/index.hpp"
#include "core/utilities/meta/enum_meta.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/data/vecfield/vecfield_component.hpp"
#include "core/numerics/differential_operators.hpp"

#include "initializer/data_provider.hpp"


namespace PHARE::core
{

// `count` closes the value range, so Constexprifier can fan out every case (see CountedEnum)
enum class HyperMode { constant, spatial, count };

struct OhmInfo
{
    double const eta;
    double const nu;
    HyperMode const hyper_mode;

    bool isResistive() const { return eta > 0.0; }
    bool isHyperResistive() const { return nu > 0.0; }

    OhmInfo static FROM(initializer::PHAREDict const& dict)
    {
        return {dict["resistivity"].template to<double>(),
                dict["hyper_resistivity"].template to<double>(),
                cppdict::get_value(dict, "hyper_mode", HyperMode::constant)};
    }
};


template<typename GridLayout>
class Ohm : public OhmInfo
{
    using Super                     = OhmInfo;
    constexpr static auto dimension = GridLayout::dimension;
    using Scalar                    = GridLayout::Scalar;

public:
    explicit Ohm(OhmInfo const& info, GridLayout const& layout)
        : Super{info}
        , layout_{layout}
    {
    }

    // n: ion charge density, F: total ion charge flux sum_s q_s F_s (moment nodes).
    // pressureTerm: -grad(Pe)/n on the E edges, from the electron pressure closure (see
    // pressure_term.hpp).
    template<typename VecField, typename Field>
    void operator()(Field const& n, VecField const& F, VecField const& pressureTerm,
                    VecField const& B, VecField const& J, VecField& Enew)
    {
        for (auto const c : {Component::X, Component::Y, Component::Z})
            if (layout_.centering(pressureTerm(c)) != layout_.centering(Enew(c)))
                throw std::runtime_error("Error - Ohm - the pressure term must be E-centred");

        // lift the resistive / hyper-resistive runtime flags into compile-time tags: the per-cell
        // E_Eq_ branches only via if constexpr, skipping the projection / laplacian when eta or
        // nu is zero without ever testing the flags inside the evalOnBox loop.
        Constexprifier{isResistive(),
                       isHyperResistive()}([&]<bool isResistiveV, bool isHyperResistiveV>() {
            solve_<isResistiveV, isHyperResistiveV>(n, F, pressureTerm, B, J, Enew);
        });
    }

private:
    GridLayout layout_;

    template<bool isResistive, bool isHyperResistive, typename VecField, typename Field>
    void solve_(Field const& n, VecField const& F, VecField const& pressureTerm, VecField const& B,
                VecField const& J, VecField& Enew) const
    {
        using Pack = OhmPack<VecField, Field>;

        auto const& [Exnew, Eynew, Eznew] = Enew();

        layout_.evalOnBox(Exnew, [&](auto&... args) {
            this->template E_Eq_<Component::X, isResistive, isHyperResistive>(
                Pack{Enew, n, F, B, J}, pressureTerm(Component::X), args...);
        });
        layout_.evalOnBox(Eynew, [&](auto&... args) {
            this->template E_Eq_<Component::Y, isResistive, isHyperResistive>(
                Pack{Enew, n, F, B, J}, pressureTerm(Component::Y), args...);
        });
        layout_.evalOnBox(Eznew, [&](auto&... args) {
            this->template E_Eq_<Component::Z, isResistive, isHyperResistive>(
                Pack{Enew, n, F, B, J}, pressureTerm(Component::Z), args...);
        });
    }

    template<typename VecField, typename Field>
    struct OhmPack
    {
        VecField& Exyz;
        Field const& n;
        VecField const &F, &B, &J;
    };


    template<auto Tag, bool isResistive, bool isHyperResistive, typename OhmPack,
             typename PressureTerm, typename... IDXs>
    void E_Eq_(OhmPack&& pack, PressureTerm const& pressureTerm, IDXs const&... ijk) const
    {
        auto const& [E, n, F, B, J] = pack;
        auto& Exyz                  = E(Tag);

        static_assert(Components::check<Tag>());

        auto E_ = ideal_<Tag>(n, F, B, J, {ijk...}) + pressureTerm(ijk...);

        if constexpr (isResistive)
            E_ += resistive_<Tag>(J, {ijk...});
        if constexpr (isHyperResistive)
            E_ += hyperresistive_<Tag>(J, B, n, {ijk...});

        Exyz(ijk...) = E_;
    }



    // The ideal term -Ve x B = -G x B / n, with G = n Ve = F - J the electron flux, written so
    // that it does no work on the electrons: sum_i G_i E_i = 0 over the periodic grid (#1577).
    // G is formed on the E nodes, B/n on the B nodes, and each product is averaged to the B node
    // and back: E_x = -P_{Bz->Ex}(wz P_{Ey->Bz} Gy) + P_{By->Ex}(wy P_{Ez->By} Gz), w = B/n.
    // The 2-point averages P_{a->b} and P_{b->a} are transposes of each other, so the operator
    // G -> E is antisymmetric. Averaging Ve and B separately to the E node and multiplying there
    // (the former form) is not, and heats the plasma from particle noise.
    // Reads F and n one node beyond the box E is computed on.
    template<auto component, typename VecField, typename Field>
    auto ideal_(Field const& n, VecField const& F, VecField const& B, VecField const& J,
                MeshIndex<dimension> index) const
    {
        auto G = [&]<auto c, auto Fq, auto Eq>(auto const& ijk) {
            return GridLayout::template project<Fq, Eq>(F(c), ijk) - J(c)(ijk);
        };
        auto Gx = [&](auto const& ijk) {
            return G.template operator()<Component::X, Scalar::Vx, Scalar::Ex>(ijk);
        };
        auto Gy = [&](auto const& ijk) {
            return G.template operator()<Component::Y, Scalar::Vy, Scalar::Ey>(ijk);
        };
        auto Gz = [&](auto const& ijk) {
            return G.template operator()<Component::Z, Scalar::Vz, Scalar::Ez>(ijk);
        };

        // B_c / n on the B_c node, times the electron flux component Gq averaged to it
        auto omegaG = [&]<auto c, auto Bq, auto Eq>(auto const& Gq, auto const& ijk) {
            auto const nOnB = GridLayout::template project<Scalar::rho, Bq>(n, ijk);
            return B(c)(ijk) / nOnB * GridLayout::template project<Eq, Bq>(Gq, ijk);
        };
        auto wxG = [&]<auto Eq>(auto const& Gq) {
            return [&](auto const& ijk) {
                return omegaG.template operator()<Component::X, Scalar::Bx, Eq>(Gq, ijk);
            };
        };
        auto wyG = [&]<auto Eq>(auto const& Gq) {
            return [&](auto const& ijk) {
                return omegaG.template operator()<Component::Y, Scalar::By, Eq>(Gq, ijk);
            };
        };
        auto wzG = [&]<auto Eq>(auto const& Gq) {
            return [&](auto const& ijk) {
                return omegaG.template operator()<Component::Z, Scalar::Bz, Eq>(Gq, ijk);
            };
        };

        if constexpr (component == Component::X)
            return -GridLayout::template project<Scalar::Bz, Scalar::Ex>(
                       wzG.template operator()<Scalar::Ey>(Gy), index)
                   + GridLayout::template project<Scalar::By, Scalar::Ex>(
                       wyG.template operator()<Scalar::Ez>(Gz), index);

        if constexpr (component == Component::Y)
            return GridLayout::template project<Scalar::Bz, Scalar::Ey>(
                       wzG.template operator()<Scalar::Ex>(Gx), index)
                   - GridLayout::template project<Scalar::Bx, Scalar::Ey>(
                       wxG.template operator()<Scalar::Ez>(Gz), index);

        if constexpr (component == Component::Z)
            return -GridLayout::template project<Scalar::By, Scalar::Ez>(
                       wyG.template operator()<Scalar::Ex>(Gx), index)
                   + GridLayout::template project<Scalar::Bx, Scalar::Ez>(
                       wxG.template operator()<Scalar::Ey>(Gy), index);
    }


    template<auto component, typename VecField>
    auto resistive_(VecField const& J, MeshIndex<VecField::dimension> index) const
    {
        auto const& Jxyx = J(component);

        if constexpr (component == Component::X)
        {
            auto const jxOnEx = GridLayout::template project<Scalar::Jx, Scalar::Ex>(Jxyx, index);
            return eta * jxOnEx;
        }

        if constexpr (component == Component::Y)
        {
            auto const jyOnEy = GridLayout::template project<Scalar::Jy, Scalar::Ey>(Jxyx, index);
            return eta * jyOnEy;
        }

        if constexpr (component == Component::Z)
        {
            auto const jzOnEz = GridLayout::template project<Scalar::Jz, Scalar::Ez>(Jxyx, index);
            return eta * jzOnEz;
        }
    }

    template<auto component, typename VecField, typename Field>
    auto hyperresistive_(VecField const& J, VecField const& B, Field const& n,
                         MeshIndex<VecField::dimension> index) const
    {
        if (hyper_mode == HyperMode::constant)
            return constant_hyperresistive_<component>(J, index);
        else if (hyper_mode == HyperMode::spatial)
            return spatial_hyperresistive_<component>(J, B, n, index);
        else // should not happen but otherwise -Wreturn-type fails with Werror
            throw std::runtime_error("Error - Ohm - unknown hyper_mode");
    }


    template<auto component, typename VecField>
    auto constant_hyperresistive_(VecField const& J, MeshIndex<VecField::dimension> index) const
    { // TODO : https://github.com/PHAREHUB/PHARE/issues/3
        return -nu * laplacian(layout_, J(component), index);
    }


    template<auto component, typename VecField, typename Field>
    auto spatial_hyperresistive_(VecField const& J, VecField const& B, Field const& n,
                                 MeshIndex<VecField::dimension> index) const
    {
        auto const lvlCoeff        = 1. / std::pow(4, layout_.levelNumber());
        auto constexpr min_density = 0.1;
        auto computeHR             = [&]<auto dst>() {
            auto const BxOnE
                = GridLayout::template project<Scalar::Bx, dst>(B(Component::X), index);
            auto const ByOnE
                = GridLayout::template project<Scalar::By, dst>(B(Component::Y), index);
            auto const BzOnE
                = GridLayout::template project<Scalar::Bz, dst>(B(Component::Z), index);
            auto const nOnE = GridLayout::template project<Scalar::rho, dst>(n, index);
            auto b          = std::sqrt(BxOnE * BxOnE + ByOnE * ByOnE + BzOnE * BzOnE);
            return -nu * (b / (nOnE + min_density) + 1) * lvlCoeff
                   * laplacian(layout_, J(component), index);
        };
        if constexpr (component == Component::X)
        {
            return computeHR.template operator()<Scalar::Ex>();
        }
        if constexpr (component == Component::Y)
        {
            return computeHR.template operator()<Scalar::Ey>();
        }
        if constexpr (component == Component::Z)
        {
            return computeHR.template operator()<Scalar::Ez>();
        }
    }
};


} // namespace PHARE::core
#endif
