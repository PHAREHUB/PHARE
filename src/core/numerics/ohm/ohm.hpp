#ifndef PHARE_OHM_HPP
#define PHARE_OHM_HPP


#include "core/utilities/index/index.hpp"
#include "core/utilities/meta/meta_utilities.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/data/vecfield/vecfield_component.hpp"
#include "core/models/quantities/hybrid_quantities.hpp"

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
    using Q                         = HybridQuantity::Scalar;

public:
    explicit Ohm(OhmInfo const& info, GridLayout const& layout)
        : Super{info}
        , layout_{layout}
    {
    }

    template<typename VecField, typename Field>
    void operator()(Field const& n, VecField const& Ve, Field const& Pe, VecField const& B,
                    VecField const& J, VecField& Enew)
    {
        // lift the resistive / hyper-resistive runtime flags into compile-time tags: the per-cell
        // E_Eq_ branches only via if constexpr, skipping the projection / laplacian when eta or
        // nu is zero without ever testing the flags inside the evalOnBox loop.
        Constexprifier{isResistive(),
                       isHyperResistive()}([&]<bool isResistiveV, bool isHyperResistiveV>() {
            solve_<isResistiveV, isHyperResistiveV>(n, Ve, Pe, B, J, Enew);
        });
    }

private:
    GridLayout layout_;

    template<bool isResistive, bool isHyperResistive, typename VecField, typename Field>
    void solve_(Field const& n, VecField const& Ve, Field const& Pe, VecField const& B,
                VecField const& J, VecField& Enew) const
    {
        using Pack = OhmPack<VecField, Field>;

        auto const& [Exnew, Eynew, Eznew] = Enew();

        layout_.evalOnBox(Exnew, [&](auto&... args) {
            this->template E_Eq_<Component::X, isResistive, isHyperResistive>(
                Pack{Enew, n, Pe, Ve, B, J}, args...);
        });
        layout_.evalOnBox(Eynew, [&](auto&... args) {
            this->template E_Eq_<Component::Y, isResistive, isHyperResistive>(
                Pack{Enew, n, Pe, Ve, B, J}, args...);
        });
        layout_.evalOnBox(Eznew, [&](auto&... args) {
            this->template E_Eq_<Component::Z, isResistive, isHyperResistive>(
                Pack{Enew, n, Pe, Ve, B, J}, args...);
        });
    }

    template<typename VecField, typename Field>
    struct OhmPack
    {
        VecField& Exyz;
        Field const &n, &Pe;
        VecField const &Ve, &B, &J;
    };


    template<auto Tag, bool isResistive, bool isHyperResistive, typename OhmPack, typename... IDXs>
    void E_Eq_(OhmPack&& pack, IDXs const&... ijk) const
    {
        auto const& [E, n, Pe, Ve, B, J] = pack;
        auto& Exyz                       = E(Tag);

        static_assert(Components::check<Tag>());

        auto E_ = ideal_<Tag>(Ve, B, {ijk...}) + pressure_<Tag>(n, Pe, {ijk...});

        if constexpr (isResistive)
            E_ += resistive_<Tag>(J, {ijk...});
        if constexpr (isHyperResistive)
            E_ += hyperresistive_<Tag>(J, B, n, {ijk...});

        Exyz(ijk...) = E_;
    }



    template<auto component, typename VecField>
    auto ideal_(VecField const& Ve, VecField const& B, MeshIndex<dimension> index) const
    {
        if constexpr (component == Component::X)
        {
            auto const& Vy = Ve(Component::Y);
            auto const& Vz = Ve(Component::Z);

            auto const& By = B(Component::Y);
            auto const& Bz = B(Component::Z);

            auto const vyOnEx = GridLayout::template project<Q::Vy, Q::Ex>(Vy, index);
            auto const vzOnEx = GridLayout::template project<Q::Vz, Q::Ex>(Vz, index);
            auto const byOnEx = GridLayout::template project<Q::By, Q::Ex>(By, index);
            auto const bzOnEx = GridLayout::template project<Q::Bz, Q::Ex>(Bz, index);

            return -vyOnEx * bzOnEx + vzOnEx * byOnEx;
        }

        if constexpr (component == Component::Y)
        {
            auto const& Vx = Ve(Component::X);
            auto const& Vz = Ve(Component::Z);
            auto const& Bx = B(Component::X);
            auto const& Bz = B(Component::Z);

            auto const vxOnEy = GridLayout::template project<Q::Vx, Q::Ey>(Vx, index);
            auto const vzOnEy = GridLayout::template project<Q::Vz, Q::Ey>(Vz, index);
            auto const bxOnEy = GridLayout::template project<Q::Bx, Q::Ey>(Bx, index);
            auto const bzOnEy = GridLayout::template project<Q::Bz, Q::Ey>(Bz, index);

            return -vzOnEy * bxOnEy + vxOnEy * bzOnEy;
        }

        if constexpr (component == Component::Z)
        {
            auto const& Vx = Ve(Component::X);
            auto const& Vy = Ve(Component::Y);
            auto const& Bx = B(Component::X);
            auto const& By = B(Component::Y);

            auto const vxOnEz = GridLayout::template project<Q::Vx, Q::Ez>(Vx, index);
            auto const vyOnEz = GridLayout::template project<Q::Vy, Q::Ez>(Vy, index);
            auto const bxOnEz = GridLayout::template project<Q::Bx, Q::Ez>(Bx, index);
            auto const byOnEz = GridLayout::template project<Q::By, Q::Ez>(By, index);

            return -vxOnEz * byOnEz + vyOnEz * bxOnEz;
        }
    }


    template<auto component, typename Field>
    auto pressure_(Field const& n, Field const& Pe, MeshIndex<Field::dimension> index) const
    {
        if constexpr (component == Component::X)
        {
            auto const nOnEx = GridLayout::template project<Q::rho, Q::Ex>(n, index);

            auto gradPOnEx = layout_.template deriv<Direction::X>(Pe, index); // TODO : issue 3391

            return -gradPOnEx / nOnEx;
        }

        else if constexpr (component == Component::Y)
        {
            if constexpr (Field::dimension >= 2)
            {
                auto const nOnEy = GridLayout::template project<Q::rho, Q::Ey>(n, index);

                auto gradPOnEy
                    = layout_.template deriv<Direction::Y>(Pe, index); // TODO : issue 3391

                return -gradPOnEy / nOnEy;
            }
            else
            {
                return 0.;
            }
        }

        else if constexpr (component == Component::Z)
        {
            if constexpr (Field::dimension >= 3)
            {
                auto const nOnEz = GridLayout::template project<Q::rho, Q::Ez>(n, index);

                auto gradPOnEz
                    = layout_.template deriv<Direction::Z>(Pe, index); // TODO : issue 3391

                return -gradPOnEz / nOnEz;
            }
            else
            {
                return 0.;
            }
        }
    }




    template<auto component, typename VecField>
    auto resistive_(VecField const& J, MeshIndex<VecField::dimension> index) const
    {
        auto const& Jxyx = J(component);

        if constexpr (component == Component::X)
        {
            auto const jxOnEx = GridLayout::template project<Q::Jx, Q::Ex>(Jxyx, index);
            return eta * jxOnEx;
        }

        if constexpr (component == Component::Y)
        {
            auto const jyOnEy = GridLayout::template project<Q::Jy, Q::Ey>(Jxyx, index);
            return eta * jyOnEy;
        }

        if constexpr (component == Component::Z)
        {
            auto const jzOnEz = GridLayout::template project<Q::Jz, Q::Ez>(Jxyx, index);
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
        return -nu * layout_.laplacian(J(component), index);
    }


    template<auto component, typename VecField, typename Field>
    auto spatial_hyperresistive_(VecField const& J, VecField const& B, Field const& n,
                                 MeshIndex<VecField::dimension> index) const
    {
        auto const lvlCoeff        = 1. / std::pow(4, layout_.levelNumber());
        auto constexpr min_density = 0.1;
        auto computeHR             = [&]<auto dst>() {
            auto const BxOnE = GridLayout::template project<Q::Bx, dst>(B(Component::X), index);
            auto const ByOnE = GridLayout::template project<Q::By, dst>(B(Component::Y), index);
            auto const BzOnE = GridLayout::template project<Q::Bz, dst>(B(Component::Z), index);
            auto const nOnE  = GridLayout::template project<Q::rho, dst>(n, index);
            auto b           = std::sqrt(BxOnE * BxOnE + ByOnE * ByOnE + BzOnE * BzOnE);
            return -nu * (b / (nOnE + min_density) + 1) * lvlCoeff
                   * layout_.laplacian(J(component), index);
        };
        if constexpr (component == Component::X)
        {
            return computeHR.template operator()<Q::Ex>();
        }
        if constexpr (component == Component::Y)
        {
            return computeHR.template operator()<Q::Ey>();
        }
        if constexpr (component == Component::Z)
        {
            return computeHR.template operator()<Q::Ez>();
        }
    }
};


} // namespace PHARE::core
#endif
