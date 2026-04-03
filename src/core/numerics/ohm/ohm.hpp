#ifndef PHARE_OHM_HPP
#define PHARE_OHM_HPP


#include "core/utilities/index/index.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/data/vecfield/vecfield_component.hpp"

#include "initializer/data_provider.hpp"


namespace PHARE::core
{

enum class HyperMode { constant, spatial };

struct OhmInfo
{
    double const eta_;
    double const nu_;
    HyperMode const hyper_mode;

    OhmInfo static FROM(initializer::PHAREDict const& dict)
    {
        return {dict["resistivity"].template to<double>(),
                dict["hyper_resistivity"].template to<double>(),
                cppdict::get_value(dict, "hyper_mode", std::string{"constant"}) == "constant"
                    ? HyperMode::constant
                    : HyperMode::spatial};
    }
};

template<typename GridLayout>
class Ohm : public OhmInfo
{
    using Super                     = OhmInfo;
    constexpr static auto dimension = GridLayout::dimension;

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
        using Pack = OhmPack<VecField, Field>;

        auto const& [Exnew, Eynew, Eznew] = Enew();

        layout_.evalOnBox(Exnew, [&](auto&... args) mutable {
            this->template E_Eq_<Component::X>(Pack{Enew, n, Pe, Ve, B, J}, args...);
        });
        layout_.evalOnBox(Eynew, [&](auto&... args) mutable {
            this->template E_Eq_<Component::Y>(Pack{Enew, n, Pe, Ve, B, J}, args...);
        });
        layout_.evalOnBox(Eznew, [&](auto&... args) mutable {
            this->template E_Eq_<Component::Z>(Pack{Enew, n, Pe, Ve, B, J}, args...);
        });
    }

private:
    GridLayout layout_;

    template<typename VecField, typename Field>
    struct OhmPack
    {
        VecField& Exyz;
        Field const &n, &Pe;
        VecField const &Ve, &B, &J;
    };


    template<auto Tag, typename OhmPack, typename... IDXs>
    void E_Eq_(OhmPack&& pack, IDXs const&... ijk) const
    {
        auto const& [E, n, Pe, Ve, B, J] = pack;
        auto& Exyz                       = E(Tag);

        static_assert(Components::check<Tag>());

        Exyz(ijk...) = ideal_<Tag>(Ve, B, {ijk...})      //
                       + pressure_<Tag>(n, Pe, {ijk...}) //
                       + resistive_<Tag>(J, {ijk...})    //
                       + hyperresistive_<Tag>(J, B, n, {ijk...});
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

            auto const vyOnEx = GridLayout::template project<GridLayout::momentsToEx>(Vy, index);
            auto const vzOnEx = GridLayout::template project<GridLayout::momentsToEx>(Vz, index);
            auto const byOnEx = GridLayout::template project<GridLayout::ByToEx>(By, index);
            auto const bzOnEx = GridLayout::template project<GridLayout::BzToEx>(Bz, index);

            return -vyOnEx * bzOnEx + vzOnEx * byOnEx;
        }

        if constexpr (component == Component::Y)
        {
            auto const& Vx = Ve(Component::X);
            auto const& Vz = Ve(Component::Z);
            auto const& Bx = B(Component::X);
            auto const& Bz = B(Component::Z);

            auto const vxOnEy = GridLayout::template project<GridLayout::momentsToEy>(Vx, index);
            auto const vzOnEy = GridLayout::template project<GridLayout::momentsToEy>(Vz, index);
            auto const bxOnEy = GridLayout::template project<GridLayout::BxToEy>(Bx, index);
            auto const bzOnEy = GridLayout::template project<GridLayout::BzToEy>(Bz, index);

            return -vzOnEy * bxOnEy + vxOnEy * bzOnEy;
        }

        if constexpr (component == Component::Z)
        {
            auto const& Vx = Ve(Component::X);
            auto const& Vy = Ve(Component::Y);
            auto const& Bx = B(Component::X);
            auto const& By = B(Component::Y);

            auto const vxOnEz = GridLayout::template project<GridLayout::momentsToEz>(Vx, index);
            auto const vyOnEz = GridLayout::template project<GridLayout::momentsToEz>(Vy, index);
            auto const bxOnEz = GridLayout::template project<GridLayout::BxToEz>(Bx, index);
            auto const byOnEz = GridLayout::template project<GridLayout::ByToEz>(By, index);

            return -vxOnEz * byOnEz + vyOnEz * bxOnEz;
        }
    }


    template<auto component, typename Field>
    auto pressure_(Field const& n, Field const& Pe, MeshIndex<Field::dimension> index) const
    {
        if constexpr (component == Component::X)
        {
            auto const nOnEx = GridLayout::template project<GridLayout::momentsToEx>(n, index);

            auto gradPOnEx = layout_.template deriv<Direction::X>(Pe, index); // TODO : issue 3391

            return -gradPOnEx / nOnEx;
        }

        else if constexpr (component == Component::Y)
        {
            if constexpr (Field::dimension >= 2)
            {
                auto const nOnEy = GridLayout::template project<GridLayout::momentsToEy>(n, index);

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
                auto const nOnEz = GridLayout::template project<GridLayout::momentsToEz>(n, index);

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
            auto const jxOnEx = GridLayout::template project<GridLayout::JxToEx>(Jxyx, index);
            return eta_ * jxOnEx;
        }

        if constexpr (component == Component::Y)
        {
            auto const jyOnEy = GridLayout::template project<GridLayout::JyToEy>(Jxyx, index);
            return eta_ * jyOnEy;
        }

        if constexpr (component == Component::Z)
        {
            auto const jzOnEz = GridLayout::template project<GridLayout::JzToEz>(Jxyx, index);
            return eta_ * jzOnEz;
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
        return -nu_ * layout_.laplacian(J(component), index);
    }


    template<auto component, typename VecField, typename Field>
    auto spatial_hyperresistive_(VecField const& J, VecField const& B, Field const& n,
                                 MeshIndex<VecField::dimension> index) const
    {
        auto const lvlCoeff        = 1. / std::pow(4, layout_->levelNumber());
        auto constexpr min_density = 0.1;
        auto computeHR             = [&]<auto BxProj, auto ByProj, auto BzProj, auto nProj>() {
            auto const BxOnE = GridLayout::template project<BxProj>(B(Component::X), index);
            auto const ByOnE = GridLayout::template project<ByProj>(B(Component::Y), index);
            auto const BzOnE = GridLayout::template project<BzProj>(B(Component::Z), index);
            auto const nOnE  = GridLayout::template project<nProj>(n, index);
            auto b           = std::sqrt(BxOnE * BxOnE + ByOnE * ByOnE + BzOnE * BzOnE);
            return -nu_ * (b / (nOnE + min_density) + 1) * lvlCoeff
                   * layout_.laplacian(J(component), index);
        };
        if constexpr (component == Component::X)
        {
            return computeHR.template operator()<GridLayout::BxToEx, GridLayout::ByToEx,
                                                 GridLayout::BzToEx, GridLayout::momentsToEx>();
        }
        if constexpr (component == Component::Y)
        {
            return computeHR.template operator()<GridLayout::BxToEy, GridLayout::ByToEy,
                                                 GridLayout::BzToEy, GridLayout::momentsToEy>();
        }
        if constexpr (component == Component::Z)
        {
            return computeHR.template operator()<GridLayout::BxToEz, GridLayout::ByToEz,
                                                 GridLayout::BzToEz, GridLayout::momentsToEz>();
        }
    }
};


} // namespace PHARE::core
#endif
