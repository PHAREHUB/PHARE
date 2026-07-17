#ifndef PHARE_OHM_HPP
#define PHARE_OHM_HPP


#include "core/data/grid/grid_tiles.hpp"
#include "core/data/tensorfield/tensorfield.hpp"
#include "core/def.hpp"
#include "core/utilities/index/index.cpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/data/vecfield/vecfield_component.hpp"

#include "initializer/data_provider.hpp"


namespace PHARE::core
{
enum class HyperMode : std::uint8_t { constant = 0, spatial, LAST };


struct OhmInfo
{
    double const eta;
    double const nu;
    HyperMode const hyper_mode = HyperMode::spatial;

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
                    VecField const& J, VecField& Enew) const _PHARE_ALL_FN_
    {
        auto const& [Exnew, Eynew, Eznew] = Enew();

        layout_.evalOnBox(
            Exnew, [] _PHARE_ALL_FN_(auto&&... args) { E_Eq_<Component::X>(args...); }, n, Ve, Pe,
            B, J, Enew, *this);
        layout_.evalOnBox(
            Eynew, [] _PHARE_ALL_FN_(auto&&... args) { E_Eq_<Component::Y>(args...); }, n, Ve, Pe,
            B, J, Enew, *this);
        layout_.evalOnBox(
            Eznew, [] _PHARE_ALL_FN_(auto&&... args) { E_Eq_<Component::Z>(args...); }, n, Ve, Pe,
            B, J, Enew, *this);
    }

private:
    GridLayout layout_; // GPU copy needs object not pointer!


    template<auto Tag>
    static void E_Eq_(auto const& ijk, auto&&... args) _PHARE_ALL_FN_
    {
        auto const& [n, Ve, Pe, B, J, E, self] = std::forward_as_tuple(args...);
        auto& Exyz                             = E(Tag);

        static_assert(Components::check<Tag>());

        Exyz(ijk) = self.template ideal_<Tag>(Ve, B, ijk)      //
                    + self.template pressure_<Tag>(n, Pe, ijk) //
                    + self.template resistive_<Tag>(J, ijk)    //
                    + self.template hyperresistive_<Tag>(J, B, n, ijk);
    }




    template<auto component, typename VecField>
    auto ideal_(VecField const& Ve, VecField const& B,
                MeshIndex<dimension> index) const _PHARE_ALL_FN_
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
    auto pressure_(Field const& n, Field const& Pe, MeshIndex<dimension> index) const _PHARE_ALL_FN_
    {
        if constexpr (component == Component::X)
        {
            auto const nOnEx = GridLayout::template project<GridLayout::momentsToEx>(n, index);

            auto gradPOnEx = layout_.template deriv<Direction::X>(Pe, index); // TODO : issue 3391

            return -gradPOnEx / nOnEx;
        }

        else if constexpr (component == Component::Y)
        {
            if constexpr (dimension >= 2)
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
            if constexpr (dimension >= 3)
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
    auto resistive_(VecField const& J, MeshIndex<dimension> index) const _PHARE_ALL_FN_
    {
        auto const& Jxyx = J(component);

        if constexpr (component == Component::X)
        {
            auto const jxOnEx = GridLayout::template project<GridLayout::JxToEx>(Jxyx, index);
            return eta * jxOnEx;
        }

        if constexpr (component == Component::Y)
        {
            auto const jyOnEy = GridLayout::template project<GridLayout::JyToEy>(Jxyx, index);
            return eta * jyOnEy;
        }

        if constexpr (component == Component::Z)
        {
            auto const jzOnEz = GridLayout::template project<GridLayout::JzToEz>(Jxyx, index);
            return eta * jzOnEz;
        }
    }

    template<auto component, typename VecField, typename Field>
    auto hyperresistive_(VecField const& J, VecField const& B, Field const& n,
                         MeshIndex<VecField::dimension> index) const _PHARE_ALL_FN_
    {
        // if compile error, fix this function
        static_assert(static_cast<std::underlying_type_t<HyperMode>>(HyperMode::LAST) == 2);

        if (hyper_mode == HyperMode::constant)
            return constant_hyperresistive_<component>(J, index);

        // else
        return spatial_hyperresistive_<component>(J, B, n, index);
    }


    template<auto component, typename VecField>
    auto constant_hyperresistive_(VecField const& J,
                                  MeshIndex<VecField::dimension> index) const _PHARE_ALL_FN_
    { // TODO : https://github.com/PHAREHUB/PHARE/issues/3
        return -nu * layout_.laplacian(J(component), index);
    }


    template<auto component, typename VecField, typename Field>
    auto spatial_hyperresistive_(VecField const& J, VecField const& B, Field const& n,
                                 MeshIndex<VecField::dimension> index) const
    {
        auto const lvlCoeff = 1. / std::pow(4, layout_.levelNumber());

        auto constexpr min_density = 0.1;
        auto computeHR             = [&]<auto BxProj, auto ByProj, auto BzProj, auto nProj>() {
            auto const BxOnE = GridLayout::template project<BxProj>(B(Component::X), index);
            auto const ByOnE = GridLayout::template project<ByProj>(B(Component::Y), index);
            auto const BzOnE = GridLayout::template project<BzProj>(B(Component::Z), index);
            auto const nOnE  = GridLayout::template project<nProj>(n, index);
            auto b           = std::sqrt(BxOnE * BxOnE + ByOnE * ByOnE + BzOnE * BzOnE);

            return -nu * (b / (nOnE + min_density) + 1) * lvlCoeff
                   * layout_.laplacian(J(component), index);
        };
        if constexpr (component == Component::X)
            return computeHR.template operator()<GridLayout::BxToEx, GridLayout::ByToEx,
                                                 GridLayout::BzToEx, GridLayout::momentsToEx>();

        if constexpr (component == Component::Y)
            return computeHR.template operator()<GridLayout::BxToEy, GridLayout::ByToEy,
                                                 GridLayout::BzToEy, GridLayout::momentsToEy>();

        if constexpr (component == Component::Z)
            return computeHR.template operator()<GridLayout::BxToEz, GridLayout::ByToEz,
                                                 GridLayout::BzToEz, GridLayout::momentsToEz>();
    }
};


class OhmSingleTransformer
{
    using info_type = OhmInfo;

    template<typename V_t>
    V_t static tt(auto& vf, auto i)
    {
        return vf.template as<V_t>([&](auto& c) { return c()[i](); });
    }

public:
    explicit OhmSingleTransformer(info_type const& info)
        : info_{info}
    {
    }

    template<typename GridLayout, typename VecField, typename Field>
    void operator()(GridLayout const& layout, Field const& n, VecField const& Ve, Field const& Pe,
                    VecField const& B, VecField const& J, VecField& Enew)
    {
        PHARE_LOG_SCOPE(2, "OhmLevelTransformer");

        // check_field(n);
        // check_field(Pe);
        // check_tensor_field(Ve, layout);
        // check_tensor_field(B);
        // check_tensor_field(J);

        if constexpr (is_field_tile_set_v<Field>)
        {
            using Tile_vt = Field::value_type::value_type;
            static_assert(is_field_v<Tile_vt>);
            using V_t = basic::TensorField<Tile_vt, 1>;

            for (std::size_t tidx = 0; tidx < n().size(); ++tidx)
            {
                auto Enw             = Enew.template as<V_t>([&](auto& c) { return c()[tidx](); });
                auto const& tile_lay = n()[tidx].layout();
                using TL             = std::remove_cvref_t<decltype(tile_lay)>;
                Ohm<TL>{info_, tile_lay}(n()[tidx](), tt<V_t>(Ve, tidx), Pe()[tidx](),
                                         tt<V_t>(B, tidx), tt<V_t>(J, tidx), Enw);
            }
            PHARE_LOG_SCOPE(2, "OhmLevelTransformer::sync");
            for (std::uint8_t i = 0; i < 3; ++i)
                Enew[i].sync_inner_ghosts();
        }
        else
        {
            Ohm<GridLayout>{info_, layout}(n, Ve, Pe, B, J, Enew);
        }

        check_tensor_field(Enew, layout);
    }


    info_type info_;
};


} // namespace PHARE::core
#endif
