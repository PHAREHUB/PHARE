#ifndef PHARE_OHM_HPP
#define PHARE_OHM_HPP


#include "core/def.hpp"
#include "core/utilities/index/index.hpp"
#include "core/utilities/meta/meta_utilities.hpp"
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
        // checked once here: the per-cell code below uses the unchecked operator[] accessors
        if (!isUsable(n, Ve, Pe, B, J, Enew))
            throw std::runtime_error("Error - Ohm - not all parameters are usable");

        // lift the resistive / hyper-resistive runtime flags into compile-time tags: the per-cell
        // E_Eq_ branches only via if constexpr, skipping the projection / laplacian when eta or
        // nu is zero without ever testing the flags inside the evalOnBox loop.
        Constexprifier{isResistive(), isHyperResistive(), hyper_mode}(
            [&]<bool isResistiveV, bool isHyperResistiveV, HyperMode hyperModeV>() {
                solve_<isResistiveV, isHyperResistiveV, hyperModeV>(n, Ve, Pe, B, J, Enew);
            });
    }

private:
    GridLayout layout_;

    template<bool isResistive, bool isHyperResistive, HyperMode hyperMode, typename VecField,
             typename Field>
    void solve_(Field const& n, VecField const& Ve, Field const& Pe, VecField const& B,
                VecField const& J, VecField& Enew) const
    {
        using Pack = OhmPack<VecField, Field>;

        auto const& [Exnew, Eynew, Eznew] = Enew();

        auto const gradPe = grad(layout_, Pe);

        layout_.evalOnBoxIndependent(Exnew, [&](auto&... args) PHARE_LAMBDA_ALWAYS_INLINE {
            this->template E_Eq_<Component::X, isResistive, isHyperResistive, hyperMode>(
                Pack{Enew, n, Ve, B, J}, gradPe.x, args...);
        });
        layout_.evalOnBoxIndependent(Eynew, [&](auto&... args) PHARE_LAMBDA_ALWAYS_INLINE {
            this->template E_Eq_<Component::Y, isResistive, isHyperResistive, hyperMode>(
                Pack{Enew, n, Ve, B, J}, gradPe.y, args...);
        });
        layout_.evalOnBoxIndependent(Eznew, [&](auto&... args) PHARE_LAMBDA_ALWAYS_INLINE {
            this->template E_Eq_<Component::Z, isResistive, isHyperResistive, hyperMode>(
                Pack{Enew, n, Ve, B, J}, gradPe.z, args...);
        });
    }

    template<typename VecField, typename Field>
    struct OhmPack
    {
        VecField& Exyz;
        Field const& n;
        VecField const &Ve, &B, &J;
    };


    template<auto Tag, bool isResistive, bool isHyperResistive, HyperMode hyperMode,
             typename OhmPack, typename GradPe, typename... IDXs>
    PHARE_ALWAYS_INLINE void E_Eq_(OhmPack&& pack, GradPe const& gradPe, IDXs const&... ijk) const
    {
        auto const& [E, n, Ve, B, J] = pack;
        auto& Exyz                   = E[Tag];

        static_assert(Components::check<Tag>());

        auto E_ = ideal_<Tag>(Ve, B, {ijk...}) + pressure_<Tag>(n, gradPe, {ijk...});

        if constexpr (isResistive)
            E_ += resistive_<Tag>(J, {ijk...});
        if constexpr (isHyperResistive)
            E_ += hyperresistive_<Tag, hyperMode>(J, B, n, {ijk...});

        Exyz(ijk...) = E_;
    }



    template<auto component, typename VecField>
    PHARE_ALWAYS_INLINE auto ideal_(VecField const& Ve, VecField const& B,
                                    MeshIndex<dimension> index) const
    {
        if constexpr (component == Component::X)
        {
            auto const& Vy = Ve[Component::Y];
            auto const& Vz = Ve[Component::Z];

            auto const& By = B[Component::Y];
            auto const& Bz = B[Component::Z];

            auto const vyOnEx = GridLayout::template project<GridLayout::momentsToEx>(Vy, index);
            auto const vzOnEx = GridLayout::template project<GridLayout::momentsToEx>(Vz, index);
            auto const byOnEx = GridLayout::template project<GridLayout::ByToEx>(By, index);
            auto const bzOnEx = GridLayout::template project<GridLayout::BzToEx>(Bz, index);

            return -vyOnEx * bzOnEx + vzOnEx * byOnEx;
        }

        if constexpr (component == Component::Y)
        {
            auto const& Vx = Ve[Component::X];
            auto const& Vz = Ve[Component::Z];
            auto const& Bx = B[Component::X];
            auto const& Bz = B[Component::Z];

            auto const vxOnEy = GridLayout::template project<GridLayout::momentsToEy>(Vx, index);
            auto const vzOnEy = GridLayout::template project<GridLayout::momentsToEy>(Vz, index);
            auto const bxOnEy = GridLayout::template project<GridLayout::BxToEy>(Bx, index);
            auto const bzOnEy = GridLayout::template project<GridLayout::BzToEy>(Bz, index);

            return -vzOnEy * bxOnEy + vxOnEy * bzOnEy;
        }

        if constexpr (component == Component::Z)
        {
            auto const& Vx = Ve[Component::X];
            auto const& Vy = Ve[Component::Y];
            auto const& Bx = B[Component::X];
            auto const& By = B[Component::Y];

            auto const vxOnEz = GridLayout::template project<GridLayout::momentsToEz>(Vx, index);
            auto const vyOnEz = GridLayout::template project<GridLayout::momentsToEz>(Vy, index);
            auto const bxOnEz = GridLayout::template project<GridLayout::BxToEz>(Bx, index);
            auto const byOnEz = GridLayout::template project<GridLayout::ByToEz>(By, index);

            return -vxOnEz * byOnEz + vyOnEz * bxOnEz;
        }
    }


    template<auto component, typename Field, typename GradPe>
    PHARE_ALWAYS_INLINE auto pressure_(Field const& n, GradPe const& gradPe,
                                       MeshIndex<Field::dimension> index) const
    {
        if constexpr (component == Component::X)
        {
            auto const nOnEx = GridLayout::template project<GridLayout::momentsToEx>(n, index);

            return -gradPe(index) / nOnEx;
        }

        else if constexpr (component == Component::Y)
        {
            if constexpr (Field::dimension >= 2)
            {
                auto const nOnEy = GridLayout::template project<GridLayout::momentsToEy>(n, index);

                return -gradPe(index) / nOnEy;
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

                return -gradPe(index) / nOnEz;
            }
            else
            {
                return 0.;
            }
        }
    }




    template<auto component, typename VecField>
    PHARE_ALWAYS_INLINE auto resistive_(VecField const& J,
                                        MeshIndex<VecField::dimension> index) const
    {
        auto const& Jxyx = J[component];

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

    template<auto component, HyperMode hyperMode, typename VecField, typename Field>
    PHARE_ALWAYS_INLINE auto hyperresistive_(VecField const& J, VecField const& B, Field const& n,
                                             MeshIndex<VecField::dimension> index) const
    {
        if constexpr (hyperMode == HyperMode::constant)
            return constant_hyperresistive_<component>(J, index);
        else
            return spatial_hyperresistive_<component>(J, B, n, index);
    }


    template<auto component, typename VecField>
    PHARE_ALWAYS_INLINE auto constant_hyperresistive_(VecField const& J,
                                                      MeshIndex<VecField::dimension> index) const
    { // TODO : https://github.com/PHAREHUB/PHARE/issues/3
        return -nu * laplacian(layout_, J[component], index);
    }


    template<auto component, typename VecField, typename Field>
    PHARE_ALWAYS_INLINE auto spatial_hyperresistive_(VecField const& J, VecField const& B,
                                                     Field const& n,
                                                     MeshIndex<VecField::dimension> index) const
    {
        // 1 / 4^level without pow: a libcall the loop vectoriser cannot hoist (errno)
        auto const lvlCoeff
            = 1. / static_cast<double>(std::size_t{1} << (2 * layout_.levelNumber()));
        auto constexpr min_density = 0.1;
        auto computeHR
            = [&]<auto BxProj, auto ByProj, auto BzProj, auto nProj>() PHARE_LAMBDA_ALWAYS_INLINE {
                  auto const BxOnE = GridLayout::template project<BxProj>(B[Component::X], index);
                  auto const ByOnE = GridLayout::template project<ByProj>(B[Component::Y], index);
                  auto const BzOnE = GridLayout::template project<BzProj>(B[Component::Z], index);
                  auto const nOnE  = GridLayout::template project<nProj>(n, index);
                  auto b           = std::sqrt(BxOnE * BxOnE + ByOnE * ByOnE + BzOnE * BzOnE);
                  return -nu * (b / (nOnE + min_density) + 1) * lvlCoeff
                         * laplacian(layout_, J[component], index);
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
