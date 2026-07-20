#ifndef PHARE_CORE_NUMERICS_GODUNOV_FLUXES_HPP
#define PHARE_CORE_NUMERICS_GODUNOV_FLUXES_HPP

#include "core/numerics/ohm/ohm.hpp"
#include "core/utilities/types.hpp"
#include "core/utilities/meta/meta_utilities.hpp"
#include "core/utilities/index/index.hpp"
#include "core/utilities/point/point.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/data/vecfield/vecfield_component.hpp"
#include "core/numerics/godunov_fluxes/godunov_utils.hpp"
#include "core/numerics/reconstructions/reconstructor.hpp"

#include "initializer/data_provider.hpp"

#include <cmath>
#include <tuple>
#include <cstddef>
#include <cstdint>
#include <variant>

namespace PHARE::core
{
template<size_t dim>
constexpr auto getDirections()
{
    if constexpr (dim == 1)
    {
        return std::make_tuple(Direction::X);
    }
    else if constexpr (dim == 2)
    {
        return std::make_tuple(Direction::X, Direction::Y);
    }
    else if constexpr (dim == 3)
    {
        return std::make_tuple(Direction::X, Direction::Y, Direction::Z);
    }
}

template<auto direction, size_t dim>
auto getGrow(int const nghosts, bool const hyperResistivity)
{
    Point<std::uint32_t, dim> p{};

    auto dir = static_cast<size_t>(direction);

    for (size_t i = 0; i < dim; ++i)
        if (i != dir)
            p[i] = nghosts;

    // add one extra layer in the direction of the flux laplacian computation. Maybe some later
    // optimisation would let us just compute for uct and have the extra layer only reconstructed
    // for j
    if (hyperResistivity)
        p[dir] += 1;

    return p;
}

struct GodunovInfo : public OhmInfo
{
    double const gamma;

    GodunovInfo static FROM(initializer::PHAREDict const& dict)
    {
        return {{OhmInfo::FROM(dict)}, dict["heat_capacity_ratio"].template to<double>()};
    }
};


template<typename GridLayout, template<typename> typename Reconstruction, typename RiemannSolver,
         typename Equations>
class Godunov : public GodunovInfo
{
    using Super                     = GodunovInfo;
    using Reconstruction_t          = Reconstruction<GridLayout>;
    using Reconstructor_t           = Reconstructor<Reconstruction_t>;
    using RiemannSolver_t           = RiemannSolver;
    constexpr static auto dimension = GridLayout::dimension;

public:
    using Info_t      = Super;
    using Equations_t = Equations;

    template<typename T>
    using Rec = Reconstruction<T>;

    constexpr static auto Hall = Equations::hall;

    explicit Godunov(GodunovInfo const& info, GridLayout const& layout)
        : Super{info}
        , layout_{layout}
        , is_resistive_{info.isResistive()}
        , is_hyper_resistive_{info.isHyperResistive()}
        , equations_{gamma}
        , riemann_{gamma}
    {
    }

    template<typename GState, typename State, typename Fluxes>
    void operator()(GState& fvm_state, auto& ct_state, State& state, Fluxes& fluxes)
    {
        constexpr auto directions = getDirections<dimension>();

        constexpr auto num_directions = std::tuple_size_v<std::decay_t<decltype(directions)>>;

        for_N<num_directions>([&](auto i) {
            constexpr Direction direction = std::get<i>(directions);

            auto fillHyperbolicFluxes = [&](auto mustSaveBtT) {
                constexpr bool mustSaveBt = decltype(mustSaveBtT)::value;
                layout_.evalOnBiggerBox(
                    fluxes.template expose_centering<direction>(),
                    getGrow<direction, dimension>(Reconstruction_t::nghosts, is_hyper_resistive_),
                    [&](auto&... indices) {
                        if constexpr (Hall)
                        {
                            auto&& [uL, uR] = Reconstructor_t::template reconstruct<direction>(
                                state, {indices...});

                            auto const& [jL, jR] = Reconstructor_t::template center_reconstruct<
                                direction, GridLayout::edgeXToCellCenter,
                                GridLayout::edgeYToCellCenter, GridLayout::edgeZToCellCenter>(
                                state.J, {indices...});

                            auto&& u      = std::forward_as_tuple(uL, uR);
                            auto const& j = std::forward_as_tuple(jL, jR);


                            auto const& [fL, fR] = for_N<2, for_N_R_mode::make_tuple>([&](auto i) {
                                return equations_.template compute<direction>(std::get<i>(u),
                                                                              std::get<i>(j));
                            });

                            fluxes.template get_dir<direction>({indices...})
                                = riemann_.template solve<direction>(uL, uR, fL, fR, jL, jR);

                            ct_state.template save<direction>(riemann_.vt, riemann_.jt,
                                                              riemann_.rhot, riemann_.uct_coefs,
                                                              {indices...});

                            // for energy ExB term
                            if constexpr (mustSaveBt)
                                save_tranverse_magnetic_field_<direction>(fvm_state, uL, uR,
                                                                          {indices...});
                        }
                        else // if (!Hall)
                        {
                            auto&& [uL, uR] = Reconstructor_t::template reconstruct<direction>(
                                state, {indices...});

                            auto&& u = std::forward_as_tuple(uL, uR);

                            auto const& [fL, fR] = for_N<2, for_N_R_mode::make_tuple>([&](auto i) {
                                return equations_.template compute<direction>(std::get<i>(u));
                            });

                            fluxes.template get_dir<direction>({indices...})
                                = riemann_.template solve<direction>(uL, uR, fL, fR);

                            if constexpr (mustSaveBt)
                            {
                                auto const& [jL, jR] = Reconstructor_t::template center_reconstruct<
                                    direction, GridLayout::edgeXToCellCenter,
                                    GridLayout::edgeYToCellCenter, GridLayout::edgeZToCellCenter>(
                                    state.J, {indices...});

                                auto const jt   = riemann_.vector_riemann_averaging(jL, jR);
                                auto const rhot = riemann_.riemann_averaging(uL.rho, uR.rho);

                                ct_state.template save<direction>(riemann_.vt, jt, rhot,
                                                                  riemann_.uct_coefs, {indices...});

                                // for energy ExB term
                                save_tranverse_magnetic_field_<direction>(fvm_state, uL, uR,
                                                                          {indices...});
                            }
                            else // ideal mhd
                                ct_state.template save<direction>(riemann_.vt, riemann_.uct_coefs,
                                                                  {indices...});
                        }
                    });
            };

            // adding resistive contributions to energy taking advantage of the already computed jt
            // fluxes for the laplacian computation. This probably doesn't need the grow as the
            // required quantities for ct are already saved.
            auto addResistiveContributions = [&](auto doResistiveT, auto doHyperT, auto hyperT) {
                constexpr bool doResistive = decltype(doResistiveT)::value;
                constexpr bool doHyper     = decltype(doHyperT)::value;
                constexpr HyperMode hyper  = decltype(hyperT)::value;
                layout_.evalOnBox(
                    fluxes.template expose_centering<direction>(), [&](auto&... indices) {
                        auto& Jt     = ct_state.template getJt<direction>();
                        auto& Bt     = getBt_<direction>(fvm_state);
                        auto F       = fluxes.template get_dir<direction>({indices...});
                        auto& F_B    = F.B;
                        auto& F_Etot = F.Etot();

                        auto const& Btidx = toPerIndexVector(Bt, {indices...});

                        if constexpr (doResistive)
                        {
                            // transverse B field components (probably a riemann operation).
                            auto const& Jtidx = toPerIndexVector(Jt, {indices...});
                            equations_.template resistive_contributions<direction>(
                                eta, Btidx, Jtidx, F_B, F_Etot);
                        }
                        if constexpr (doHyper)
                        {
                            auto const vecLaplJ
                                = transverse_laplacian_<direction>(Jt, {indices...});

                            if constexpr (hyper == HyperMode::constant)
                                return constant_hyperresistive_<direction>(Btidx, vecLaplJ, F_B,
                                                                           F_Etot);
                            else // HyperMode::spatial
                            {
                                auto const& Bn = toPerIndexVector(state.B, {indices...});
                                auto const& rhot
                                    = ct_state.template getRhot<direction>()(indices...);

                                return spatial_hyperresistive_<direction>(Btidx, Bn, vecLaplJ, rhot,
                                                                          F_B, F_Etot);
                            }
                        }
                    });
            };

            // loop at compile time to avoid runtime checks in compute loops
            std::visit(
                [&](auto isResistiveTag, auto isHyperResistiveTag, auto hyperModeTag) {
                    constexpr bool isResistive      = decltype(isResistiveTag)::value;
                    constexpr bool isHyperResistive = decltype(isHyperResistiveTag)::value;

                    fillHyperbolicFluxes(std::bool_constant < isResistive || isHyperResistive > {});

                    if constexpr (isResistive || isHyperResistive)
                        addResistiveContributions(isResistiveTag, isHyperResistiveTag,
                                                  hyperModeTag);
                },
                asBoolConstant(is_resistive_), asBoolConstant(is_hyper_resistive_),
                asEnumConstant(hyper_mode));
        });
    }

    bool isResistive() const { return is_resistive_; }
    bool isHyperResistive() const { return is_hyper_resistive_; }

private:
    template<auto direction>
    auto save_tranverse_magnetic_field_(auto& fvm_state, auto const& uL, auto const& uR,
                                        MeshIndex<dimension> idx)
    {
        auto Bidx = riemann_.vector_riemann_averaging(uL.B, uR.B);

        auto& Bt = getBt_<direction>(fvm_state);

        Bt(Component::X)(idx) = Bidx.x;
        Bt(Component::Y)(idx) = Bidx.y;
        Bt(Component::Z)(idx) = Bidx.z;
    }

    template<auto direction>
    auto& getBt_(auto& fvm_state) const
    {
        return fvm_state.template getBt<direction>();
    }

    template<auto direction>
    void constant_hyperresistive_(auto const& Bt, auto const& vecLaplJ, auto& F_B,
                                  auto& F_Etot) const
    {
        equations_.template resistive_contributions<direction>(-nu, Bt, vecLaplJ, F_B, F_Etot);
    }

    template<auto direction>
    void spatial_hyperresistive_(auto const& Bt, auto const& B, auto const& vecLaplJ,
                                 auto const& rhot, auto& F_B, auto& F_Etot) const
    {
        auto minMeshSize = [&]() {
            auto const meshSize = layout_.meshSize();
            if constexpr (dimension == 1)
                return meshSize[0];
            else if constexpr (dimension == 2)
                return std::min({meshSize[0], meshSize[1]});
            else
                return std::min({meshSize[0], meshSize[1], meshSize[2]});
        }();


        auto computeHR = [&](auto Bx, auto By, auto Bz) {
            auto b          = std::sqrt(Bx * Bx + By * By + Bz * Bz);
            auto const coef = -nu * minMeshSize * minMeshSize * (b / rhot + 1);
            equations_.template resistive_contributions<direction>(coef, Bt, vecLaplJ, F_B, F_Etot);
        };

        if constexpr (direction == Direction::X)
        {
            auto const Bx = B.x; // normal component
            auto const By = Bt.y;
            auto const Bz = Bt.z;

            computeHR(Bx, By, Bz);
        }
        else if constexpr (direction == Direction::Y)
        {
            auto const Bx = Bt.x;
            auto const By = B.y; // normal component
            auto const Bz = Bt.z;

            computeHR(Bx, By, Bz);
        }
        else if constexpr (direction == Direction::Z)
        {
            auto const Bx = Bt.x;
            auto const By = Bt.y;
            auto const Bz = B.z; // normal component

            computeHR(Bx, By, Bz);
        }
    }

    template<auto direction>
    auto transverse_laplacian_(auto const& Jt, MeshIndex<dimension> index) const
    {
        if constexpr (direction == Direction::X)
        {
            auto const JyLapl = layout_.laplacian(Jt(Component::Y), index);
            auto const JzLapl = layout_.laplacian(Jt(Component::Z), index);
            return PerIndexVector<double>{std::nan(""), JyLapl, JzLapl};
        }
        else if constexpr (direction == Direction::Y)
        {
            auto const JxLapl = layout_.laplacian(Jt(Component::X), index);
            auto const JzLapl = layout_.laplacian(Jt(Component::Z), index);
            return PerIndexVector<double>{JxLapl, std::nan(""), JzLapl};
        }
        else if constexpr (direction == Direction::Z)
        {
            auto const JxLapl = layout_.laplacian(Jt(Component::X), index);
            auto const JyLapl = layout_.laplacian(Jt(Component::Y), index);
            return PerIndexVector<double>{JxLapl, JyLapl, std::nan("")};
        }
    }


    GridLayout layout_;
    bool const is_resistive_;
    bool const is_hyper_resistive_;
    Equations equations_;
    RiemannSolver_t riemann_;
};

} // namespace PHARE::core

#endif
