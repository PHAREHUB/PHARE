#ifndef PHARE_CORE_NUMERICS_BASE_MHD_TIMESTEPPER_HPP
#define PHARE_CORE_NUMERICS_BASE_MHD_TIMESTEPPER_HPP

#include "initializer/data_provider.hpp"
#include "core/numerics/godunov_fluxes/godunov_utils.hpp"

namespace PHARE::solver
{
template<typename MHDModel>
class BaseMHDTimestepper
{
    using FieldT      = typename MHDModel::field_type;
    using VecFieldT   = typename MHDModel::vecfield_type;
    using GridLayoutT = typename MHDModel::gridlayout_type;

public:
    BaseMHDTimestepper(PHARE::initializer::PHAREDict const& dict)
        : butcherFluxes_{{"timeRho_fx", core::MHDQuantity::Scalar::ScalarFlux_x},
                         {"timeRhoV_fx", core::MHDQuantity::Vector::VecFlux_x},
                         {"timeB_fx", core::MHDQuantity::Vector::VecFlux_x},
                         {"timeEtot_fx", core::MHDQuantity::Scalar::ScalarFlux_x},

                         {"timeRho_fy", core::MHDQuantity::Scalar::ScalarFlux_y},
                         {"timeRhoV_fy", core::MHDQuantity::Vector::VecFlux_y},
                         {"timeB_fy", core::MHDQuantity::Vector::VecFlux_y},
                         {"timeEtot_fy", core::MHDQuantity::Scalar::ScalarFlux_y},

                         {"timeRho_fz", core::MHDQuantity::Scalar::ScalarFlux_z},
                         {"timeRhoV_fz", core::MHDQuantity::Vector::VecFlux_z},
                         {"timeB_fz", core::MHDQuantity::Vector::VecFlux_z},
                         {"timeEtot_fz", core::MHDQuantity::Scalar::ScalarFlux_z}}
        , butcherE_{"timeE", core::MHDQuantity::Vector::E}
        , butcherSources_{{"timeB1_source", core::MHDQuantity::Vector::B},
                          {"timeEtot_source", core::MHDQuantity::Scalar::Etot1}}
    {
    }

    void registerResources(MHDModel& model)
    {
        model.resourcesManager->registerResources(butcherFluxes_);
        model.resourcesManager->registerResources(butcherE_);
        model.resourcesManager->registerResources(butcherSources_);
    }

    void allocate(MHDModel& model, auto& patch, double const allocateTime) const
    {
        model.resourcesManager->allocate(butcherFluxes_, patch, allocateTime);
        model.resourcesManager->allocate(butcherE_, patch, allocateTime);
        model.resourcesManager->allocate(butcherSources_, patch, allocateTime);
    }

    void fillMessengerInfo(auto& info) const {}

    NO_DISCARD auto getCompileTimeResourcesViewList()
    {
        return std::forward_as_tuple(butcherFluxes_, butcherE_, butcherSources_);
    }

    NO_DISCARD auto getCompileTimeResourcesViewList() const
    {
        return std::forward_as_tuple(butcherFluxes_, butcherE_, butcherSources_);
    }

    auto exposeFluxes() { return std::forward_as_tuple(butcherFluxes_, butcherE_); }

    auto exposeFluxes() const { return std::forward_as_tuple(butcherFluxes_, butcherE_); }

    auto& exposeSources() { return butcherSources_; }

    auto const& exposeSources() const { return butcherSources_; }

protected:
    void resetButcherFluxes_(MHDModel& model, auto& level)
    {
        for (auto& patch : level)
        {
            auto const& layout = amr::layoutFromPatch<GridLayoutT>(*patch);
            auto _ = model.resourcesManager->setOnPatch(*patch, butcherFluxes_, butcherE_);

            evalFluxesOnGhostBox(
                layout, [&](auto& left, auto const&... args) mutable { left(args...) = 0.0; },
                butcherFluxes_);

            layout.evalOnGhostBox(butcherE_(core::Component::X), [&](auto const&... args) mutable {
                butcherE_(core::Component::X)(args...) = 0.0;
            });

            layout.evalOnGhostBox(butcherE_(core::Component::Y), [&](auto const&... args) mutable {
                butcherE_(core::Component::Y)(args...) = 0.0;
            });

            layout.evalOnGhostBox(butcherE_(core::Component::Z), [&](auto const&... args) mutable {
                butcherE_(core::Component::Z)(args...) = 0.0;
            });
        }
    }

    void accumulateButcherFluxes_(MHDModel& model, auto& E, auto& fluxes, auto& level,
                                  double const coef = 1.0)
    {
        for (auto& patch : level)
        {
            auto const& layout = amr::layoutFromPatch<GridLayoutT>(*patch);
            auto _
                = model.resourcesManager->setOnPatch(*patch, butcherFluxes_, butcherE_, fluxes, E);

            evalFluxesOnGhostBox(
                layout,
                [&](auto& left, auto const& right, auto const&... args) mutable {
                    left(args...) += right(args...) * coef;
                },
                butcherFluxes_, fluxes);


            layout.evalOnGhostBox(butcherE_(core::Component::X), [&](auto const&... args) mutable {
                butcherE_(core::Component::X)(args...) += E(core::Component::X)(args...) * coef;
            });

            layout.evalOnGhostBox(butcherE_(core::Component::Y), [&](auto const&... args) mutable {
                butcherE_(core::Component::Y)(args...) += E(core::Component::Y)(args...) * coef;
            });

            layout.evalOnGhostBox(butcherE_(core::Component::Z), [&](auto const&... args) mutable {
                butcherE_(core::Component::Z)(args...) += E(core::Component::Z)(args...) * coef;
            });
        }
    }

    // Body-source Butcher accumulation, parallel to fluxes/E but over the DOMAIN box (evalOnBox):
    // sources are volumetric, not fluxes through faces, so they take no part in the coarse-fine
    // flux correction. reset then accumulate the per-stage `sources` with the same RK weight used
    // for the fluxes; the final Butcher-flux Euler applies butcherSources_.
    void resetButcherSources_(MHDModel& model, auto& level)
    {
        for (auto& patch : level)
        {
            auto const& layout = amr::layoutFromPatch<GridLayoutT>(*patch);
            auto _             = model.resourcesManager->setOnPatch(*patch, butcherSources_);

            for (auto const& c : {core::Component::X, core::Component::Y, core::Component::Z})
            {
                auto& src = butcherSources_.B1_source(c);
                layout.evalOnBox(src, [&](auto const&... args) mutable { src(args...) = 0.0; });
            }
            layout.evalOnBox(butcherSources_.Etot_source, [&](auto const&... args) mutable {
                butcherSources_.Etot_source(args...) = 0.0;
            });
        }
    }

    void accumulateButcherSources_(MHDModel& model, auto& sources, auto& level,
                                   double const coef = 1.0)
    {
        for (auto& patch : level)
        {
            auto const& layout = amr::layoutFromPatch<GridLayoutT>(*patch);
            auto _ = model.resourcesManager->setOnPatch(*patch, butcherSources_, sources);

            for (auto const& c : {core::Component::X, core::Component::Y, core::Component::Z})
            {
                auto& dst       = butcherSources_.B1_source(c);
                auto const& src = sources.B1_source(c);
                layout.evalOnBox(
                    dst, [&](auto const&... args) mutable { dst(args...) += src(args...) * coef; });
            }
            layout.evalOnBox(butcherSources_.Etot_source, [&](auto const&... args) mutable {
                butcherSources_.Etot_source(args...) += sources.Etot_source(args...) * coef;
            });
        }
    }

    core::AllFluxes<FieldT, VecFieldT> butcherFluxes_;
    VecFieldT butcherE_;
    core::MHDSources<FieldT, VecFieldT> butcherSources_;
};
} // namespace PHARE::solver

#endif
