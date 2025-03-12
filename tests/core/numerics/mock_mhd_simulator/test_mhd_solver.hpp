#ifndef PHARE_TESTS_CORE_NUMERICS_TEST_MHD_SOLVER_FIXTURES_HPP
#define PHARE_TESTS_CORE_NUMERICS_TEST_MHD_SOLVER_FIXTURES_HPP

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "amr/physical_models/physical_model.hpp"
#include "core/data/vecfield/vecfield.hpp"
#include "core/data/vecfield/vecfield_component.hpp"
#include "core/utilities/point/point.hpp"
#include "initializer/data_provider.hpp"

#include "highfive/H5File.hpp"
#include "highfive/H5DataSet.hpp"
#include "highfive/H5DataSpace.hpp"

#include "phare_core.hpp"
#include "amr/types/amr_types.hpp"
#include "amr/messengers/messenger.hpp"
#include "amr/wrappers/hierarchy.hpp"
#include "amr/solvers/solver_mhd.hpp"
#include "amr/solvers/solver.hpp"
#include "amr/solvers/solver_mhd_model_view.hpp"
#include "amr/physical_models/physical_model.hpp"

#include "core/data/vecfield/vecfield.hpp"
#include "core/data/vecfield/vecfield_component.hpp"
#include "core/utilities/point/point.hpp"
#include "core/mhd/mhd_quantities.hpp"
#include "core/data/grid/gridlayoutimplyee_mhd.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "amr/solvers/time_integrator/tvdrk2_integrator.hpp"
#include "amr/solvers/time_integrator/tvdrk3_integrator.hpp"
#include "core/numerics/godunov_fluxes/godunov_fluxes.hpp"

#include "tests/core/data/field/test_field_fixtures_mhd.hpp"
#include "tests/core/data/gridlayout/test_gridlayout.hpp"
#include "tests/core/data/mhd_state/test_mhd_state_fixtures.hpp"
#include "tests/core/data/field/test_usable_field_fixtures_mhd.hpp"
#include "tests/core/data/vecfield/test_vecfield_fixtures_mhd.hpp"

using MHDQuantity = PHARE::core::MHDQuantity;
using Direction   = PHARE::core::Direction;

struct DummyHierarchy
{
    auto getPatchLevel(std::size_t lvl) const
    {
        int a     = 1;
        int* ptrA = &a;
        return ptrA;
    };
};


struct DummyTypes
{
    using patch_t     = int;
    using level_t     = int;
    using hierarchy_t = DummyHierarchy;
};

struct DummyResourcesManager
{
    void registerResources(auto& resource) {}

    void allocate(auto& resource, auto& patch, double const allocateTime) const {}
};

template<std::size_t dim, std::size_t order>
struct DummyMHDModel : public PHARE::solver::IPhysicalModel<DummyTypes>
{
    static constexpr auto dimension = dim;
    using FieldMHD                  = PHARE::core::FieldMHD<dimension>;
    using VecFieldMHD               = PHARE::core::VecField<FieldMHD, MHDQuantity>;
    using YeeLayout_t               = PHARE::core::GridLayoutImplYeeMHD<dimension, order>;
    using GridLayout_t              = PHARE::core::GridLayout<YeeLayout_t>;

    using field_type                 = FieldMHD;
    using vecfield_type              = VecFieldMHD;
    using state_type                 = PHARE::core::MHDState<vecfield_type>;
    using gridlayout_type            = GridLayout_t;
    static constexpr auto model_name = "mhd_model";

    DummyMHDModel(GridLayout_t const& layout, PHARE::initializer::PHAREDict const& dict)
        : PHARE::solver::IPhysicalModel<DummyTypes>(model_name)
        , usablestate{layout, dict["state"]}
        , state{usablestate.super()}
        , resourcesManager{std::make_shared<DummyResourcesManager>()}
    {
        state.initialize(layout);
    }

    void initialize(level_t& level) override {}

    void allocate(patch_t& patch, double const allocateTime) override {}

    void fillMessengerInfo(std::unique_ptr<PHARE::amr::IMessengerInfo> const& info) const override
    {
    }

    PHARE::core::UsableMHDState<dimension> usablestate;
    PHARE::core::MHDState<VecFieldMHD>& state;
    std::shared_ptr<DummyResourcesManager> resourcesManager;
};

template<std::size_t dimension, std::size_t order>
struct DummyModelView : public PHARE::solver::ISolverModelView
{
    using YeeLayout_t  = PHARE::core::GridLayoutImplYeeMHD<dimension, order>;
    using GridLayout_t = PHARE::core::GridLayout<YeeLayout_t>;

    using FieldMHD    = PHARE::core::FieldMHD<dimension>;
    using VecFieldMHD = PHARE::core::VecField<FieldMHD, MHDQuantity>;

    DummyModelView(GridLayout_t& layout, DummyMHDModel<dimension, order>& model)
        : model_{model}
    {
        layouts.push_back(&layout);
    }

    auto& model() { return model_; }
    auto& model() const { return model_; }

    std::vector<GridLayout_t*> layouts;

    DummyMHDModel<dimension, order>& model_;
};

template<std::size_t dimension, std::size_t order>
class DummyMessenger : public PHARE::amr::IMessenger<PHARE::solver::IPhysicalModel<DummyTypes>>
{
    using YeeLayout_t    = PHARE::core::GridLayoutImplYeeMHD<dimension, order>;
    using GridLayout_t   = PHARE::core::GridLayout<YeeLayout_t>;
    using IPhysicalModel = PHARE::solver::IPhysicalModel<DummyTypes>;
    using level_t        = DummyTypes::level_t;

    std::string name() override { return "DummyMessenger"; }

    std::unique_ptr<PHARE::amr::IMessengerInfo> emptyInfoFromCoarser() override
    {
        return std::make_unique<PHARE::amr::IMessengerInfo>();
    }

    std::unique_ptr<PHARE::amr::IMessengerInfo> emptyInfoFromFiner() override
    {
        return std::make_unique<PHARE::amr::IMessengerInfo>();
    }

    void allocate(SAMRAI::hier::Patch& patch, double const allocateTime) const override {}

    void registerQuantities(std::unique_ptr<PHARE::amr::IMessengerInfo> fromCoarserInfo,
                            std::unique_ptr<PHARE::amr::IMessengerInfo> fromFinerInfo) override
    {
    }

    void registerLevel(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                       int level) override
    {
    }

    void regrid(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                int const levelNumber, std::shared_ptr<SAMRAI::hier::PatchLevel> const& oldLevel,
                IPhysicalModel& model, double const initDataTime) override
    {
    }

    void initLevel(IPhysicalModel& model, SAMRAI::hier::PatchLevel& level,
                   double const initDataTime) override
    {
    }

    void firstStep(IPhysicalModel& model, SAMRAI::hier::PatchLevel& level,
                   std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                   double const currentTime, double const prevCoarserTime,
                   double const newCoarserTime) override
    {
    }

    void lastStep(IPhysicalModel& model, SAMRAI::hier::PatchLevel& level) override {}

    void prepareStep(IPhysicalModel& model, SAMRAI::hier::PatchLevel& level,
                     double currentTime) override
    {
    }

    void fillRootGhosts(IPhysicalModel& model, SAMRAI::hier::PatchLevel& level,
                        double const initDataTime) override
    {
    }

    void synchronize(SAMRAI::hier::PatchLevel& level) override {}

    void postSynchronize(IPhysicalModel& model, SAMRAI::hier::PatchLevel& level,
                         double const time) override
    {
    }

    std::string fineModelName() const override { return "mhd_model"; }
    std::string coarseModelName() const override { return "mhd_model"; }

    template<typename Field>
    void fillGhosts_(Field& F)
    {
        auto const centering = layout_.centering(F.physicalQuantity());
        auto const isDualX   = static_cast<std::uint32_t>(centering[0]);

        auto const phStartX = layout_.physicalStartIndex(F, Direction::X);
        auto const phEndX   = layout_.physicalEndIndex(F, Direction::X);

        auto const nghost = layout_.nbrGhosts();

        if constexpr (dimension == 1)
        {
            for (auto g = 1u; g <= nghost; ++g)
            {
                F(phStartX - g) = F(phEndX - g + isDualX);   // left ghost
                F(phEndX + g)   = F(phStartX + g - isDualX); // right ghost
            }
        }
        else if constexpr (dimension >= 2)
        {
            auto const [ix0, ix1] = layout_.physicalStartToEnd(F, Direction::X);
            auto const [iy0, iy1] = layout_.physicalStartToEnd(F, Direction::Y);

            auto const isDualY = static_cast<std::uint32_t>(centering[1]);

            auto const phStartY = layout_.physicalStartIndex(F, Direction::Y);
            auto const phEndY   = layout_.physicalEndIndex(F, Direction::Y);

            if constexpr (dimension == 2)
            {
                for (auto g = 1u; g <= nghost; ++g)
                {
                    for (auto iy = iy0; iy <= iy1; ++iy)
                    {
                        F(phStartX - g, iy) = F(phEndX - g + isDualX,
                                                iy); // left ghost
                        F(phEndX + g, iy)   = F(phStartX + g - isDualX,
                                                iy); // right ghost
                    }
                    for (auto ix = ix0; ix <= ix1; ++ix)
                    {
                        F(ix, phStartY - g) = F(ix, phEndY - g + isDualY);   // bottom ghost
                        F(ix, phEndY + g)   = F(ix, phStartY + g - isDualY); // top ghost
                    }
                }
                // corners
                for (auto g1 = 1u; g1 <= nghost; ++g1)
                {
                    for (auto g2 = 1u; g2 <= nghost; ++g2)
                    {
                        F(phStartX - g1, phStartY - g2) = F(phEndX - g1 + isDualX,
                                                            phEndY - g2 + isDualY); // bottom left
                        F(phEndX + g1, phStartY - g2)
                            = F(phStartX + g1 - isDualX, phEndY - g2 + isDualY); // bottom right
                        F(phStartX - g1, phEndY + g2) = F(phEndX - g1 + isDualX,
                                                          phStartY + g2 - isDualY); // top left
                        F(phEndX + g1, phEndY + g2)   = F(phStartX + g1 - isDualX,
                                                          phStartY + g2 - isDualY); // top right
                    }
                }
            }
            else if constexpr (dimension == 3)
            {
                auto const [ix0, ix1] = layout_.physicalStartToEnd(F, Direction::X);
                auto const [iy0, iy1] = layout_.physicalStartToEnd(F, Direction::Y);
                auto const [iz0, iz1] = layout_.physicalStartToEnd(F, Direction::Z);

                auto const isDualZ = static_cast<std::uint32_t>(centering[2]);

                auto const phStartZ = layout_.physicalStartIndex(F, Direction::Z);
                auto const phEndZ   = layout_.physicalEndIndex(F, Direction::Z);

                for (auto g = 1u; g <= nghost; ++g)
                {
                    for (auto ix = ix0; ix <= ix1; ++ix)
                    {
                        for (auto iy = iy0; iy <= iy1; ++iy)
                        {
                            F(ix, iy, phStartZ - g) = F(ix, iy,
                                                        phEndZ - g + isDualZ); // bottom ghost
                            F(ix, iy, phEndZ + g)   = F(ix, iy,
                                                        phStartZ + g - isDualZ); // top ghost
                        }
                        for (auto iz = iz0; iz <= iz1; ++iz)
                        {
                            F(ix, phStartY - g, iz) = F(ix, phEndY - g + isDualY,
                                                        iz); // front ghost
                            F(ix, phEndY + g, iz)   = F(ix, phStartY + g - isDualY,
                                                        iz); // back ghost
                        }
                    }
                    for (auto iy = iy0; iy <= iy1; ++iy)
                    {
                        for (auto iz = iz0; iz <= iz1; ++iz)
                        {
                            F(phStartX - g, iy, iz) = F(phEndX - g + isDualX, iy,
                                                        iz); // left ghost
                            F(phEndX + g, iy, iz)   = F(phStartX + g - isDualX, iy,
                                                        iz); // right ghost
                        }
                    }
                }
                // corners
                for (auto g1 = 1u; g1 <= nghost; ++g1)
                {
                    for (auto g2 = 1u; g2 <= nghost; ++g2)
                    {
                        for (auto g3 = 1u; g3 <= nghost; ++g3)
                        {
                            F(phStartX - g1, phStartY - g2, phStartZ - g3)
                                = F(phEndX - g1 + isDualX, phEndY - g2 + isDualY,
                                    phEndZ - g3 + isDualZ); // left front bottom
                            F(phEndX + g1, phStartY - g2, phStartZ - g3)
                                = F(phStartX + g1 - isDualX, phEndY - g2 + isDualY,
                                    phEndZ - g3 + isDualZ); // right front bottom
                            F(phStartX - g1, phEndY + g2, phStartZ - g3)
                                = F(phEndX - g1 + isDualX, phStartY + g2 - isDualY,
                                    phEndZ - g3 + isDualZ); // left top front
                            F(phStartX - g1, phStartY - g2, phEndZ + g3)
                                = F(phEndX - g1 + isDualX, phEndY - g2 + isDualY,
                                    phStartZ + g3 - isDualZ); // left bottom back
                            F(phEndX + g1, phStartY - g2, phEndZ + g3)
                                = F(phStartX + g1 - isDualX, phEndY - g2 + isDualY,
                                    phStartZ + g3 - isDualZ); // right bottom back
                            F(phStartX - g1, phEndY + g2, phEndZ + g3)
                                = F(phEndX - g1 + isDualX, phStartY + g2 - isDualY,
                                    phStartZ + g3 - isDualZ); // left top back
                            F(phEndX + g1, phEndY + g2, phStartZ - g3)
                                = F(phStartX + g1 - isDualX, phStartY + g2 - isDualY,
                                    phEndZ - g3 + isDualZ); // right top front
                            F(phEndX + g1, phEndY + g2, phEndZ + g3)
                                = F(phStartX + g1 - isDualX, phStartY + g2 - isDualY,
                                    phStartZ + g3 - isDualZ); // right top back
                        }
                    }
                }
            }
        }
    }

    GridLayout_t layout_;

public:
    DummyMessenger(GridLayout_t const& layout)
        : layout_{layout}
    {
    }

    template<typename State>
    void fillMomentsGhosts(State& state, level_t& level, double const newTime)
    {
        fillGhosts_(state.rho);
        fillGhosts_(state.V(PHARE::core::Component::X));
        fillGhosts_(state.V(PHARE::core::Component::Y));
        fillGhosts_(state.V(PHARE::core::Component::Z));
        fillGhosts_(state.P);
    }

    template<typename Field>
    void fillMomentGhosts(Field& F, level_t& level, double const newTime)
    {
        fillGhosts_(F);
    }

    template<typename VecField>
    void fillMagneticGhosts(VecField& B, level_t& level, double const newTime)
    {
        auto& Bx = B(PHARE::core::Component::X);
        auto& By = B(PHARE::core::Component::Y);
        auto& Bz = B(PHARE::core::Component::Z);

        fillGhosts_(Bx);
        fillGhosts_(By);
        fillGhosts_(Bz);
    }

    template<typename VecField>
    void fillCurrentGhosts(VecField& J, level_t& level, double const newTime)
    {
        auto& Jx = J(PHARE::core::Component::X);
        auto& Jy = J(PHARE::core::Component::Y);
        auto& Jz = J(PHARE::core::Component::Z);

        fillGhosts_(Jx);
        fillGhosts_(Jy);
        fillGhosts_(Jz);
    }

    template<typename VecField>
    void fillElectricGhosts(VecField& E, level_t& level, double const newTime)
    {
        auto& Ex = E(PHARE::core::Component::X);
        auto& Ey = E(PHARE::core::Component::Y);
        auto& Ez = E(PHARE::core::Component::Z);

        fillGhosts_(Ex);
        fillGhosts_(Ey);
        fillGhosts_(Ez);
    }

    template<typename VecField>
    void fillMagneticFluxGhosts(VecField F_B, level_t& level, double const newTime)
    {
        auto& F_Bx = F_B(PHARE::core::Component::X);
        auto& F_By = F_B(PHARE::core::Component::Y);
        auto& F_Bz = F_B(PHARE::core::Component::Z);

        fillGhosts_(F_Bx);
        fillGhosts_(F_By);
        fillGhosts_(F_Bz);
    }
};

template<std::size_t dimension, std::size_t order, typename TimeIntegrator,
         template<typename> typename FVMethodStrategy, typename MHDModel>
class RessourceSetter
{
public:
    using YeeLayout_t  = PHARE::core::GridLayoutImplYeeMHD<dimension, order>;
    using GridLayout_t = PHARE::core::GridLayout<YeeLayout_t>;

    RessourceSetter(GridLayout_t const& layout)
        : rho_fx{"rho_fx", layout, MHDQuantity::Scalar::ScalarFlux_x}
        , rhoV_fx{"rhoV_fx", layout, MHDQuantity::Vector::VecFlux_x}
        , B_fx{"B_fx", layout, MHDQuantity::Vector::VecFlux_x}
        , Etot_fx{"Etot_fx", layout, MHDQuantity::Scalar::ScalarFlux_x}

        , rho_fy{"rho_fy", layout, MHDQuantity::Scalar::ScalarFlux_y}
        , rhoV_fy{"rhoV_fy", layout, MHDQuantity::Vector::VecFlux_y}
        , B_fy{"B_fy", layout, MHDQuantity::Vector::VecFlux_y}
        , Etot_fy{"Etot_fy", layout, MHDQuantity::Scalar::ScalarFlux_y}

        , rho_fz{"rho_fz", layout, MHDQuantity::Scalar::ScalarFlux_z}
        , rhoV_fz{"rhoV_fz", layout, MHDQuantity::Vector::VecFlux_z}
        , B_fz{"B_fz", layout, MHDQuantity::Vector::VecFlux_z}
        , Etot_fz{"Etot_fz", layout, MHDQuantity::Scalar::ScalarFlux_z}

        , usablestate1{layout, "state1"}
        , usablestate2{layout, "state2"}
    {
    }

    template<typename Solver>
    void operator()(Solver& solver)
    {
        auto&& [fluxes, evolve] = solver.getCompileTimeResourcesViewList();

        rho_fx.set_on(fluxes.rho_fx);
        rhoV_fx.set_on(fluxes.rhoV_fx);
        B_fx.set_on(fluxes.B_fx);
        Etot_fx.set_on(fluxes.Etot_fx);

        rho_fy.set_on(fluxes.rho_fy);
        rhoV_fy.set_on(fluxes.rhoV_fy);
        B_fy.set_on(fluxes.B_fy);
        Etot_fy.set_on(fluxes.Etot_fy);

        rho_fz.set_on(fluxes.rho_fz);
        rhoV_fz.set_on(fluxes.rhoV_fz);
        B_fz.set_on(fluxes.B_fz);
        Etot_fz.set_on(fluxes.Etot_fz);


        if constexpr (std::is_same_v<TimeIntegrator,
                                     PHARE::solver::TVDRK3Integrator<FVMethodStrategy, MHDModel>>)
        {
            auto&& [s1, s2] = evolve.getCompileTimeResourcesViewList();

            usablestate1.set_on(s1);
            usablestate2.set_on(s2);
        }
        else if constexpr (std::is_same_v<TimeIntegrator, PHARE::solver::TVDRK2Integrator<
                                                              FVMethodStrategy, MHDModel>>)
        {
            auto&& [s1] = evolve.getCompileTimeResourcesViewList();

            usablestate1.set_on(s1);
        }
    }

private:
    PHARE::core::UsableFieldMHD<dimension> rho_fx;
    PHARE::core::UsableVecFieldMHD<dimension> rhoV_fx;
    PHARE::core::UsableVecFieldMHD<dimension> B_fx;
    PHARE::core::UsableFieldMHD<dimension> Etot_fx;

    PHARE::core::UsableFieldMHD<dimension> rho_fy;
    PHARE::core::UsableVecFieldMHD<dimension> rhoV_fy;
    PHARE::core::UsableVecFieldMHD<dimension> B_fy;
    PHARE::core::UsableFieldMHD<dimension> Etot_fy;

    PHARE::core::UsableFieldMHD<dimension> rho_fz;
    PHARE::core::UsableVecFieldMHD<dimension> rhoV_fz;
    PHARE::core::UsableVecFieldMHD<dimension> B_fz;
    PHARE::core::UsableFieldMHD<dimension> Etot_fz;

    PHARE::core::UsableMHDState<dimension> usablestate1;
    PHARE::core::UsableMHDState<dimension> usablestate2;
};

class H5writer
{
public:
    template<typename Field, typename GridLayout>
    static void writeField(Field const& field, GridLayout const& layout, HighFive::File& h5file,
                           std::string const& groupName, std::string const& fieldName)
    {
        auto group
            = h5file.exist(groupName) ? h5file.getGroup(groupName) : h5file.createGroup(groupName);

        auto sizes = layout.allocSize(field.physicalQuantity());

        std::vector<std::size_t> h5_sizes(sizes.begin(), sizes.end());

        auto dataset = group.createDataSet<double>(fieldName, HighFive::DataSpace(h5_sizes));

        dataset.write_raw(field.data());
    }

    static std::string timeToGroupName(double time)
    {
        std::stringstream ss;
        ss << std::fixed << std::setprecision(6) << time;
        return ss.str();
    }
};

template<std::size_t dim, std::size_t ord,
         template<template<typename> typename, typename> typename TimeIntegrator,
         template<typename, typename> typename Reconstruction, typename SlopeLimiter,
         template<typename, bool> typename RiemannSolver,
         template<bool, bool, bool> typename Equations, bool Hall, bool Resistivity,
         bool HyperResistivity>
class MHDMockSimulator
{
public:
    static constexpr std::size_t dimension = dim;
    static constexpr std::size_t order     = ord;

    MHDMockSimulator(PHARE::initializer::PHAREDict const& dict)
        : dt_{dict["time_step"].template to<double>()}
        , final_time_{dict["final_time"].template to<double>()}
        , layout_{[&]() {
            auto initData = initLayout(dict);
            return GridLayout_t(initData.meshSize, initData.nbrCells, initData.origin);
        }()}
        , model_{layout_, dict}
        , dummy_view_{layout_, model_}
        , dummy_hierachy_()
        , dummy_messenger_{layout_}
        , TestMHDSolver_{dict}
        , ressource_setter_{layout_}
        , to_conservative{dict["to_conservative"]}
        , to_primitive{dict["to_primitive"]}
    {
        ressource_setter_(TestMHDSolver_);
    }

    void advance(std::string const& filename, int const dumpfrequency)
    {
        HighFive::File h5file(filename, HighFive::File::Overwrite);

        double time = 0.0;

        std::string initialGroup = H5writer::timeToGroupName(time);

        to_conservative(dummy_view_.layouts, dummy_view_.model().state);

        H5writer::writeField(dummy_view_.model().state.rho, layout_, h5file, initialGroup, "rho");
        H5writer::writeField(dummy_view_.model().state.V(PHARE::core::Component::X), layout_,
                             h5file, initialGroup, "vx");
        H5writer::writeField(dummy_view_.model().state.V(PHARE::core::Component::Y), layout_,
                             h5file, initialGroup, "vy");
        H5writer::writeField(dummy_view_.model().state.V(PHARE::core::Component::Z), layout_,
                             h5file, initialGroup, "vz");
        H5writer::writeField(dummy_view_.model().state.B(PHARE::core::Component::X), layout_,
                             h5file, initialGroup, "bx");
        H5writer::writeField(dummy_view_.model().state.B(PHARE::core::Component::Y), layout_,
                             h5file, initialGroup, "by");
        H5writer::writeField(dummy_view_.model().state.B(PHARE::core::Component::Z), layout_,
                             h5file, initialGroup, "bz");
        H5writer::writeField(dummy_view_.model().state.P, layout_, h5file, initialGroup, "p");

        H5writer::writeField(dummy_view_.model().state.rhoV(PHARE::core::Component::X), layout_,
                             h5file, initialGroup, "rhovx");
        H5writer::writeField(dummy_view_.model().state.rhoV(PHARE::core::Component::Y), layout_,
                             h5file, initialGroup, "rhovy");
        H5writer::writeField(dummy_view_.model().state.rhoV(PHARE::core::Component::Z), layout_,
                             h5file, initialGroup, "rhovz");
        H5writer::writeField(dummy_view_.model().state.Etot, layout_, h5file, initialGroup, "etot");

        H5writer::writeField(dummy_view_.model().state.J(PHARE::core::Component::X), layout_,
                             h5file, initialGroup, "jx");
        H5writer::writeField(dummy_view_.model().state.J(PHARE::core::Component::Y), layout_,
                             h5file, initialGroup, "jy");
        H5writer::writeField(dummy_view_.model().state.J(PHARE::core::Component::Z), layout_,
                             h5file, initialGroup, "jz");
        H5writer::writeField(dummy_view_.model().state.E(PHARE::core::Component::X), layout_,
                             h5file, initialGroup, "ex");
        H5writer::writeField(dummy_view_.model().state.E(PHARE::core::Component::Y), layout_,
                             h5file, initialGroup, "ey");
        H5writer::writeField(dummy_view_.model().state.E(PHARE::core::Component::Z), layout_,
                             h5file, initialGroup, "ez");

        int step = 1;

        while (time < final_time_)
        {
            TestMHDSolver_.advanceLevel(dummy_hierachy_, 1, dummy_view_, dummy_messenger_, time,
                                        time + dt_);
            time += dt_;
            std::cout << time << std::endl;

            if (step % dumpfrequency == 0 || time >= final_time_)
            {
                std::string currentGroup = H5writer::timeToGroupName(time);

                to_primitive(dummy_view_.layouts, dummy_view_.model().state);

                H5writer::writeField(dummy_view_.model().state.rho, layout_, h5file, currentGroup,
                                     "rho");
                H5writer::writeField(dummy_view_.model().state.V(PHARE::core::Component::X),
                                     layout_, h5file, currentGroup, "vx");
                H5writer::writeField(dummy_view_.model().state.V(PHARE::core::Component::Y),
                                     layout_, h5file, currentGroup, "vy");
                H5writer::writeField(dummy_view_.model().state.V(PHARE::core::Component::Z),
                                     layout_, h5file, currentGroup, "vz");
                H5writer::writeField(dummy_view_.model().state.B(PHARE::core::Component::X),
                                     layout_, h5file, currentGroup, "bx");
                H5writer::writeField(dummy_view_.model().state.B(PHARE::core::Component::Y),
                                     layout_, h5file, currentGroup, "by");
                H5writer::writeField(dummy_view_.model().state.B(PHARE::core::Component::Z),
                                     layout_, h5file, currentGroup, "bz");
                H5writer::writeField(dummy_view_.model().state.P, layout_, h5file, currentGroup,
                                     "p");

                H5writer::writeField(dummy_view_.model().state.rhoV(PHARE::core::Component::X),
                                     layout_, h5file, currentGroup, "rhovx");
                H5writer::writeField(dummy_view_.model().state.rhoV(PHARE::core::Component::Y),
                                     layout_, h5file, currentGroup, "rhovy");
                H5writer::writeField(dummy_view_.model().state.rhoV(PHARE::core::Component::Z),
                                     layout_, h5file, currentGroup, "rhovz");
                H5writer::writeField(dummy_view_.model().state.Etot, layout_, h5file, currentGroup,
                                     "etot");

                H5writer::writeField(dummy_view_.model().state.J(PHARE::core::Component::X),
                                     layout_, h5file, currentGroup, "jx");
                H5writer::writeField(dummy_view_.model().state.J(PHARE::core::Component::Y),
                                     layout_, h5file, currentGroup, "jy");
                H5writer::writeField(dummy_view_.model().state.J(PHARE::core::Component::Z),
                                     layout_, h5file, currentGroup, "jz");
                H5writer::writeField(dummy_view_.model().state.E(PHARE::core::Component::X),
                                     layout_, h5file, currentGroup, "ex");
                H5writer::writeField(dummy_view_.model().state.E(PHARE::core::Component::Y),
                                     layout_, h5file, currentGroup, "ey");
                H5writer::writeField(dummy_view_.model().state.E(PHARE::core::Component::Z),
                                     layout_, h5file, currentGroup, "ez");

                if constexpr (0)
                {
                    auto&& [_rho_fx, _rhoV_fx, _B_fx, _Etot_fx, _rho_fy, _rhoV_fy, _B_fy, _Etot_fy,
                            _rho_fz, _rhoV_fz, _B_fz, _Etot_fz, _evolve]
                        = TestMHDSolver_.getCompileTimeResourcesViewList();

                    H5writer::writeField(_rho_fx, layout_, h5file, currentGroup, "rho_x");
                    H5writer::writeField(_rhoV_fx(PHARE::core::Component::X), layout_, h5file,
                                         currentGroup, "rhovx_x");
                    H5writer::writeField(_rhoV_fx(PHARE::core::Component::Y), layout_, h5file,
                                         currentGroup, "rhovy_x");
                    H5writer::writeField(_rhoV_fx(PHARE::core::Component::Z), layout_, h5file,
                                         currentGroup, "rhovz_x");
                    H5writer::writeField(_B_fx(PHARE::core::Component::X), layout_, h5file,
                                         currentGroup, "bx_x");
                    H5writer::writeField(_B_fx(PHARE::core::Component::Y), layout_, h5file,
                                         currentGroup, "by_x");
                    H5writer::writeField(_B_fx(PHARE::core::Component::Z), layout_, h5file,
                                         currentGroup, "bz_x");
                    H5writer::writeField(_Etot_fx, layout_, h5file, currentGroup, "etot_x");

                    H5writer::writeField(_rho_fy, layout_, h5file, currentGroup, "rho_y");
                    H5writer::writeField(_rhoV_fy(PHARE::core::Component::X), layout_, h5file,
                                         currentGroup, "rhovx_y");
                    H5writer::writeField(_rhoV_fy(PHARE::core::Component::Y), layout_, h5file,
                                         currentGroup, "rhovy_y");
                    H5writer::writeField(_rhoV_fy(PHARE::core::Component::Z), layout_, h5file,
                                         currentGroup, "rhovz_y");
                    H5writer::writeField(_B_fy(PHARE::core::Component::X), layout_, h5file,
                                         currentGroup, "bx_y");
                    H5writer::writeField(_B_fy(PHARE::core::Component::Y), layout_, h5file,
                                         currentGroup, "by_y");
                    H5writer::writeField(_B_fy(PHARE::core::Component::Z), layout_, h5file,
                                         currentGroup, "bz_y");
                    H5writer::writeField(_Etot_fy, layout_, h5file, currentGroup, "etot_y");

                    auto&& [s1] = _evolve.getCompileTimeResourcesViewList();

                    H5writer::writeField(s1.rho, layout_, h5file, currentGroup, "rho1");
                    H5writer::writeField(s1.V(PHARE::core::Component::X), layout_, h5file,
                                         currentGroup, "vx1");
                    H5writer::writeField(s1.V(PHARE::core::Component::Y), layout_, h5file,
                                         currentGroup, "vy1");
                    H5writer::writeField(s1.V(PHARE::core::Component::Z), layout_, h5file,
                                         currentGroup, "vz1");
                    H5writer::writeField(s1.B(PHARE::core::Component::X), layout_, h5file,
                                         currentGroup, "bx1");
                    H5writer::writeField(s1.B(PHARE::core::Component::Y), layout_, h5file,
                                         currentGroup, "by1");
                    H5writer::writeField(s1.B(PHARE::core::Component::Z), layout_, h5file,
                                         currentGroup, "bz1");
                    H5writer::writeField(s1.P, layout_, h5file, currentGroup, "p1");

                    H5writer::writeField(s1.rhoV(PHARE::core::Component::X), layout_, h5file,
                                         currentGroup, "rhovx1");
                    H5writer::writeField(s1.rhoV(PHARE::core::Component::Y), layout_, h5file,
                                         currentGroup, "rhovy1");
                    H5writer::writeField(s1.rhoV(PHARE::core::Component::Z), layout_, h5file,
                                         currentGroup, "rhovz1");
                    H5writer::writeField(s1.Etot, layout_, h5file, currentGroup, "etot1");
                }

                to_conservative(dummy_view_.layouts, dummy_view_.model().state);
            }
            step++;
        }
    }

private:
    using YeeLayout_t  = PHARE::core::GridLayoutImplYeeMHD<dimension, order>;
    using GridLayout_t = PHARE::core::GridLayout<YeeLayout_t>;
    using Model_t      = DummyMHDModel<dimension, order>;
    using ModelView_t  = DummyModelView<dimension, order>;
    using Messenger_t  = DummyMessenger<dimension, order>;

    using Equations_t = Equations<Hall, Resistivity, HyperResistivity>;

    template<typename Layout>
    using RiemannSolver_t = RiemannSolver<Layout, Hall>;

    template<typename Layout>
    using Reconstruction_t = Reconstruction<Layout, SlopeLimiter>;

    template<typename Layout>
    using FVMethodStrategy
        = PHARE::core::Godunov<Layout, Reconstruction_t, RiemannSolver_t, Equations_t>;

    using TimeIntegrator_t = TimeIntegrator<FVMethodStrategy, Model_t>;

    using RessourceSetter_t
        = RessourceSetter<dimension, order, TimeIntegrator_t, FVMethodStrategy, Model_t>;

    using ToConservativeConverter_t
        = PHARE::solver::Dispatchers<GridLayout_t>::ToConservativeConverter_t;

    using ToPrimitiveConverter_t = PHARE::solver::Dispatchers<GridLayout_t>::ToPrimitiveConverter_t;

    double dt_;
    double final_time_;

    GridLayout_t layout_;
    Model_t model_;
    ModelView_t dummy_view_;
    DummyHierarchy dummy_hierachy_;
    Messenger_t dummy_messenger_;
    PHARE::solver::SolverMHD<Model_t, DummyTypes, TimeIntegrator_t, Messenger_t, ModelView_t>
        TestMHDSolver_;
    RessourceSetter_t ressource_setter_;

    ToConservativeConverter_t to_conservative;
    ToPrimitiveConverter_t to_primitive;

    struct LayoutInitData
    {
        std::array<double, dimension> meshSize;
        std::array<std::uint32_t, dimension> nbrCells;
        PHARE::core::Point<double, dimension> origin;
    };

    LayoutInitData initLayout(PHARE::initializer::PHAREDict const& dict) const
    {
        LayoutInitData data;

        if constexpr (dimension == 1)
        {
            data.meshSize = {dict["mesh_size"]["x"].template to<double>()};
            data.nbrCells = {static_cast<std::uint32_t>(dict["nbr_cells"]["x"].template to<int>())};
            data.origin   = {dict["origin"]["x"].template to<double>()};
        }
        else if constexpr (dimension == 2)
        {
            data.meshSize = {dict["mesh_size"]["x"].template to<double>(),
                             dict["mesh_size"]["y"].template to<double>()};
            data.nbrCells = {static_cast<std::uint32_t>(dict["nbr_cells"]["x"].template to<int>()),
                             static_cast<std::uint32_t>(dict["nbr_cells"]["y"].template to<int>())};
            data.origin   = {dict["origin"]["x"].template to<double>(),
                             dict["origin"]["y"].template to<double>()};
        }
        else if constexpr (dimension == 3)
        {
            data.meshSize = {dict["mesh_size"]["x"].template to<double>(),
                             dict["mesh_size"]["y"].template to<double>(),
                             dict["mesh_size"]["z"].template to<double>()};
            data.nbrCells = {static_cast<std::uint32_t>(dict["nbr_cells"]["x"].template to<int>()),
                             static_cast<std::uint32_t>(dict["nbr_cells"]["y"].template to<int>()),
                             static_cast<std::uint32_t>(dict["nbr_cells"]["z"].template to<int>())};
            data.origin   = {dict["origin"]["x"].template to<double>(),
                             dict["origin"]["y"].template to<double>(),
                             dict["origin"]["z"].template to<double>()};
        }

        return data;
    }
};

#endif // PHARE_TESTS_CORE_NUMERICS_TEST_MHD_SOLVER_FIXTURES_HPP
