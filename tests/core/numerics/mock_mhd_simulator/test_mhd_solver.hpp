#ifndef PHARE_TESTS_CORE_NUMERICS_TEST_MHD_SOLVER_FIXTURES_HPP
#define PHARE_TESTS_CORE_NUMERICS_TEST_MHD_SOLVER_FIXTURES_HPP

#include <amr/physical_models/physical_model.hpp>
#include <core/data/vecfield/vecfield.hpp>
#include <core/data/vecfield/vecfield_component.hpp>
#include <core/numerics/godunov_fluxes/godunov_fluxes.hpp>
#include <core/utilities/point/point.hpp>
#include <cstddef>
#include <cstdint>
#include <initializer/data_provider.hpp>
#include <string>
#include <vector>

#include "highfive/H5File.hpp"
#include "highfive/H5DataSet.hpp"
#include "highfive/H5DataSpace.hpp"

#include "phare_core.hpp"
#include "amr/types/amr_types.hpp"
#include "amr/messengers/messenger.hpp"
#include "amr/wrappers/hierarchy.hpp"
#include "amr/solvers/solver_mhd.hpp"
#include "amr/solvers/solver.hpp"

#include "core/mhd/mhd_quantities.hpp"
#include "core/data/grid/gridlayoutimplyee_mhd.hpp"
#include "core/numerics/boundary_condition/boundary_condition.hpp"

#include "tests/core/data/field/test_field_fixtures_mhd.hpp"
#include "tests/core/data/gridlayout/test_gridlayout.hpp"
#include "tests/core/data/mhd_state/test_mhd_state_fixtures.hpp"
#include "tests/core/data/field/test_usable_field_fixtures_mhd.hpp"
#include "tests/core/data/vecfield/test_vecfield_fixtures_mhd.hpp"

using MHDQuantity = PHARE::core::MHDQuantity;
using Direction   = PHARE::core::Direction;


template<std::size_t dimension, std::size_t order>
struct DummyModelViewConstructor
{
    using YeeLayout_t     = PHARE::core::GridLayoutImplYeeMHD<dimension, order>;
    using GridLayout_t    = PHARE::core::GridLayout<YeeLayout_t>;
    using GodunovFluxes_t = PHARE::core::GodunovFluxes<GridLayout_t>;

    DummyModelViewConstructor(GridLayout_t const& layout)
        : rho_x{"rho_x", layout, MHDQuantity::Scalar::ScalarFlux_x}
        , rhoV_x{"rhoV_x", layout, MHDQuantity::Vector::VecFlux_x}
        , B_x{"B_x", layout, MHDQuantity::Vector::VecFlux_x}
        , Etot_x{"Etot_x", layout, MHDQuantity::Scalar::ScalarFlux_x}

        , rho_y{"rho_y", layout, MHDQuantity::Scalar::ScalarFlux_y}
        , rhoV_y{"rhoV_y", layout, MHDQuantity::Vector::VecFlux_y}
        , B_y{"B_y", layout, MHDQuantity::Vector::VecFlux_y}
        , Etot_y{"Etot_y", layout, MHDQuantity::Scalar::ScalarFlux_y}

        , rho_z{"rho_z", layout, MHDQuantity::Scalar::ScalarFlux_z}
        , rhoV_z{"rhoV_z", layout, MHDQuantity::Vector::VecFlux_z}
        , B_z{"B_z", layout, MHDQuantity::Vector::VecFlux_z}
        , Etot_z{"Etot_z", layout, MHDQuantity::Scalar::ScalarFlux_z}

        , layouts{layout}
    {
    }

    PHARE::core::UsableFieldMHD<dimension> rho_x;
    PHARE::core::UsableVecFieldMHD<dimension> rhoV_x;
    PHARE::core::UsableVecFieldMHD<dimension> B_x;
    PHARE::core::UsableFieldMHD<dimension> Etot_x;

    PHARE::core::UsableFieldMHD<dimension> rho_y;
    PHARE::core::UsableVecFieldMHD<dimension> rhoV_y;
    PHARE::core::UsableVecFieldMHD<dimension> B_y;
    PHARE::core::UsableFieldMHD<dimension> Etot_y;

    PHARE::core::UsableFieldMHD<dimension> rho_z;
    PHARE::core::UsableVecFieldMHD<dimension> rhoV_z;
    PHARE::core::UsableVecFieldMHD<dimension> B_z;
    PHARE::core::UsableFieldMHD<dimension> Etot_z;

    GridLayout_t layouts;
};


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
    using patch_t     = PHARE::amr::SAMRAI_Types::patch_t;
    using level_t     = int;
    using hierarchy_t = DummyHierarchy;
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
    using gridlayout_type            = GridLayout_t;
    static constexpr auto model_name = "mhd_model";

    DummyMHDModel(GridLayout_t const& layout, PHARE::initializer::PHAREDict const& dict)
        : PHARE::solver::IPhysicalModel<DummyTypes>(model_name)
        , usablestate{layout, dict["state"]}
        , usablestate1{layout, dict["state1"]}
        , usablestate2{layout, dict["state2"]}
        , state{usablestate.super()}
        , state1{usablestate1.super()}
        , state2{usablestate2.super()}
    {
        state.initialize(layout);
    }

    void initialize(level_t& level) override {}

    void allocate(patch_t& patch, double const allocateTime) override {}

    void fillMessengerInfo(std::unique_ptr<PHARE::amr::IMessengerInfo> const& info) const override
    {
    }

    PHARE::core::UsableMHDState<dimension> usablestate;
    PHARE::core::UsableMHDState<dimension> usablestate1;
    PHARE::core::UsableMHDState<dimension> usablestate2;
    PHARE::core::MHDState<VecFieldMHD>& state;
    PHARE::core::MHDState<VecFieldMHD>& state1;
    PHARE::core::MHDState<VecFieldMHD>& state2;
};

template<std::size_t dimension, std::size_t order>
struct DummyModelView : public PHARE::solver::ISolverModelView
{
    using YeeLayout_t     = PHARE::core::GridLayoutImplYeeMHD<dimension, order>;
    using GridLayout_t    = PHARE::core::GridLayout<YeeLayout_t>;
    using GodunovFluxes_t = PHARE::core::GodunovFluxes<GridLayout_t>;

    using FieldMHD    = PHARE::core::FieldMHD<dimension>;
    using VecFieldMHD = PHARE::core::VecField<FieldMHD, MHDQuantity>;

    DummyModelView(DummyModelViewConstructor<dimension, order>& construct,
                   DummyMHDModel<dimension, order>& _model)
        : model_{_model}
    {
        rho.push_back(&model_.state.rho);
        V.push_back(&model_.state.V);
        B.push_back(&model_.state.B);
        P.push_back(&model_.state.P);
        J.push_back(&model_.state.J);
        E.push_back(&model_.state.E);
        rhoV.push_back(&model_.state.rhoV);
        Etot.push_back(&model_.state.Etot);

        rho_x.push_back(&construct.rho_x.super());
        rhoV_x.push_back(&construct.rhoV_x.super());
        B_x.push_back(&construct.B_x.super());
        Etot_x.push_back(&construct.Etot_x.super());

        rho_y.push_back(&construct.rho_y.super());
        rhoV_y.push_back(&construct.rhoV_y.super());
        B_y.push_back(&construct.B_y.super());
        Etot_y.push_back(&construct.Etot_y.super());

        rho_z.push_back(&construct.rho_z.super());
        rhoV_z.push_back(&construct.rhoV_z.super());
        B_z.push_back(&construct.B_z.super());
        Etot_z.push_back(&construct.Etot_z.super());

        rho1.push_back(&model_.state1.rho);
        V1.push_back(&model_.state1.V);
        B1.push_back(&model_.state1.B);
        P1.push_back(&model_.state1.P);
        J1.push_back(&model_.state1.J);
        E1.push_back(&model_.state1.E);
        rhoV1.push_back(&model_.state1.rhoV);
        Etot1.push_back(&model_.state1.Etot);

        rho2.push_back(&model_.state2.rho);
        V2.push_back(&model_.state2.V);
        B2.push_back(&model_.state2.B);
        P2.push_back(&model_.state2.P);
        J2.push_back(&model_.state2.J);
        E2.push_back(&model_.state2.E);
        rhoV2.push_back(&model_.state2.rhoV);
        Etot2.push_back(&model_.state2.Etot);

        layouts.push_back(&construct.layouts);
    }

    auto& model() { return model_; }
    auto& model() const { return model_; }

    std::vector<FieldMHD*> rho;
    std::vector<VecFieldMHD*> V;
    std::vector<VecFieldMHD*> B;
    std::vector<FieldMHD*> P;
    std::vector<VecFieldMHD*> J;
    std::vector<VecFieldMHD*> E;
    std::vector<VecFieldMHD*> rhoV;
    std::vector<FieldMHD*> Etot;

    std::vector<FieldMHD*> rho_x;
    std::vector<VecFieldMHD*> rhoV_x;
    std::vector<VecFieldMHD*> B_x;
    std::vector<FieldMHD*> Etot_x;

    std::vector<FieldMHD*> rho_y;
    std::vector<VecFieldMHD*> rhoV_y;
    std::vector<VecFieldMHD*> B_y;
    std::vector<FieldMHD*> Etot_y;

    std::vector<FieldMHD*> rho_z;
    std::vector<VecFieldMHD*> rhoV_z;
    std::vector<VecFieldMHD*> B_z;
    std::vector<FieldMHD*> Etot_z;

    std::vector<FieldMHD*> rho1;
    std::vector<VecFieldMHD*> V1;
    std::vector<VecFieldMHD*> B1;
    std::vector<FieldMHD*> P1;
    std::vector<VecFieldMHD*> J1;
    std::vector<VecFieldMHD*> E1;
    std::vector<VecFieldMHD*> rhoV1;
    std::vector<FieldMHD*> Etot1;

    std::vector<FieldMHD*> rho2;
    std::vector<VecFieldMHD*> V2;
    std::vector<VecFieldMHD*> B2;
    std::vector<FieldMHD*> P2;
    std::vector<VecFieldMHD*> J2;
    std::vector<VecFieldMHD*> E2;
    std::vector<VecFieldMHD*> rhoV2;
    std::vector<FieldMHD*> Etot2;

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

    void registerLevel(const std::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
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
                   const std::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
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
        auto centering = layout_.centering(F.physicalQuantity());
        auto isDualX   = static_cast<std::uint32_t>(centering[0]);

        auto phStartX = layout_.physicalStartIndex(F, Direction::X);
        auto phEndX   = layout_.physicalEndIndex(F, Direction::X);

        auto nghost = layout_.nbrGhosts();

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

            auto isDualY = static_cast<std::uint32_t>(centering[1]);

            auto phStartY = layout_.physicalStartIndex(F, Direction::Y);
            auto phEndY   = layout_.physicalEndIndex(F, Direction::Y);

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

                auto isDualZ = static_cast<std::uint32_t>(centering[2]);

                auto phStartZ = layout_.physicalStartIndex(F, Direction::Z);
                auto phEndZ   = layout_.physicalEndIndex(F, Direction::Z);

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
    void fillMagneticFluxGhosts(std::vector<VecField*> viewB, level_t& level, double const newTime)
    {
        for (auto& ptrB : viewB)
        {
            auto& B = *ptrB;

            auto& Bx = B(PHARE::core::Component::X);
            auto& By = B(PHARE::core::Component::Y);
            auto& Bz = B(PHARE::core::Component::Z);

            fillGhosts_(Bx);
            fillGhosts_(By);
            fillGhosts_(Bz);
        }
    }
};

class H5writer
{
public:
    template<typename Field, typename GridLayout>
    static void writeField(Field const& field, GridLayout const& layout, HighFive::File& h5file,
                           const std::string& groupName, const std::string& fieldName)
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

template<std::size_t dim, std::size_t ord>
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
        , dummy_view_construct_(layout_)
        , dummy_view_{dummy_view_construct_, model_}
        , dummy_hierachy_()
        , dummy_messenger_{layout_}
        , TestMHDSolver_{dict}
    {
    }

    void advance(const std::string& filename, const int dumpfrequency)
    {
        HighFive::File h5file(filename, HighFive::File::Overwrite);

        double time = 0.0;

        std::string initialGroup = H5writer::timeToGroupName(time);

        H5writer::writeField(model_.state.rho, layout_, h5file, initialGroup, "rho");
        H5writer::writeField(model_.state.V(PHARE::core::Component::X), layout_, h5file,
                             initialGroup, "vx");
        H5writer::writeField(model_.state.V(PHARE::core::Component::Y), layout_, h5file,
                             initialGroup, "vy");
        H5writer::writeField(model_.state.V(PHARE::core::Component::Z), layout_, h5file,
                             initialGroup, "vz");
        H5writer::writeField(model_.state.B(PHARE::core::Component::X), layout_, h5file,
                             initialGroup, "bx");
        H5writer::writeField(model_.state.B(PHARE::core::Component::Y), layout_, h5file,
                             initialGroup, "by");
        H5writer::writeField(model_.state.B(PHARE::core::Component::Z), layout_, h5file,
                             initialGroup, "bz");
        H5writer::writeField(model_.state.P, layout_, h5file, initialGroup, "p");

        H5writer::writeField(model_.state.rhoV(PHARE::core::Component::X), layout_, h5file,
                             initialGroup, "rhovx");
        H5writer::writeField(model_.state.rhoV(PHARE::core::Component::Y), layout_, h5file,
                             initialGroup, "rhovy");
        H5writer::writeField(model_.state.rhoV(PHARE::core::Component::Z), layout_, h5file,
                             initialGroup, "rhovz");
        H5writer::writeField(model_.state.Etot, layout_, h5file, initialGroup, "etot");

        H5writer::writeField(model_.state.J(PHARE::core::Component::X), layout_, h5file,
                             initialGroup, "jx");
        H5writer::writeField(model_.state.J(PHARE::core::Component::Y), layout_, h5file,
                             initialGroup, "jy");
        H5writer::writeField(model_.state.J(PHARE::core::Component::Z), layout_, h5file,
                             initialGroup, "jz");
        H5writer::writeField(model_.state.E(PHARE::core::Component::X), layout_, h5file,
                             initialGroup, "ex");
        H5writer::writeField(model_.state.E(PHARE::core::Component::Y), layout_, h5file,
                             initialGroup, "ey");
        H5writer::writeField(model_.state.E(PHARE::core::Component::Z), layout_, h5file,
                             initialGroup, "ez");

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

                H5writer::writeField(model_.state.rho, layout_, h5file, currentGroup, "rho");
                H5writer::writeField(model_.state.V(PHARE::core::Component::X), layout_, h5file,
                                     currentGroup, "vx");
                H5writer::writeField(model_.state.V(PHARE::core::Component::Y), layout_, h5file,
                                     currentGroup, "vy");
                H5writer::writeField(model_.state.V(PHARE::core::Component::Z), layout_, h5file,
                                     currentGroup, "vz");
                H5writer::writeField(model_.state.B(PHARE::core::Component::X), layout_, h5file,
                                     currentGroup, "bx");
                H5writer::writeField(model_.state.B(PHARE::core::Component::Y), layout_, h5file,
                                     currentGroup, "by");
                H5writer::writeField(model_.state.B(PHARE::core::Component::Z), layout_, h5file,
                                     currentGroup, "bz");
                H5writer::writeField(model_.state.P, layout_, h5file, currentGroup, "p");

                H5writer::writeField(model_.state.rhoV(PHARE::core::Component::X), layout_, h5file,
                                     currentGroup, "rhovx");
                H5writer::writeField(model_.state.rhoV(PHARE::core::Component::Y), layout_, h5file,
                                     currentGroup, "rhovy");
                H5writer::writeField(model_.state.rhoV(PHARE::core::Component::Z), layout_, h5file,
                                     currentGroup, "rhovz");
                H5writer::writeField(model_.state.Etot, layout_, h5file, currentGroup, "etot");

                H5writer::writeField(model_.state.J(PHARE::core::Component::X), layout_, h5file,
                                     currentGroup, "jx");
                H5writer::writeField(model_.state.J(PHARE::core::Component::Y), layout_, h5file,
                                     currentGroup, "jy");
                H5writer::writeField(model_.state.J(PHARE::core::Component::Z), layout_, h5file,
                                     currentGroup, "jz");
                H5writer::writeField(model_.state.E(PHARE::core::Component::X), layout_, h5file,
                                     currentGroup, "ex");
                H5writer::writeField(model_.state.E(PHARE::core::Component::Y), layout_, h5file,
                                     currentGroup, "ey");
                H5writer::writeField(model_.state.E(PHARE::core::Component::Z), layout_, h5file,
                                     currentGroup, "ez");
            }
            step++;
        }
    }

private:
    using YeeLayout_t            = PHARE::core::GridLayoutImplYeeMHD<dimension, order>;
    using GridLayout_t           = PHARE::core::GridLayout<YeeLayout_t>;
    using Model_t                = DummyMHDModel<dimension, order>;
    using ModelViewConstructor_t = DummyModelViewConstructor<dimension, order>;
    using ModelView_t            = DummyModelView<dimension, order>;
    using Messenger_t            = DummyMessenger<dimension, order>;

    double dt_;
    double final_time_;

    GridLayout_t layout_;
    Model_t model_;
    ModelViewConstructor_t dummy_view_construct_;
    ModelView_t dummy_view_;
    DummyHierarchy dummy_hierachy_;
    Messenger_t dummy_messenger_;
    PHARE::solver::SolverMHD<Model_t, DummyTypes, Messenger_t, ModelView_t> TestMHDSolver_;

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
