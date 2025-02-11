#include "tests/core/data/mhd_state/init_functions.hpp"
#include "gtest/gtest.h"

#include "core/numerics/time_integrator/euler_integrator.hpp"
#include "core/numerics/time_integrator/tvdrk2_integrator.hpp"
#include "core/numerics/time_integrator/tvdrk3_integrator.hpp"

#include "core/numerics/reconstructions/constant.hpp"
#include "core/numerics/reconstructions/linear.hpp"
#include "core/numerics/reconstructions/weno3.hpp"

#include "core/numerics/slope_limiters/min_mod.hpp"
#include "core/numerics/slope_limiters/van_leer.hpp"

#include "core/numerics/riemann_solvers/rusanov.hpp"
#include "core/numerics/riemann_solvers/hll.hpp"

#include "core/numerics/MHD_equations/MHD_equations.hpp"

#include "tests/core/numerics/mock_mhd_simulator/test_mhd_solver.hpp"


template<auto dim>
PHARE::core::PHAREDict getDict()
{
    using initfunc = PHARE::core::InitFunction<dim>;

    PHARE::core::PHAREDict dict;
    dict["time_step"]  = 0.1;
    dict["final_time"] = 1.0;

    dict["mesh_size"]["x"] = 0.1;
    dict["nbr_cells"]["x"] = 100;
    dict["origin"]["x"]    = 0.0;


    if constexpr (dim > 1)
    {
        dict["mesh_size"]["y"] = 0.1;
        dict["nbr_cells"]["y"] = 100;
        dict["origin"]["y"]    = 0.0;

        if constexpr (dim > 2)
        {
            dict["mesh_size"]["z"] = 0.1;
            dict["nbr_cells"]["z"] = 100;
            dict["origin"]["z"]    = 0.0;
        }
    }

    dict["fv_method"]["resistivity"]               = 0.0;
    dict["fv_method"]["hyper_resistivity"]         = 0.0;
    dict["fv_method"]["heat_capacity_ratio"]       = 5.0 / 3.0;
    dict["fv_euler"]["heat_capacity_ratio"]        = 5.0 / 3.0;
    dict["to_primitive"]["heat_capacity_ratio"]    = 5.0 / 3.0;
    dict["to_conservative"]["heat_capacity_ratio"] = 5.0 / 3.0;

    dict["state"]["name"] = std::string("state");


    if constexpr (dim == 1)
    {
        using namespace PHARE::initializer::test_fn::func_1d;

        dict["state"]["density"]["initializer"] = static_cast<initfunc>(density);

        dict["state"]["velocity"]["initializer"]["x_component"] = static_cast<initfunc>(vx);
        dict["state"]["velocity"]["initializer"]["y_component"] = static_cast<initfunc>(vy);
        dict["state"]["velocity"]["initializer"]["z_component"] = static_cast<initfunc>(vz);

        dict["state"]["magnetic"]["initializer"]["x_component"] = static_cast<initfunc>(bx);
        dict["state"]["magnetic"]["initializer"]["y_component"] = static_cast<initfunc>(by);
        dict["state"]["magnetic"]["initializer"]["z_component"] = static_cast<initfunc>(bz);

        dict["state"]["pressure"]["initializer"] = static_cast<initfunc>(pressure);
    }

    if constexpr (dim == 2)
    {
        using namespace PHARE::initializer::test_fn::func_2d;

        dict["state"]["density"]["initializer"] = static_cast<initfunc>(density);

        dict["state"]["velocity"]["initializer"]["x_component"] = static_cast<initfunc>(vx);
        dict["state"]["velocity"]["initializer"]["y_component"] = static_cast<initfunc>(vy);
        dict["state"]["velocity"]["initializer"]["z_component"] = static_cast<initfunc>(vz);

        dict["state"]["magnetic"]["initializer"]["x_component"] = static_cast<initfunc>(bx);
        dict["state"]["magnetic"]["initializer"]["y_component"] = static_cast<initfunc>(by);
        dict["state"]["magnetic"]["initializer"]["z_component"] = static_cast<initfunc>(bz);

        dict["state"]["pressure"]["initializer"] = static_cast<initfunc>(pressure);
    }

    return dict;
}

TEST(TestSolver, Advance)
{
    using namespace PHARE::core;

    static bool constexpr hall              = true;
    static bool constexpr resistivity       = true;
    static bool constexpr hyper_resistivity = true;
    static auto constexpr dim               = 1;
    static auto constexpr ord               = 2;

    auto dict = getDict<dim>();
    MHDMockSimulator<dim, ord, TVDRK2Integrator, LinearReconstruction, VanLeerLimiter, HLL,
                     MHDEquations, hall, resistivity, hyper_resistivity>
        sim{dict};

    ASSERT_NO_THROW(sim.advance("test_mhd_state_1d.h5", 1));
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
