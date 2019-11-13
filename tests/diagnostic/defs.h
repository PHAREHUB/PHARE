
#ifndef PHARE_TEST_DIAGNOSTIC_DEF
#define PHARE_TEST_DIAGNOSTIC_DEF

#include "core/data/grid/gridlayout.h"
#include "core/data/grid/gridlayoutimplyee.h"
#include "core/data/vecfield/vecfield.h"
#include "core/data/ions/ion_population/ion_population.h"
#include "core/data/ions/ions.h"
#include "core/data/electromag/electromag.h"
#include "solver/physical_models/hybrid_model.h"
#include "amr/messengers/hybrid_messenger.h"
#include "amr/messengers/hybrid_hybrid_messenger_strategy.h"
#include "amr/resources_manager/resources_manager.h"

#include "func.h"


using namespace PHARE::amr;
using namespace PHARE::core;
using namespace PHARE::solver;

namespace PHARE_test
{
namespace _1d
{
    static constexpr std::size_t dim_         = 1;
    static constexpr std::size_t interpOrder_ = 1;

    using GridImplYeeT  = GridLayoutImplYee<dim_, interpOrder_>;
    using GridYeeT      = GridLayout<GridImplYeeT>;
    using FieldT        = Field<NdArrayVector1D<>, HybridQuantity::Scalar>;
    using VecFieldT     = VecField<NdArrayVector1D<>, HybridQuantity>;
    using IonsPopT      = IonPopulation<ParticleArray<dim_>, VecFieldT, GridYeeT>;
    using IonsT         = Ions<IonsPopT, GridYeeT>;
    using ElectromagT   = Electromag<VecFieldT>;
    using HybridModelT  = HybridModel<GridYeeT, ElectromagT, IonsT, SAMRAI_Types>;
    using HybridHybridT = HybridHybridMessengerStrategy<HybridModelT, IPhysicalModel<SAMRAI_Types>>;
    using HybridMessengerT  = HybridMessenger<HybridModelT, IPhysicalModel<SAMRAI_Types>>;
    using ResourcesManagerT = ResourcesManager<GridYeeT>;

    using ScalarFunctionT = PHARE::initializer::ScalarFunction<1>;

    PHARE::initializer::PHAREDict CREATE_IONS_DICT()
    {
        PHARE::initializer::PHAREDict dict;
        dict["ions"]["name"]                                 = std::string{"ions"};
        dict["ions"]["nbrPopulations"]                       = int{2};
        dict["ions"]["pop0"]["name"]                         = std::string{"protons"};
        dict["ions"]["pop0"]["mass"]                         = 1.;
        dict["ions"]["pop0"]["particle_initializer"]["name"] = std::string{"maxwellian"};
        dict["ions"]["pop0"]["particle_initializer"]["density"]
            = static_cast<ScalarFunctionT>(density);

        dict["ions"]["pop0"]["particle_initializer"]["bulk_velocity_x"]
            = static_cast<ScalarFunctionT>(vx);

        dict["ions"]["pop0"]["particle_initializer"]["bulk_velocity_y"]
            = static_cast<ScalarFunctionT>(vy);

        dict["ions"]["pop0"]["particle_initializer"]["bulk_velocity_z"]
            = static_cast<ScalarFunctionT>(vz);


        dict["ions"]["pop0"]["particle_initializer"]["thermal_velocity_x"]
            = static_cast<ScalarFunctionT>(vthx);

        dict["ions"]["pop0"]["particle_initializer"]["thermal_velocity_y"]
            = static_cast<ScalarFunctionT>(vthy);

        dict["ions"]["pop0"]["particle_initializer"]["thermal_velocity_z"]
            = static_cast<ScalarFunctionT>(vthz);


        dict["ions"]["pop0"]["particle_initializer"]["nbr_part_per_cell"] = int{100};
        dict["ions"]["pop0"]["particle_initializer"]["charge"]            = -1.;
        dict["ions"]["pop0"]["particle_initializer"]["basis"] = std::string{"cartesian"};

        dict["ions"]["pop1"]["name"]                         = std::string{"alpha"};
        dict["ions"]["pop1"]["mass"]                         = 1.;
        dict["ions"]["pop1"]["particle_initializer"]["name"] = std::string{"maxwellian"};
        dict["ions"]["pop1"]["particle_initializer"]["density"]
            = static_cast<ScalarFunctionT>(density);

        dict["ions"]["pop1"]["particle_initializer"]["bulk_velocity_x"]
            = static_cast<ScalarFunctionT>(vx);

        dict["ions"]["pop1"]["particle_initializer"]["bulk_velocity_y"]
            = static_cast<ScalarFunctionT>(vy);

        dict["ions"]["pop1"]["particle_initializer"]["bulk_velocity_z"]
            = static_cast<ScalarFunctionT>(vz);


        dict["ions"]["pop1"]["particle_initializer"]["thermal_velocity_x"]
            = static_cast<ScalarFunctionT>(vthx);

        dict["ions"]["pop1"]["particle_initializer"]["thermal_velocity_y"]
            = static_cast<ScalarFunctionT>(vthy);

        dict["ions"]["pop1"]["particle_initializer"]["thermal_velocity_z"]
            = static_cast<ScalarFunctionT>(vthz);


        dict["ions"]["pop1"]["particle_initializer"]["nbr_part_per_cell"] = int{100};
        dict["ions"]["pop1"]["particle_initializer"]["charge"]            = -1.;
        dict["ions"]["pop1"]["particle_initializer"]["basis"] = std::string{"cartesian"};

        dict["electromag"]["name"]             = std::string{"EM"};
        dict["electromag"]["electric"]["name"] = std::string{"E"};
        dict["electromag"]["magnetic"]["name"] = std::string{"B"};

        dict["electromag"]["electric"]["initializer"]["x_component"]
            = static_cast<ScalarFunctionT>(ex);
        dict["electromag"]["electric"]["initializer"]["y_component"]
            = static_cast<ScalarFunctionT>(ey);
        dict["electromag"]["electric"]["initializer"]["z_component"]
            = static_cast<ScalarFunctionT>(ez);

        dict["electromag"]["magnetic"]["initializer"]["x_component"]
            = static_cast<ScalarFunctionT>(bx);
        dict["electromag"]["magnetic"]["initializer"]["y_component"]
            = static_cast<ScalarFunctionT>(by);
        dict["electromag"]["magnetic"]["initializer"]["z_component"]
            = static_cast<ScalarFunctionT>(bz);

        return dict;
    }

} // namespace _1d
} // namespace PHARE_test

#endif /*PHARE_TEST_DIAGNOSTIC_DEF*/
