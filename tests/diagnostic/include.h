
#ifndef PHARE_TEST_DIAGNOSTIC_INCLUDE
#define PHARE_TEST_DIAGNOSTIC_INCLUDE

#include <map>
#include <memory>
#include <string>
#include <vector>
#include <iostream>

#include <SAMRAI/algs/TimeRefinementIntegrator.h>
#include <SAMRAI/algs/TimeRefinementLevelStrategy.h>
#include <SAMRAI/geom/CartesianGridGeometry.h>
#include <SAMRAI/hier/CoarsenOperator.h>
#include <SAMRAI/hier/RefineOperator.h>
#include <SAMRAI/mesh/GriddingAlgorithm.h>
#include <SAMRAI/mesh/StandardTagAndInitStrategy.h>
#include <SAMRAI/mesh/StandardTagAndInitialize.h>
#include <SAMRAI/tbox/InputManager.h>
#include <SAMRAI/tbox/Logger.h>
#include <SAMRAI/tbox/SAMRAIManager.h>
#include <SAMRAI/tbox/SAMRAI_MPI.h>
#include <SAMRAI/xfer/CoarsenAlgorithm.h>
#include <SAMRAI/xfer/RefineAlgorithm.h>
#include <SAMRAI/mesh/BergerRigoutsos.h>
#include <SAMRAI/mesh/TreeLoadBalancer.h>

#include "amr/types/amr_types.h"
#include "core/data/electromag/electromag.h"
#include "core/data/grid/gridlayout.h"
#include "core/data/grid/gridlayout_impl.h"
#include "core/data/grid/gridlayoutimplyee.h"
#include "core/data/ions/ion_population/ion_population.h"
#include "core/data/ions/ions.h"
#include "core/data/ions/particle_initializers/maxwellian_particle_initializer.h"
#include "core/data/ndarray/ndarray_vector.h"
#include "core/data/particles/particle_array.h"
#include "core/data/vecfield/vecfield.h"
#include "initializer/data_provider.h"
#include "core/hybrid/hybrid_quantities.h"
#include "solver/messenger_registration.h"
#include "amr/messengers/hybrid_messenger.h"
#include "amr/messengers/messenger_factory.h"
#include "amr/resources_manager/resources_manager.h"
#include "solver/physical_models/hybrid_model.h"
#include "solver/physical_models/mhd_model.h"
#include "solver/solvers/solver_mhd.h"
#include "solver/solvers/solver_ppc.h"

#endif /*PHARE_TEST_DIAGNOSTIC_INCLUDE*/
