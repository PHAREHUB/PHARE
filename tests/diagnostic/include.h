
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
#include <SAMRAI/mesh/ChopAndPackLoadBalancer.h>
#include <SAMRAI/mesh/GriddingAlgorithm.h>
#include <SAMRAI/mesh/StandardTagAndInitStrategy.h>
#include <SAMRAI/mesh/StandardTagAndInitialize.h>
#include <SAMRAI/mesh/TileClustering.h>
#include <SAMRAI/tbox/InputManager.h>
#include <SAMRAI/tbox/Logger.h>
#include <SAMRAI/tbox/SAMRAIManager.h>
#include <SAMRAI/tbox/SAMRAI_MPI.h>
#include <SAMRAI/xfer/CoarsenAlgorithm.h>
#include <SAMRAI/xfer/RefineAlgorithm.h>

<<<<<<< HEAD
=======
#include "gmock/gmock.h"
#include "gtest/gtest.h"

>>>>>>> Enable Hi5 MPI if HDF5_IS_PARALLEL
#include "types/amr_types.h"
#include "data/electromag/electromag.h"
#include "data/grid/gridlayout.h"
#include "data/grid/gridlayout_impl.h"
#include "data/grid/gridlayoutimplyee.h"
#include "data/ions/ion_population/ion_population.h"
#include "data/ions/ions.h"
#include "data/ions/particle_initializers/maxwellian_particle_initializer.h"
#include "data/ndarray/ndarray_vector.h"
#include "data/particles/particle_array.h"
#include "data/vecfield/vecfield.h"
#include "data_provider.h"
#include "hybrid/hybrid_quantities.h"
#include "messenger_registration.h"
#include "messengers/hybrid_messenger.h"
#include "messengers/messenger_factory.h"
#include "physical_models/hybrid_model.h"
#include "physical_models/mhd_model.h"
#include "resources_manager/resources_manager.h"
#include "solvers/solver_mhd.h"
#include "solvers/solver_ppc.h"

#endif /*PHARE_TEST_DIAGNOSTIC_INCLUDE*/
