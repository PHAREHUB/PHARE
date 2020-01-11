
#ifndef PHARE_SIMULATOR_INCLUDE_H
#define PHARE_SIMULATOR_INCLUDE_H

#include "initializer/python_data_provider.h"
// intended blank HAVE_SYS_TIMES_His defined by samrai
#include "amr/types/amr_types.h"

#include "core/models/physical_state.h"
#include "core/data/electromag/electromag.h"
#include "core/data/grid/gridlayout.h"
#include "core/data/grid/gridlayoutimplyee.h"
#include "core/data/ions/ion_population/ion_population.h"
#include "core/data/ions/ions.h"
#include "core/data/ions/particle_initializers/maxwellian_particle_initializer.h"
#include "core/data/ndarray/ndarray_vector.h"
#include "core/data/particles/particle_array.h"
#include "core/data/vecfield/vecfield.h"
#include "core/models/physical_state.h"
#include "core/utilities/meta/meta_utilities.h"
#include "core/utilities/algorithm.h"

#include "initializer/data_provider.h"

#include "amr/messengers/messenger_factory.h"

#include "solver/level_initializer/level_initializer.h"
#include "solver/level_initializer/level_initializer_factory.h"
#include "solver/multiphysics_integrator.h"
#include "solver/physical_models/hybrid_model.h"
#include "solver/physical_models/mhd_model.h"
#include "solver/physical_models/physical_model.h"
#include "solver/solvers/solver.h"
#include "solver/solvers/solver_mhd.h"
#include "solver/solvers/solver_ppc.h"

#include <SAMRAI/algs/TimeRefinementIntegrator.h>
#include <SAMRAI/geom/CartesianGridGeometry.h>
#include <SAMRAI/hier/Box.h>
#include <SAMRAI/hier/BoxContainer.h>
#include <SAMRAI/hier/IntVector.h>
#include <SAMRAI/hier/PatchHierarchy.h>
#include <SAMRAI/mesh/BergerRigoutsos.h>
#include <SAMRAI/mesh/GriddingAlgorithm.h>
#include <SAMRAI/mesh/StandardTagAndInitialize.h>
#include <SAMRAI/mesh/TreeLoadBalancer.h>
#include <SAMRAI/tbox/Database.h>
#include <SAMRAI/tbox/DatabaseBox.h>
#include <SAMRAI/tbox/InputManager.h>
#include <SAMRAI/tbox/MemoryDatabase.h>

#include <algorithm>
#include <iterator>
#include <memory>
#include <sstream>

#endif /*PHARE_SIMULATOR_INCLUDE_H*/
