# PHARE



## Initialization

### Running PHARE
The executable is launched like this :

```
./phare ini_file.[extension]
```

where ‘.extension’ is ‘.py’ or any other initializer format that PHARE knows. The prefered format is .py, that is a python script that defines all necessary parameters and functions to initialize a hybrid simulation.


### Root level initialization

In PHARE the  root level is a hybrid model. The initialization step consists in loading particles and electromagnetic fields on it. In SAMRAI, this step is done by the `StandardTagAndInitialize` object, which delegates the job to its abstract `StandardTagAndInitStrategy`. In PHARE, this strategy is derived in `MultiPhysicsIntegrator` (that also derives from `TimeRefinementLevelStrategy`). The `MultiPhysicsIntegrator` knows which physical model exists at each level, and in particular the one at the root level. In PHARE, the `MultiPhysicsIntegrator` always says the root level is a HybridModel. When SAMRAI calls the virtual function `initializeLevelData()`, and when the level number equals the root number, it will
initialize the hybrid model, i.e. load particles in Ions IonsPopulation and load Electromagnetic fields.


### Initializing Ions

#### Construction
The Ions constructor is called by the HybridModel constructor, when the `MultiPhysicsIntegrator` is constructed.
The Ions constructor takes an IonInitializer from which it takes, for each of its IonPopulation:

- the mass of an ion of that population
- the charge of an ion of that population
- the name of that population
- an abstract ParticleInitializer of that population

and give them to each IonPopulation. This information will be used by IonsPopulation when initializing the level data.



#### Loading particles
The loading of the particles is not done at the construction of the Ions object. Indeed, the Ions is just
an empty shell that offers a convenient interface to manipulate particles and moments, but it does not own
the resources. Therefore, loading of particles can only be done while looping over the patch hierarchy. This occurs
in `MultiPhysicsIntegrator::initializeLevelData()`.

```
for level in patchHierarchy:
  for patch in level:    
    resourcesGuard(ions, patch):
    layout = layoutFromPatch(patch);
    ions.loadParticles(layout);
```

In the above loop, the ions has its resource pointers set to the resources by the resourcesGuard. From that moment
it can access memory and thus load particles.

The particle loading depends on the `ParticleInitializer` each population received. Ions::loadParticles()
basically calls `IonsPopulation::loadParticles()` which itself calls `ParticleInitializer::loadParticles()`.
This last methods knows how to load particles in a particle array, provided it is given a GridLayout.
The layout associated with the patch is obtained from layoutFromPatch().

The ParticleInitializer mush load particles only in the physical domain of the patch. These particles will be loaded
in the `interiorParticles` particle array of the IonPopulation.

When the particle initializer is finished, the IonPopulation will not have:

- its `ghostParticles` loaded.
- its `coarseToFine_n` and `coarseToFine_np1` loaded.

`ghostParticles` are particles in the ghost particle region that extends over a neighbor patch
of the same level. These particles should *not* be loaded by the particle initializer since they are clones
of the particles that exist in the first cells inside the neighbor patch. These particles are obtained by applying
a `refineSchedule()` that communicates particles from the interior of patches of the same level to the ghost region of a patch.

`coarseToFine_n` and `coarseToFine_np1` are particles buffers that contain the particles that are split from underlying
coarse particles, at t_n and t_{n+1} of the coarse patch, respectively. These buffers will be filled by the coarser level
at the beginning and end of its time step, so the ParticleInitializer should *not* touch them.


#### Moments initialization

There is no need to do something specific for moments.



#### FluidParticleInitializer
There are several ways to load particle distributions. The most usual method is to load particles in locally
Maxwellian distributions which moments (n,V,T) are prescribed.



### Initializing Electromag
Electromagnetic fields used in PHARE are, like Ions, contained in the HybridModel object. Like the Ions, they will be
built with an initializer, called ElectromagInitializer. A loop on the hierarchy will give access to the patch  

TODO...



### Python initialization
PHARE uses pybind11 to directly run the python functions that return the value of density, magnetic field etc. necessary to initialize th root level.



TODO...



## Advancing the solution

Advancing a solution in SAMRAI is done by the `TimeRefinementIntegrator` with the refined time stepping option. The `TimeRefinementIntegrator` actually calls an abstract strategy `TimeRefinementLevelStrategy` to advance the solution on each level. In PHARE the `TimeRefinementLevelStrategy` is the `MultiPhysicsIntegrator`.


### The MultiPhysicsIntegrator

The `MultiPhysicsIntegrator` knows:

- which `PhysicalModel` is solved at each level
- which `Transaction` is used between levels
- which `Solver` is used at each level


### advanceLevel()

The `TimeRefinementIntegrator` will call the the method `MultiPhysicsIntegrator::avdanceLevel()` at each time step to advance the solution.
The function is given the following arguments:

- the level
- the patch hierarchy
- current_time : the time at the beginning of the step
- new_time : the time after the step to perform
- first_step : a bool that is true if the step to be done is the first of the substeps
- last_step : a bool that is true if the step to be done is the last step before synchronizing with next coarser level
- regrid_advance : a bool that is true only if the advanceLevel() is called in a reggriding procedure (richardson extrapolation), in PHARE it is always false.

When it is called, the `MultiPhysicsIntegrator` must :

- check which model is used for the given level
- which Transaction to use.
- which solver to used for the given level

once the solver is determined, its function Solver::advanceLevel() is called and given the Transaction and Model to use. In PHARE the model will always be a `HybridModel`, and the Transaction a `HybridTransaction`.



### Transactions



### HybridSolverPPC::advanceLevel()

The `HybridSolverPPC` is a `Solver` that solves hybrid equations using a `HybridModel` to access data.










