# AMR package

The PHARE package amr is the part of the code that depends on the specific AMR library used, here [SAMRAI](https://github.com/LLNL/SAMRAI).
The code in the amr package aims at two objectives:

- **Specializing SAMRAI types** used to manipulate data on the patch hierarchy, such as PatchData, Variable, BoxGeometry, Overlap, refine/time/coarsen operators, etc. This allows SAMRAI to manipulate our data containers implemented in the [core](md_core.html) package. This part demands a high expertise on the SAMRAI library. The main components here are:

    - for Fields:
        - [FieldData](@ref PHARE::FieldData)
        - [FieldDataFactory] (@ref PHARE::FieldDataFactory)
        - [FieldOverlap](@ref PHARE::FieldOverlap)
        - [FieldVariable](@ref PHARE::FieldVariable)
        - [FieldRefineOperator](@ref PHARE::FieldRefineOperator)
        - [FieldCoarsenOperator](@ref PHARE::FieldCoarsenOperator)
        - [TimeInterpolateOperator](@ref PHARE::FieldLinearTimeInterpolate)
    - for Particles
        - [ParticlesData](@ref PHARE::ParticlesData)
        - [ParticlesVariable](@ref PHARE::ParticlesVariable)
        - [ParticlesDataFactory](@ref PHARE::ParticlesDataFactory)
        - [ParticlesRefineOperator](@ref PHARE::ParticlesRefineOperator)



- **Providing the developer with a high level interface** for manipulating PHARE data types on a patch hierarchy, and communicating data between different refinement levels, without having to manipulate low-level SAMRAI objects. The main components are:

    - [ResourcesManager](@ref PHARE::ResourcesManager) is an interface for using PHARE data containers defined on patches without having to manipulate PatchData IDs, the variable database, etc.
    - [PhysicalState](@ref PHARE::PhysicalState) : is an abstraction of an object encapsulating a collection of PHARE data types that all together define the state of a model. For instance an MHDState would contain the density, bulk velocity, thermal pressure and magnetic field variables.
    - [IPhysicalModel](@ref PHARE::IPhysicalModel) contains a PhysicalState and a ResourcesManager. The PHARE variables in the PhysicalState are registered to the model ResourcesManager. This Model ResourceManager will also be used to register any variable (from Solver and Messenger) that are associated with this PhysicalModel. For instance a SolverPPC or a HybridMessenger will allocate, register and access their data in the HybridModel ResourcesManager.
    - [MultiPhysicsIntegrator](@ref PHARE::MultiPhysicsIntegrator) is high level object that inherits from the SAMRAI::algs::TimeRefinementLevelStrategy and the SAMRAI::mesh::StandardTagAndInitStrategy. Called from the SAMRAI system, this object has pointers to [IPhysicalModel](@ref IPhysicalModel), [ISolver](@ref ISolver) and [IMessenger](@ref IMessenger) to use at any specific level in the hierarchy, but does not know their concrete types. A [MultiPhysicsIntegrator](@ref PHARE::MultiPhysicsIntegrator) does know know what "MHD", or "Hybrid" means. It delegates level initializations to Initializers and Messengers, and solving a step to the Solvers.
    - [ISolver](@ref PHARE::ISolver) is a class that represents an abstract solver used by the [MultiPhysicsIntegrator](@ref PHARE::MultiPhysicsIntegrator) to advance a model from one time to the next.
    - [IMessenger](@ref PHARE::IMessenger) is an interface used by the [MultiPhysicsIntegrator](@ref PHARE::MultiPhysicsIntegrator) and Solvers to perform all actions requiring to communicate data between refinement levels and within a level.



