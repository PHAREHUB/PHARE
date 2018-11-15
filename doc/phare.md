# PHARE Developer documentation


PHARE is a Parallel Hybrid code with Adaptive Mesh REfinement.

## Projet components

The PHARE project is divided into:

- [core](core.md) : package which gathers all the source code completely independent of the AMR library (here SAMRAI). In particular [core](core.md) will contain all numerical, physics, I/O etc.
- [amr](amr.md) : package that depends on the [core](core.md) package and on the AMR library SAMRAI. [AMR](amr.md) contains low-level code that specializes SAMRAI abstract clases (PatchData, RefineOperator, etc.) and code that provides a high-level interface to easily manipulate an abstract AMR hierarchy of plasma models (Transactions, MultiPhysicsIntegrator, etc.) and hides SAMRAI details (e.g. ResourcesManager).
- **tests** (phare/tests) :  groups all tests under the same directory hierarchy as the source code.

## Build

PHARE uses CMake to be built.


## Tests

PHARE uses GoogleTest for its test suite.


## Documentation

PHARE uses Doxygen to generate the developer documentation. It uses Sphinx to generate the user documentation in a separate repository called `pharead`.

- The developer documentation (this document)
- The user documentation (repository `pharead`)


## PHARESEE : PHARE visualization
- The visualization toolkit (repository `pharesee`)

## PAHREIN : python input API

- The python input kit (repository `pharein`)


## Dependencies

The main dependency of PHARE is the library SAMRAI. PHARE can build with SAMRAI in several fashions:

- If SAMRAI is pre-installed on the developer's machine and can be found by CMake, the CMake configuration will link this version
- or SAMRAI is not found by CMake then it is cloned as a subproject and compile with PHARE.
