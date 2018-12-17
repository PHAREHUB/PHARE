# Core package

The core library contains source code that is independant of the specific AMR library (here SAMRAI). Core groups:

- [data](core_data.md), containing code for data containers and spatial discretization (quantiy layout), more specifically:
	-  [electromag](electromag_package.md) defines the Electromag class
	-  [field](field_package.md) defines the Field container which stores gridded data independently from its layout on the grid
	-  [grid](grid_package.md) defines the GridLayout object, which codes how quantities are discretized on grids
	-  [ions](ions.md) defines all classes related with Ion particles, populations and moments
	-  [ndarray](ndarray.md) defines generic multidimensional array containers from which Field object inherit
	-  [particles](particles.md) defines classes related to particles and particle containers
	-  VecField defines a vector field
-  [numerics](numerics.md) groups all code related to solving equations numerically
	- Ampere
	- Faraday
	- Interpolator
	- Pusher
	- BoundaryCondition
-  [utilities](utilities.md) groups various utility functions
-  [hybrid](hybrid_package.md) defines constants and code related to hybrid model equations and quantitie
