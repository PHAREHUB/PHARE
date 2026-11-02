# overview of the different concepts/parts of the tiling/gpu/threading branch,



## concepts now:

    PatchLevel
     Patch
      Cell
       Fields
       Particles

    Patch ghost values
      Cross Patch communications
      MPI synchronisations










<br /><br /><br /><br /><br /><br /><br /><br /><br /><br /><br /><br />
<br /><br /><br /><br /><br /><br /><br /><br /><br /><br /><br /><br />






















## concepts tiled:

    PatchLevel
     Patch
      Tile
       Cell
        Fields
        Particles

    Patch ghost values
      Cross Patch communications
      MPI synchronisations

    Tile ghost values
      Cross Tile communcations





<br /><br /><br /><br /><br /><br /><br /><br /><br /><br /><br /><br />
<br /><br /><br /><br /><br /><br /><br /><br /><br /><br /><br /><br />








































## discussing where we are


    cpu tiles
      1d init/advance tests
      3d updater/boris/sync
<br />


    gpu tiles
      3d updater/boris/sync
<br />


https://github.com/PhilipDeegan/PHARE/blob/d20bd65ea118d910110a6f616a4729a193977f5d/tests/functional/harris/harris_2d.py#L216

https://github.com/PhilipDeegan/PHARE/actions

<br /><br /><br /><br /><br /><br /><br /><br /><br /><br /><br /><br />
<br /><br /><br /><br /><br /><br /><br /><br /><br /><br /><br /><br />



































## what remains to be done

    Tile new divb fix
    2d+ validations (possible numeric issue)
    MPI validations (minor impl for pack/unpack)
    Missing GPU impls
    Optimizations
    Werrors




<br /><br /><br /><br /><br /><br /><br /><br /><br /><br /><br /><br />
<br /><br /><br /><br /><br /><br /><br /><br /><br /><br /><br /><br />












































##  what are the different parts (no details atm),

    ParticleArrays
      Arbitrary Storage and Layouts
      Allows for case by case Impls/Optimizations

https://github.com/PhilipDeegan/PHARE/blob/d20bd65ea118d910110a6f616a4729a193977f5d/src/core/data/particles/particle_array_def.hpp#L17-L27

https://github.com/PhilipDeegan/PHARE/blob/d20bd65ea118d910110a6f616a4729a193977f5d/src/core/def/phare_config.hpp#L59-L63

<br />

    TileSets
      Maps cells on tiles

    GridTileSet
      Maps fields to tiles

    ParticleArrayTileSet
      Maps ParticleArrays to tiles

<br />

    CPU Threading
      ThreadPools for a given task
      Multi-Stage threaded job executor for updater/boris
       Can also use threadpools

<br />

    GPU
      Multi-Stage threaded job executor for updater/boris
       Kernel per patch, block per tile, blocksize==warpsize, iterates

<br />

    Atomics
      Particle To Mesh
      Particle Leaving/Entering tiles

<br />

    Multi-Stage threaded job executor
      Implicit synchronisations between jobs
       Such that, the next job only runs when the previous is finished

       Why? Data Integrity.
        Example: We don't want to inject level ghost particles into
                  the domain while the domain is still being moved.

      Example steps:
        1 Move all pops for patch
        2 Synchronise across tiles
        3 Synchronise Domain/Ghosts
        4 Interpolate



<br /><br /><br /><br /><br /><br /><br /><br /><br /><br /><br /><br />
<br /><br /><br /><br /><br /><br /><br /><br /><br /><br /><br /><br />










































## their completion level and decide for a plan to discuss them in detail.

  `what remains to be done` - `where we are` == `completion level`

<br />

  plan to make a plan...


<br /><br /><br /><br /><br /><br /><br /><br /><br /><br /><br /><br />
<br /><br /><br /><br /><br /><br /><br /><br /><br /><br /><br /><br />

