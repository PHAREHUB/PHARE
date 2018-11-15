<style>
	.markdown-body {
#		box-sizing: border-box;
#		min-width: 200px;
#		max-width: 980px;
#		margin: 0 auto;
#		padding: 45px;
	}

.section
{
 color:red;
}
		
		.subsection
		{
		text-decoration:underline;
		}
	}
</style>


<article class="markdown-body">



# Data communications in PHARE for hybrid hierarchies

Transaction of data between patches can occur when

- two neighboring patches exchange ghost node information
- a new refinement level is created in the hierarchy
- a level needs its coarse to fine boundary to be filled.
- level data is synchronized with the next coarser level


Level initialization takes place when: 

- a new refinement level is created, increasing the depth of the patch hierarchy. In this case, because the level did not exist before, all the data is obtained by refining the solution existing on the next coarser level.
- there is a *regridding*, i.e. an existing level is reshaped. Here, data is obtained from the old level where it has an intersection with the new level, and by refinement of the next coarser solution where it does not.





## <span class="section">Data to be communicated</span>

This section describes which data is communicated and when.

### General to any solver

- **model E** :
	- refined in initializeLevelData(), filled in whole level interior, patch ghosts and level borders
	- refined in `solver::advanceLevel()` when updating model, need to fill ghosts : `transaction.fillElectricGhosts()`
	- coarsened in `HybridTransaction::synchronize()`

- **model B** : 
	- in initializeLevelData(), filled in whole level interior, patch ghosts and level borders
	- in `solver::advanceLevel()` when updating model, need to fill ghosts : `transaction.fillMagneticGhosts()`
	- coarsened in HybridTransaction::synchronize()
	
- **model coarseToFineOld of the ith population**:
	- at initializeLevelData(). HybridTransaction fills these particles so to be able to use t=n and t=n+1 refined particles in coarse to fine boundaries and get moments on incomplete nodes close to level boundaries. Note that HybridTransaction need to fill `coarseToFineOld` with particles moved from `coarseToFineNew` at the end of last_step

- **model coarseToFineNew of the ith population** :
	- at first step in MultiPhysicsIntegrator::advanceLevel(), HybridTransaction, taking a model, fills `coarseToFineNew` particles with refined particles so that it can later compute moments from `coarseToFineNew` and `coarseToFineOld`

- **model coarseToFine of the ith population** :
	- each time coarseToFineOld is filled, it is copied into coarseToFine. No dedicated schedule needed.

- **Transaction Eold** : 
	- is used in communication in HybridTransaction::fillElectricGhosts() together with model E to be time interpolated at finer level ghost nodes.
	- copied from model E after solver::advanceLevel()

- **Transaction Bold** : 
	- is used in communication in HybridTransaction::fillElectricGhosts() together with model B to be time interpolated at finer level ghost nodes.
	- copied from model B after solver::advanceLevel()

- **model Ions total density and bulk velocity**:
	- coarsened to next coarser level ion total density and bulk velocity



### Specific to SolverPPC


- **SolverPPC Epred** :
	- in predictor 1 and 2 after ohm(). Needs E model n+1 and E transaction n on coarser

- **SolverPPC Bpred** : 
	- in predictor 1 and 2 after faraday(). Needs E model n+1 and E transaction n on coarser







## <span class="section"> Particle Communications</span>
### <span class="subsection">ParticlesData reminder</span>


Let's remind that a `ParticleData` has 5 `ParticleArray`

```
{
	ParticleArray* interior;
	ParticleArray* ghosts;
	ParticleArray* coarseToFine;
	ParticleArray* coarseToFineOld;
	ParticleArray* coarseToFineNew;
};
```

A ParticleData encapsulate all particle data required by a Patch.


See [ParticleData](particledata.md) for more details.

### <span class="subsection">Level Initialization</span>

#### Initializing a new level

This operation consists in creating a brand new finest level in the hierarchy. In SAMRAI, this operation occurs in `initializeLevelData()` when the `oldLevel` pointer is `nullptr`. At the end of this procedure, the interior of the level patches, their ghost region need to be filled, whether it is inside the level or at the coarse to fine boundary.

There are three important communications that need to happen:

1. Fill the interior of level patches from refined particles lying in the next coarser level
2. Fill the level border with refined particles lying in the next coarser level. Refined particles need to be put in the `coarseToFineOld` particle array since they are needed in `advanceLevel()` and are not available anymore on the coarser level at this point.
3. Fill the ghost regions of the level patches that overlap the interior of the level with clones of particles lying in the interior of neighbor patches.

The three steps are represented in the figure below

![](newLevelParticles.png)


##### FAQ: 

- Why don't we fill both interior and level borders with refined particles at the same time, we could do that using a PatchLevelBorderAndInteriorFillPattern?
	- The reason is that we need to use the particle refine operator to refine particles, and that operator needs to know in which `PaticleArray` of the `ParticlesData` it must put refined particles in. The way it knows that is with a template parameter so that it is compile time. Therefore a particle refine operator can only put particles in one `ParticleArray` and thus we need to do interior and borders separately.

- Why can't we fill the ghost regions of the level patches together with interior particles, i.e. from refining the next coarser level particles?
	- The reason is that particles in patch ghost regions in the level interior must be clones of particles in neighboring patches for the simulation to work. Filling these regions by refining particles from the next coarser level as we do for interior particles will only work if the splitting of coarse particles is deterministic. Since we don't know the splitting properties in this part of the code, we assume it is not deterministic and fill ghost particles by exchanging particles with neighbor level patches.



##### SAMRAI schedules and algorithms

1. Filling the interior of level patches can be done with the schedule 6 that takes the `level`, the `next_coarser_level` and the `hierarchy`, with the `PatchLevelInteriorFillPattern`. The particle refine operator will be templated with the Interior enum.
2. Filling the coarse to fine region can be done with the same schedule 6, but this time with `PatchLevelBorderFillPattern`. The particle refine operator will be templated with the `CoarseToFineOld` enum.
3. Filling the patch ghost regions can be done with the single level schedule, no need for a fill pattern.


Variables  | Schedule overload       | Space interpolator template arg | PatchLevelFillPattern
 :---:     | :---:                   | :---:                           | :---:
 Pop_i    | level ; ncl ; hierarchy  | Interior                        |  PatchLevelInteriorFillPattern
 Pop_i     | level ; ncl ; hierarchy | CoarseToFineOld                 |  PatchLevelBorderFillPattern
 Pop_i     | level                   | no interp                       |  none



#### Regridding

Regridding means that a level `i` is reshaped and the new level must be filled with particles. In this case, we want to fill the patches of the new level using particles from the old one it overlaps the new one, and refined particles from the coarser one elsewhere. This operation occurs in `initializeLevelData()` when `oldLevel` is not `nullptr`. At the end of this procedure, the interior of the level patches, their ghost region need to be filled, whether it is inside the level or at the coarse to fine boundary.


The steps are the same as for the intialization of a new level except now we can take from the old level instead of getting data exclusively from the next coarser.


This is described in the following figure:


![](regriddingParticles.png)



##### SAMRAI schedules and algorithms

1. Filling the interior of level patches can be done with the schedule 8 that takes the `dst_level`, the `src_level` the `next_coarser_level` and the `hierarchy`, with the `PatchLevelInteriorFillPattern`. The particle refine operator will be templated with the Interior enum.
2. Filling the coarse to fine region can be done with the same schedule 8, but this time with `PatchLevelBorderFillPattern`. The particle refine operator will be templated with the `CoarseToFineOld` enum.
3. Filling the patch ghost regions can be done with the single level schedule, no need for a fill pattern.


Variables  | prototype                                | Space interpolator template arg | PatchLevelFillPattern
 :---:     | :---:                                    | :---:                           | :---:
 Pop_i     | level                                    | no interp                       |  none
 Pop_i     | level dst ; level src ;  ncl ; hierarchy | Interior                        | PatchLevelInteriorFillPattern
 Pop_i     | level dst ; level src ;  ncl ; hierarchy | CoarseToFineOld                 | PatchLevelBorderFillPattern




### <span class = "subsection">Communications during advancement of the solution </span>

These are communications of particles that are needed to get ghost particles. 


#### Ghost particle region

These particles are needed to get complete density and flux estimates on the interior and ghost nodes of the moment grids. 


![](ghosts_width_particles.png)

The figure above shows primal nodes of 1D patches for interpolation orders 1, 2 and 3 around the patch boundary. The green nodes are nodes that receive density and fluxes from particles that are all in the patch interior. Blue nodes are domain patch nodes that lack contributions from particles that would exist outside the patch. Red nodes are ghost nodes.

For each order, the blue shape represents the shape factor of the particle that is epsilon from the patch boundary *inside* the patch. This particle contributes interp_order + 1 nodes around its position. At orders 2 and 3 these particles need one ghost node to project their contribution to moments. At order 1 it requires no ghost.

PHARE only needs (so far) the density to defined on patch nodes because the density is only involved in the calculation of the electric field which is dual in its component direction. The purple band represents the region in which one needs to have ghost particles. The boundaries of this region are defined by the patch boundary on one side and the position of the farthest particle contributing to the first node of the patch, on the other side. This definition depends on the interpolation order and is equal to (interp_order+1)/2.The light red shape represents the shape factor of the farthest particle that contributes to the first patch node (first blue node from the left).

The total number of ghost nodes must therefore be large enough to let all particles in the ghost region defined above project on the grid. In other words, there should be interp_order ghost cells.


InterpOrder | max loading distance (purple region) (fine dx unit)   | number of needed primal ghost nodes |
 :---:      | :---:                                                 | :---:                               |
 1          | 1                                                     | 1                                   |
 2          | (2+1)/2 = 1.5                                         | 2                                   |
 3          | (3+1)/2 = 2                                           | 3                                   |


This ghost region is the one to fill with "ghosts" and "refined particles".


#### Ghost communications

We need to get two kind of ghost particles:

- ghost particles that are copies of neighbor patch interior particles located in the patch ghost box
- refined particles that fall in the coarse to fine level boundaries

![](ghost_exchange.png)


The above figure represents two patches (solid line) ad their ghost zone (dashed lines). Two particles are represented. The filled and empty circles represent their position at t and t+dt, respectively. The red particle is moving close to the border of the black patch that overlaps the ghost box of the blue patch. The black particle is already in this region at time t, but is leaving the patch at time t+dt. 

Same level ghost communication consists in sending particles of the interior `ParticleArray` of a `ParticlesData` to a neighbor `ParticleData` ghost ParticleArray if they are also in that `ParticlesData` ghost box. For instance at time t, the communication woud copy the black particle of the black patch into the blue ghost particle array. At t+dt the communication would send the red particle of the black patch into the ghost particle array of the blue patch (while the black would have left it in the mean time, entering the interior of the patch).


**SAMRAI schedule involved**

Variables  | prototype                                | Space interpolator template arg | PatchLevelFillPattern
 :---:     | :---:                                    | :---:                           | :---:
 Pop_i     | level                                    | no interp                       |  none





Coarse to fine particle transactions are needed to fill the ghost zone of patches found at the level boundary (see the region between the dashed and solid lines in the above figure about regridding). 

At first_step, the solver needs to get particles in `coarseToFineNew` so to fill the moments on the level incomplete nodes. Note that at the moment the code is in advanceLevel(), `coarseToFineOld` ParticleArray is already filled, either from level initialization (see above) or from moving `coarseToFineNew` particles in it at the end of last time step.


**SAMRAI schedule involved**

Variables  | prototype                                | Space interpolator template arg | PatchLevelFillPattern
 :---:     | :---:                                    | :---:                           | :---:
 Pop_i     | level ; ncl ; hierarchy                  | Border                          |  PatchLevelBorderFillPattern


This operation is only needed at first step of the subcycle.
 
 
 


## <span class="section">Field Refinement</span>

### <span class="subsection">Electromagnetic field refinement</span>


#### Level Initialization



At level initialization, the model magnetic field always needs to be initialized, and everywhere, on interior, ghost and coarse to fine nodes. The electric field however may not be needed. If the solver starts with the Ohm's law (E = ...) then no need to initialize it. If however it starts with faraday or pushing particles, then it is needed.

The Transactions need to have their internal electric and magnetic field initialized. These fields represent the t=n state of the level that is used by finer level to get their ghost nodes at times between t=n et t=n+dt_coarse.

Solvers may need additional variables to be initialized. For instance, solvers may have internal copies of the electromagnetic field at past times.


None of the schedules need a `PatchLevelFillPattern`.


   Data Owner |  Variables | schedule prototype      | Optional ?| Purpose                                              
     :---:    | :---:      |     :---:               | :---:     | :---:                                                
  model       |  Ex,Ey,Ez  | level ; ncl ; hierarchy | yes       | have model E on grid                                 
  model       |  Bx,By,Bz  | level ; ncl ; hierarchy | no        | have model B on grid                                 
  Transaction |  Ex,Ey,Ez  | level ; ncl ; hierarchy | no        | hold E at t=n for time interpolation                 
  Transaction |  Bx,By,Bz  | level ; ncl ; hierarchy | no        | hold B at t=n for time interpolation                 
  Solver      |  Ex,Ey,Ez  | level ; ncl ; hierarchy | yes       | solver internal var. that needs init. ex: E at t=n-1 
  Solver      |  Bx,By,Bz  | level ; ncl ; hierarchy | yes       | solver internal var. that needs init. ex: B at t=n-1
     

#### Level Regridding

At regridding, essentially the same variables need to be communicated to the new level. The only difference is that now we need to use the schedule that takes the `dst_level`, the `src_level` and the `next_coarser_level`, to copy fields from oldLevel where it overlaps the new one.



### <span class="subsection">Moments refinement</span>

There is no need to refine moments in a hybrid code.


## <span class="section">Field Coarsering</span>

### <span class="subsection">Electromagnetic field coarsening</span>

### <span class="subsection">Moments coarsening</span>




</article>











