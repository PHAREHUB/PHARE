# Field Coarsening

## definition
The coarsening of a field is the operation that is used to project the solution obtained on a fine level onto its parent (coarser).

We focus on the coarsening that preserves the centering.


## Refinement factor
The refinement factor (RF) between two levels is a vector of dimension 1, 2 or 3 of which each component is the refinement factor of the associated dimension. For instance RF = (2,4,5) means the mesh size is multiplied by 2, 4 and 5 when passing from the fine to the next coarser level, in the x, y and z dimensions, respectively. The RF can vary in the hierarchy, i.e. can be (2,3,5) between levels 2 and 3, but (2,1,2) between levels 3 and 4.


## The interpolation loop

The interpolation loop is the loop that calculates the value at a coarse index, from a linear combination of values on the fine grid.
In each dimension, we loop over the fine points.



```
for (i = 0; shift_i.size; ++i)
{
  for (j = 0; shift_j.size; ++j)
  {
    for (k=0; shift_k.size; ++k)
    {
      zInterp = weights_k[k] * f(i_start + shift_i[i],
                                       j_start + shift_j[j],
                                       k_start + shift_k[k]);      
    }      
    yInterp = zInterp * weights_j[j];      
  }
  coarseValue = yInterp * weights_i[u];

}

```

Doing so allows to decouple the problem of doing the interpolation, which depends on the dimensionality of the system, from the problem of finding the fine indexes and associated weights, which is a 1D problem.

We now focus on finding the fine indexes and their associated weights. Since this is the same problem in each direction, we thus reduce here to 1D considerations. The following figure summarizes all configurations we can find, in a given direction.


### Finding fine indexes

The goal here is to loop on the coarse mesh and find all find indexes  that contribute to the interpolation.


<img src="coarsening.png" width="70%" height="70%"/>



Two parameters are important:

- the parity of the RF
- the centering of the quantity

there are thus 4 distinct cases:

- Case A : even RF and primal index
- Case B : odd RF and primal index
- Case C : even RF and dual index
- Case D : odd RF and primal index


In **case A**, the coarse primal node falls on top of a fine primal node which AMR index is just RF\*coarse_index. There are RF fine cells in between coarse cell. Fine indexes coarse_index*RF + i, with i in [-RF/2, ..., 0, ..., RF/2] can contribute to the coarse_index. The fine nodes RF\*coarse_index -/+ RF/2 contribute to both RF(coarse_index+1) and RF\*(coarse_index-1) because RF being even, there is an odd number of fine primal nodes between each coarse primal node, and thus RF/2 falls in the middle. In this case, we need the coarse_index*RF fine node and up to the shared middle node.


In **case B**: the situation is very similar to case A, except here the RF being odd, there is an even number of primal indexes in between two coarse primal nodes. Therefore the fine node at coarse_index\*RF + RF/2 is not shared between coarse_index and coarse_index+1.


In **case C**, where we have an even RF and dual indexes, the coarse dual index falls on top of a fine primal index. There are therefore at minimum the two fine dual nodes around that fine primal nodes to consider. Note that here the fine node  coarse_index\*RF is the first fine dual node from the left within the coarse cell. There are RF fine dual nodes per coarse cell indexed 0, 1, 2, ... RF-1. We need the two nodes around dual coarse nodes and add nodes down to coarse_index*RF and up to coarse_index\*RF + RF-1.


In **case  D** the situation is very similar to case C, except not there is an odd number of fine cells within a coarse cell, i.e. an odd number of fine dual nodes in between two coarse primal nodes. The coarse dual node therefore falls onto a fine dual node that is exactly at coarse_index\*RF + RF/2. We need to take coarse_index\*RF+RF/2 and progressively add nodes on the left and right down to coarse_index\*RF and up to coarse_index\*RF + RF



### Weights

The weight calculation depends on the parity of the number of weights. For an even number of weights, we take :

- 1/x for coarse_index\*RF + RF/2 and coarse_index\*RF + RF/2-1
- 1/(i*x) for the i-th  point appart from those two points

Given that the sum of the weights must be equal to 1, we have x = 2\*sum(1/i) with i from 1 to NbrPoints/2. The coef 2 comes from the symmetry.


For an odd number of weights, we have:

- 1/x for the middle point coarse_index\*RF
- 1/((i+1)\*x) for all points around that middle point

Because the sum of the weights is equal to 1, we have :

x = 1 + 2\*sum(1/(i+1)) with i = 1..(n-1)/2

### Number of points necessary for correct coarsening operation.
In PHARE we will limit the refinment ratio to 10, so that we will only need at least 3 ghost node for primal centered data,
and 1 ghost for dual centered data. (by at least: max of nbrGhost with the given inteprOrder and 3(for primal, and 1 for dual))

