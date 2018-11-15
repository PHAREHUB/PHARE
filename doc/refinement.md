# Refinement

The refinement consists in looping over the fine indexes and for each of them calculating the linear combination of the values of the two surrounding nodes of the same type (primal or dual). Note that the weights of the linear combination only depend on the position of the fine index relatively to the left and right surrounding coarse nodes. Therefore, these weights can be pre-computed, and the refinement just consists, for a given fine index, in finding which of the pre-computed weights to use, i.e. finding where is the fine index located in the coarse cell.

For a refinement ratio R, there are R couples of weights (leftWeight, rightWeight) to be computed, for the R fine indexes within the coarse cell. We store these pre-computed weights.


For primal quantities, fine indexes values are linearly interpolated between the two surrounding coarse primal nodes.
There is no difference between even and odd refinement ratios. 

For dual quantities, three coarse dual nodes are needed to get values on all dual fine nodes.


For the refinment, at least 1 ghost node is needed to perform the interpolation , which means that we will need ( max 1, nbrGhost for a given centering).


<img src="refinement.png" width="70%" height="70%"/>




