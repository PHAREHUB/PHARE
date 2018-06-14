# This scripts aims at generate some signal on two levels (cos, sin, other functions)
# and then to perform a coarsening on a given region


import numpy as np
import sys
import os

def sin(x,domainSize):
    return np.sin((2*np.pi/domainSize)*x);

def cos(x, domainSize):
    return np.cos((2*np.pi/domainSize)*x)

functionList = [cos,sin]

def test_coarsen_field_1d_yee_linear(path):

    nbrCellXCoarse = 40
    # dualAllocSizeCoarse = 42
    # primalAllocSizeCoarse = 43
    meshSizeCoarse = 0.2

    nbrCellXFine = 20
    # dualAllocSizeFine = 22
    # primalAllocSizeFine = 23
    meshSizeFine = 0.1


    domainSize = nbrCellXCoarse * meshSizeCoarse

    coarseXStartIndex = 5
    coarseXEndIndex = coarseXStartIndex + nbrCellXCoarse

    # in coarse index space, the fine path start at 9
    # the ratio is two, so the fine path start index
    # in fine index space is 18

    fineXStartIndex  =  18
    fineXEndIndex = fineXStartIndex + nbrCellXFine

    ghostWidth = 5
    ratio = 2

    #
    # for dual it is: extend the box of ghostWidth element in left and right
    # for primal it is : extend the box of ghostWidth element in left and right, and add 1 to the right
    #
    x_dual_fine = meshSizeFine * np.arange(fineXStartIndex-ghostWidth, fineXEndIndex+ghostWidth) + meshSizeFine/2.
    x_primal_fine = meshSizeFine * np.arange(fineXStartIndex-ghostWidth, fineXEndIndex+ghostWidth+1  )

    x_dual_coarse = meshSizeCoarse * np.arange(coarseXStartIndex-ghostWidth, coarseXEndIndex+ghostWidth  ) + meshSizeCoarse/2.
    x_primal_coarse = meshSizeCoarse * np.arange(coarseXStartIndex-ghostWidth, coarseXEndIndex+ghostWidth + 1  )


    isFirstFunction = True

    value_dual_fine_Data = []
    value_primal_fine_Data = []

    value_dual_coarse_Data = []
    value_primal_coarse_Data = []


    coarsed_dual_corse_data = []
    coarsed_primal_coarse_data = []

    for function in functionList:
        value_dual_fine = function(x_dual_fine,domainSize)
        value_primal_fine = function(x_primal_fine,domainSize)

        value_dual_coarse = function(x_dual_coarse,domainSize)
        value_primal_coarse = function(x_primal_coarse,domainSize)

        value_dual_after_coarse = np.copy(value_dual_coarse)
        value_primal_after_coarse = np.copy(value_primal_coarse)

        # let's set an arbitrary portion of the overlap
        coarseIndexDual = np.arange(10, 15 )
        coarseIndexPrimal = np.arange(10,16)


        for index in coarseIndexDual:
            # coarse index to fine index
            fineIndex = index*ratio

            fineIndex -= (fineXStartIndex - ghostWidth) # AMRToLocal

            coarseLocalIndex = index - (coarseXStartIndex - ghostWidth)

            value_dual_after_coarse[coarseLocalIndex] = value_dual_fine[fineIndex] * 0.5 \
                                                      + value_dual_fine[fineIndex+1] * 0.5

        for index in coarseIndexPrimal:
            # coarse index to fine index
            fineIndex = index*ratio

            fineIndex -= (fineXStartIndex - ghostWidth) # AMRToLocal
            coarseLocalIndex = index - (coarseXStartIndex - ghostWidth)

            value_primal_after_coarse[coarseLocalIndex] = value_primal_fine[fineIndex-1] * 0.25 \
                                                          + value_primal_fine[fineIndex] * 0.5 \
                                                          + value_primal_fine[ fineIndex+1 ] * 0.25


        if (isFirstFunction):
            value_dual_fine_Data = np.copy(value_dual_fine)
            value_primal_fine_Data = np.copy(value_primal_fine)

            value_dual_coarse_Data =  np.copy(value_dual_coarse)
            value_primal_coarse_Data =  np.copy(value_primal_coarse)

            coarsed_dual_corse_data = np.copy(value_dual_after_coarse)
            coarsed_primal_coarse_data = np.copy(value_primal_after_coarse)

            isFirstFunction = False

        else:
            value_dual_fine_Data = np.vstack((value_dual_fine_Data,value_dual_fine))
            value_primal_fine_Data = np.vstack((value_primal_fine_Data,value_primal_fine))

            value_dual_coarse_Data =  np.vstack((value_dual_coarse_Data,value_dual_coarse))
            value_primal_coarse_Data =  np.vstack((value_primal_coarse_Data,value_primal_coarse))

            coarsed_dual_corse_data = np.vstack((coarsed_dual_corse_data, value_dual_after_coarse))
            coarsed_primal_coarse_data = np.vstack((coarsed_primal_coarse_data, value_primal_after_coarse))


    filename_dual_fine = "dual_fine_original1d.txt"
    filename_primal_fine = "primal_fine_original1d.txt"

    filename_dual_coarse = "dual_coarse_original1d.txt"
    filename_primal_coarse = "primal_coarse_original1d.txt"

    filename_dual_coarse_coarsed = "dual_coarse_linear_coarsed_1d.txt"
    filename_primal_coarse_coarsed = "primal_coarse_linear_coarsed_1d.txt"


    np.savetxt(os.path.join(path, filename_dual_fine),value_dual_fine_Data, delimiter=" " )
    np.savetxt(os.path.join(path, filename_primal_fine),value_primal_fine_Data, delimiter=" " )

    np.savetxt(os.path.join(path, filename_dual_coarse),value_dual_coarse_Data, delimiter=" " )
    np.savetxt(os.path.join(path, filename_primal_coarse),value_primal_coarse_Data, delimiter=" " )

    np.savetxt(os.path.join(path, filename_dual_coarse_coarsed), coarsed_dual_corse_data, delimiter=" ");
    np.savetxt(os.path.join(path, filename_primal_coarse_coarsed), coarsed_primal_coarse_data, delimiter=" ")


def main(path="./"):

    if len(sys.argv) > 1:
        path=sys.argv[1]

    test_coarsen_field_1d_yee_linear(path)


if __name__ == '__main__':
    main()

