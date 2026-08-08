#!/usr/bin/env python3

import allocSizes
import cellCenteredCoordinates
import fieldNodeCoordinates
import gridIndexing
import test_deriv
import test_laplacian
import test_linear_combinations_yee


def main():
    allocSizes.main("./")
    gridIndexing.main("./")
    cellCenteredCoordinates.main("./")
    fieldNodeCoordinates.main("./")
    test_deriv.main("./")
    test_laplacian.main("./")
    test_linear_combinations_yee.main("./")


if __name__ == "__main__":
    main()
