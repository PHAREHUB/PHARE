#!/usr/bin/env python3

import allocSizes
import gridIndexing
import cellCenteredCoordinates
import fieldNodeCoordinates
import deriv

def main():
    allocSizes.main('./')
    gridIndexing.main('./')
    cellCenteredCoordinates.main('./')
    fieldNodeCoordinates.main('./')
    deriv.main('./')



if __name__ == "__main__":
    main()
